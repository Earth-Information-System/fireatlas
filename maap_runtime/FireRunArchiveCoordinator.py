import argparse
import s3fs
import glob
import dask.config
import geopandas as gpd
import datetime as dt

# @TODO add maap-py to the environment if needed, e.g.:
# pip install "git+https://github.com/MAAP-Project/maap-py.git@develop"        
# from maap.maap import MAAP

from functools import partial
from dask.distributed import Client
from fireatlas.FireMain import Fire_Forward
from fireatlas.FireRunDaskCoordinator import (
    job_data_update_checker,
    job_preprocess_region,
    job_preprocess_region_t, 
    get_timesteps_needing_region_t_processing
)
from fireatlas.FireConsts import YAML_ABS_PATH
from fireatlas.FireIO import copy_from_local_to_s3, s3_log_destination_path, s3_config_path
from fireatlas.FireTime import dt2t, t2dt, t_nb
from fireatlas.postprocess import allfires_filepath, allpixels_filepath, get_t_of_last_allfires_run
from fireatlas.utils import timed
from fireatlas import settings
from fireatlas.FireLog import logger

dask.config.set({'logging.distributed': 'error'})

# NOTE: this expects credentials to be resolvable globally
# via boto3/botocore common resolution paths
fs = s3fs.S3FileSystem(config_kwargs={"max_pool_connections": 10})

# @TODO where should region defintions be stored, and, should they be copied from s3 
# before settings initialization? 
# e.g. (from run_dps_cli.sh) copy_s3_object "s3://maap-ops-workspace/shared/gsfc_landslides/FEDSpreprocessed/${regnm}/.env" ../fireatlas/.env


def main(run_name):

    config_path = s3_config_path(run_name)
    if not fs.exists(config_path):
        raise FileNotFoundError(f"Run configuration file {config_path} does not exist on S3. "
                                "Please ensure the run configuration is uploaded to S3 before running this script.")


    # copy config file from s3 to local 
    fs.get(config_path, YAML_ABS_PATH)
    settings.__init__() # re-load settings with run_config.yaml
    logger.info(f"Finished loading settings from {config_path}")
    logger.info(settings.model_dump())

    # @TODO update this error message with correct location
    if not (settings.RUN_NAME & settings.TST & settings.TED):
        raise ValueError("Run parameters are not defined in run_config.yaml. "
        "To use this script, you must define the full run parameters and settings in " 
        " FEDSinput/region_definitions/run_config.yaml.")
    else: 
        # parse TST and TED 
        tst = settings.TST
        ted = settings.TED
        regnm = settings.RUN_NAME #@TODO set up the region definition correctly
        if settings.REGION_SHAPEFILE is not None: 
            region = [regnm, settings.REGION_SHAPEFILE]
        elif settings.REGION_BBOX is not None: 
            region = [regnm, settings.REGION_BBOX]
        else:
            raise ValueError("No region shape found. Did you set settings.REGION_SHAPEFILE or"
                             "settings.REGION_BBOX in run_config.yaml.")
    
    gpd.show_versions() # for debugging 

    logger.info(f"------------- Starting full run from {tst=} to {ted=} -------------")

    client = Client(n_workers=settings.N_DASK_WORKERS)
    logger.info(f"dask workers = {len(client.cluster.workers)}")

    # run the first two jobs in parallel
    data_update_futures = job_data_update_checker(client, tst, ted)
    region_future = client.submit(job_preprocess_region, region)
    
    # block until data update is complete
    client.gather(data_update_futures)
    
    data_upload_futures = client.map(
        partial(copy_from_local_to_s3, fs=fs),
        glob.glob(f"{settings.LOCAL_PATH}/{settings.PREPROCESSED_DIR}/*/*.txt")
    )
    # block until half-day timesteps and region are on s3
    timed(client.gather, text=f"Dask upload of {len(data_upload_futures) + 1} files")([*data_upload_futures, region_future])

    logger.info("------------- Done with preprocessing t -------------")

    # then run all region-plus-t in parallel that need it
    timesteps_needing_processing = get_timesteps_needing_region_t_processing(
        tst, ted, region, force=True
    )
    region_and_t_futures = client.map(
        partial(job_preprocess_region_t, region=region),
        timesteps_needing_processing
    )
    # block until preprocessing is complete
    client.gather(region_and_t_futures)
    
    logger.info("------------- Done with preprocessing region + t -------------")
    
    # run fire forward for the next settings.ARCHIVE_RUN_JOB_SIZE days 
    # produce and upload only allfires and allpixels 

    # get t of latest allfires/allpixels if any 
    t_saved = get_t_of_last_allfires_run(tst, ted, region=region, location=settings.READ_LOCATION)
    if t_saved is not None: 
        run_tst = t_nb(t_saved)
    else: 
        run_tst = tst 
    run_ted = min(ted, t2dt(run_tst) + dt.timedelta(days=settings.ARCHIVE_RUN_JOB_SIZE)) 
    run_ted = dt2t(run_ted)
    
    logger.info(f"------------- Running Fire_Forward for {run_tst=} to {run_ted=} -------------")

    try: 
        Fire_Forward(tst=run_tst, ted=run_ted, region=region, restart=False)
        copy_from_local_to_s3(allpixels_filepath(run_tst, run_ted, region, location="local"), fs=fs)
        copy_from_local_to_s3(allfires_filepath(run_tst, run_ted, region, location="local"), fs=fs)
    except KeyError as e: 
        logger.warning(f"Fire_Forward has already run. {e}")

    logger.info(f"------------- Done running Fire_Forward for {run_tst=} to {run_ted=} -------------")
    
    if run_ted < ted:

        print(f"*************** Mock submitting next job for {t_nb(run_ted)} to {ted} ****************")
        # @TODO actually submit next job
        logger.info("------------- Submitted next job to DPS -------------")
    else:
        logger.info("------------- Full run completed -------------")

    fs.put_file(settings.LOG_FILEPATH, s3_log_destination_path(region[0]))

    return  

if __name__ == "__main__": 
    """
    CLI script for coordinating retrospective "archive" runs of 
    FEDS on the MAAP DPS platform. 

    If the user-specified time range is greater than 
    settings.ARCHIVE_RUN_JOB_SIZE days, the run will be split into multiple 
    jobs. Each job reads the latest checkpoint data, does any neccessary preprocessing, 
    then runs Fire_Forward on the region for settings.ARCHIVE_RUN_JOB_SIZE days.
    Upon completetion, if we still haven't reached the end of the time range, 
    this script will submit another job for the next chunk of time. 

    Parameters
    ----------
    run_name : str 
        Name of the region to run FEDS over, e.g. "ArchiveCONUS"

    """

    parser = argparse.ArgumentParser()
    parser.add_argument("run_name", type=str, help="Name of the run definition to execute")
    args = parser.parse_args()

    main(args.run_id)

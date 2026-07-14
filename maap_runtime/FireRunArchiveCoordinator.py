import argparse
import os
import subprocess
import s3fs
import glob
import dask.config
import geopandas as gpd
import datetime as dt
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
from fireatlas.FireIO import (
    copy_from_local_to_s3, 
    copy_from_local_to_veda_s3, 
    s3_log_destination_path, s3_config_path, 
    s3_metadata_destination_path
)
from fireatlas.FireTime import dt2t, t2dt, t_nb
from fireatlas.postprocess import (
    all_dir, 
    allfires_filepath, 
    allpixels_filepath, 
    combined_lf_perims_nifc_join, 
    find_largefires, 
    get_t_of_last_allfires_run, 
    read_allfires_gdf, 
    read_allpixels, 
    save_large_fires_layers, 
    save_large_fires_nplist, 
    save_snapshots
)
from fireatlas.utils import timed
from fireatlas import settings
from fireatlas.FireLog import logger, write_run_metadata
from maap.maap import MAAP

dask.config.set({'logging.distributed': 'error'})

# NOTE: this expects credentials to be resolvable globally
# via boto3/botocore common resolution paths
fs = s3fs.S3FileSystem(config_kwargs={"max_pool_connections": 10})

def main(run_name, copy_to_veda=False):

    wallclock_start = dt.datetime.now()

    if os.path.exists(YAML_ABS_PATH): 
        logger.info(f"run_config.yaml file found at {YAML_ABS_PATH}. Including settings overrides.")
    else: 
        logger.info(f"run_config.yaml NOT found at {YAML_ABS_PATH}.")
    
    logger.info(settings.model_dump())

    if settings.RUN_NAME is None or settings.TST is None:
        raise ValueError("Run parameters are not defined in run_config.yaml. "
        "To use this script, you must define the full run parameters and settings in " 
        " FEDSinput/run_definitions/{run_name}/run_config.yaml.")
    else: 
        # parse TST and TED 
        tst = settings.TST

        if settings.TED is not None: 
            ted = settings.TED
        else: 
            # if no end time set, use current time (for NRT runs) 
            ctime = dt.datetime.now(tz=dt.timezone.utc)
            if ctime.hour >= 18: 
                ampm = 'PM' 
            else: 
                ampm = 'AM'
            ted = [ctime.year, ctime.month, ctime.day, ampm]
        
        
        regnm = settings.RUN_NAME 
        if settings.REGION_SHAPEFILE is not None: 
            region = [regnm, settings.REGION_SHAPEFILE]
        elif settings.REGION_BBOX is not None: 
            region = [regnm, settings.REGION_BBOX]
        else:
            raise ValueError("No region shape found. Did you set settings.REGION_SHAPEFILE or"
                             "settings.REGION_BBOX in run_config.yaml?")
    
    gpd.show_versions() # for debugging 

    # log commit hash of current fireatlas version
    try:
        logger.info("fireatlas current branch: " + subprocess.check_output(["git", "rev-parse", "--abbrev-ref", "HEAD"],stderr=subprocess.DEVNULL, text=True).strip())
        logger.info("commit hash of fireatlas version used: " + subprocess.check_output(["git", "rev-parse", "HEAD"],stderr=subprocess.DEVNULL, text=True).strip())
    except: 
        pass 

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

    # Pick up where we left off if possible; calculate size of next run chunk 
    if t_saved is None or t2dt(t_saved) <= t2dt(tst) or t2dt(t_saved) > t2dt(ted): 
        run_tst = tst 
    else: 
        run_tst = t_saved 
    run_ted = dt2t(min(t2dt(ted), t2dt(run_tst) + dt.timedelta(days=settings.ARCHIVE_RUN_JOB_SIZE)))
    
    logger.info(f"------------- Running Fire_Forward for {run_tst=} to {run_ted=} -------------")

    try: 
        allfires, allpixels, t_saved = Fire_Forward(tst=run_tst, ted=run_ted, region=region, restart=False)
        allfires_gdf = allfires.gdf
        copy_from_local_to_s3(allpixels_filepath(run_tst, run_ted, region, location="local"), fs=fs)
        copy_from_local_to_s3(allfires_filepath(run_tst, run_ted, region, location="local"), fs=fs)
    except KeyError as e: 
        logger.warning(f"Fire_Forward has already run. {e}")
        allpixels = read_allpixels(tst, ted, region)
        allfires_gdf = read_allfires_gdf(tst, ted, region)
        t_saved = run_ted

    logger.info(f"------------- Done running Fire_Forward for {run_tst=} to {run_ted=} -------------")
    
    if t2dt(run_ted) < t2dt(ted):

        logger.info(f"------------- Submitting next job for {t_nb(run_ted)} to {ted} -------------")
        
        maap = MAAP(maap_host='api.maap-project.org')
        job = maap.submitJob(
            identifier=f"job-eis-feds-archive:1.5.2",
            algo_id="eis-feds-archive",
            version="1.5.2",
            username="zbecker", 
            queue="maap-dps-eis-worker-128gb",
            run_id=run_name
        )

        logger.info(f"------------- Submitted next job to DPS. Submission status: {job.status} -------------")

    else:
        # all done with run: do postprocessing 

        if t_saved is not None: 
            snapshot_tst = t_saved
        else: 
            snapshot_tst = run_tst
        snapshot_futures = save_snapshots(allfires_gdf, region, snapshot_tst, ted, client=client)
        large_fires = find_largefires(allfires_gdf)
        save_large_fires_nplist(allpixels, region, large_fires, tst)
        save_large_fires_layers(allfires_gdf, region, large_fires, tst, ted, client=client)

        client.gather(snapshot_futures)

        # If flag matching flat set, add overlaps with this year's NIFC incidents to 
        # CombinedLargefire/lf_perimeter.fgb for ted only. 
        if settings.DO_NIFC_MATCHING:
            logger.info("Started NIFC matching")
            combined_lf_perims_nifc_join(tst, ted, region, active_only=True, time_filter=None)
            logger.info("Finished NIFC matching")


        # take all fire forward output and upload all outputs in parallel
        data_dir = all_dir(tst, region, location="local")
        fgb_s3_upload_futures = client.map(
            partial(copy_from_local_to_s3, fs=fs),
            glob.glob(os.path.join(data_dir, "*", "*", "*.fgb"))
        )
        # block until everything is uploaded
        timed(client.gather, text=f"Dask upload of {len(fgb_s3_upload_futures)} files")(fgb_s3_upload_futures)

        if copy_to_veda:
            # take latest fire forward output and upload to VEDA S3 in parallel
            fgb_veda_upload_futures = client.map(
                partial(copy_from_local_to_veda_s3, fs=fs, regnm=region[0]),
                glob.glob(os.path.join(data_dir, "*", f"{ted[0]}{ted[1]:02}{ted[2]:02}{ted[3]}", "*.fgb"))
            )
            timed(client.gather, text=f"Dask upload of {len(fgb_veda_upload_futures)} files")(fgb_veda_upload_futures)

        logger.info("------------- Full run completed -------------")

    wallclock_end = dt.datetime.now() 
    run_duration = wallclock_end - wallclock_start 
    logger.info(f"This job completed in {str(run_duration)}")
    
    client.close()
    
    # write environment metadata and copy to s3
    metadata_path = write_run_metadata()
    fs.put_file(metadata_path, s3_metadata_destination_path(region[0])) 
    logger.info(f"Copied environment information to {s3_metadata_destination_path(region[0])}")
    
    # finally, copy log file to s3
    fs.put_file(settings.LOG_FILEPATH, s3_log_destination_path(region[0], run_ted))

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

    main(args.run_name)

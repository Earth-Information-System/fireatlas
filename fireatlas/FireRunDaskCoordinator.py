import json
import argparse
import os
import glob
from itertools import chain

import s3fs

from typing import Tuple
from dask.distributed import Client
from dask import delayed
from dask.delayed import Delayed
from datetime import datetime

from fireatlas.FireMain import Fire_Forward
from fireatlas.FireTypes import Region, TimeStep
from fireatlas.utils import timed
from fireatlas.postprocess import (
    save_allpixels,
    save_allfires_gdf,
    save_snapshots,
    find_largefires,
    save_large_fires_layers,
    save_large_fires_nplist,
)
from fireatlas.preprocess import preprocess_region_t, preprocess_region
from fireatlas.DataCheckUpdate import update_VNP14IMGTDL, update_VJ114IMGTDL
from fireatlas.FireIO import copy_from_local_to_s3
from fireatlas.FireTime import t_generator
from fireatlas.FireLog import logger
from fireatlas import settings

# NOTE: the current eis queue 64gb is set up as an
# AWS r5.2xlarge instance with 64gb RAM and 8 CPUs
MAX_WORKERS = 6

# NOTE: this expects credentials to be resolvable globally
# via boto3/botocore common resolution paths
fs = s3fs.S3FileSystem(config_kwargs={"max_pool_connections": 10})


def validate_json(s):
    try:
        return json.loads(s)
    except ValueError:
        raise argparse.ArgumentTypeError("Not a valid JSON string")



def concurrent_copy_from_local_to_s3(eventual_results: Tuple[Delayed], local_filepath: str):
    copy_from_local_to_s3(local_filepath, fs)


def job_fire_forward(eventual_results: Tuple[Delayed], region: Region, tst: TimeStep, ted: TimeStep):
    logger.info(f"Running code for {region[0]} from {tst} to {ted} with source {settings.FIRE_SOURCE}")

    allfires, allpixels = Fire_Forward(tst=tst, ted=ted, region=region, restart=False)
    allpixels_filepath = save_allpixels(allpixels, tst, ted, region)
    allfires_filepath = save_allfires_gdf(allfires.gdf, tst, ted, region)

    copy_from_local_to_s3(allpixels_filepath, fs)
    copy_from_local_to_s3(allfires_filepath, fs)

    save_snapshots(allfires.gdf, region, tst, ted)

    large_fires = find_largefires(allfires.gdf)
    save_large_fires_nplist(allpixels, region, large_fires, tst)
    save_large_fires_layers(allfires.gdf, region, large_fires, tst)


def job_preprocess_region_t(
        eventual_results: Tuple[Delayed],
        region: Region, t: TimeStep):
    logger.info(f"Running preprocessing code for {region[0]} at {t=} with source {settings.FIRE_SOURCE}")
    filepath = preprocess_region_t(t, region=region)
    copy_from_local_to_s3(filepath, fs)


def job_preprocess_region(region: Region):
    filepath = preprocess_region(region)
    copy_from_local_to_s3(filepath, fs)


def job_data_update_checker():
    """"""
    try:
        # Download SUOMI-NPP
        update_VNP14IMGTDL()
        # Download NOAA-20
        update_VJ114IMGTDL()
    except Exception as exc:
        logger.exception(exc)

@timed
def Run(region: Region, tst: TimeStep, ted: TimeStep):

    ctime = datetime.now()
    if tst in (None, "", []):  # if no start is given, run from beginning of year
        tst = [ctime.year, 1, 1, 'AM']

    if ted in (None, "", []):  # if no end time is given, set it as the most recent time
        if ctime.hour >= 18:
            ampm = 'PM'
        else:
            ampm = 'AM'
        ted = [ctime.year, ctime.month, ctime.day, ampm]

    list_of_timesteps = list(t_generator(tst, ted))
    dask_client = Client(n_workers=MAX_WORKERS)
    logger.info(f"dask workers = {len(dask_client.cluster.workers)}")

    # run the first two jobs in parallel
    data_input_results = delayed(job_data_update_checker)()
    region_results = delayed(job_preprocess_region)(region)

    # block on the first two jobs, then uploads raw satellite files from `job_data_update_checker` in parallel
    data_upload_results = [
        delayed(concurrent_copy_from_local_to_s3)([data_input_results, region_results], local_filepath)
        for local_filepath in glob.glob(f"{settings.LOCAL_PATH}/{settings.PREPROCESSED_DIR}/*/*.txt")
    ]

    # blocks on data upload results, then runs all region-plus-t in parallel
    region_and_t_results = [
        delayed(job_preprocess_region_t)(data_upload_results, region, t)
        for t in list_of_timesteps
    ]

    # blocks on region_and_t_results and then runs fire forward algorithm (which cannot be run in parallel)
    fire_forward_results = delayed(job_fire_forward)(region_and_t_results, region, tst, ted)

    # blocks on fire forward algorithm run and then uploads all snapshots/largefire outputs in parallel
    data_dir = os.path.join(settings.LOCAL_PATH, settings.OUTPUT_DIR, region[0], str(tst[0]))
    fgb_upload_results = [
        delayed(concurrent_copy_from_local_to_s3)([fire_forward_results,], local_filepath)
        for local_filepath in list(chain(
            glob.glob(os.path.join(data_dir, "Snapshot", "*", "*.fgb")),
            glob.glob(os.path.join(data_dir, "Largefire", "*", "*.fgb"))
        ))
    ]

    # take the DAG and run it
    dag = delayed(lambda x: x)(fgb_upload_results)
    dag.compute()
    dask_client.close()


if __name__ == "__main__":
    """coordinating all jobs
    
    Example:
    python3 FireRunDaskCoordinator.py --regnm="CONUS" --bbox="[-126,24,-61,49]" --tst="[2023,6,1,\"AM\"]" --ted="[2023,9,1,\"AM\"]"
    """
    
    parser = argparse.ArgumentParser()
    parser = argparse.ArgumentParser()
    parser.add_argument("--regnm", type=str)
    parser.add_argument("--bbox", type=validate_json)
    parser.add_argument("--tst", type=validate_json)
    parser.add_argument("--ted", type=validate_json)
    args = parser.parse_args()
    Run([args.regnm, args.bbox], args.tst, args.ted)

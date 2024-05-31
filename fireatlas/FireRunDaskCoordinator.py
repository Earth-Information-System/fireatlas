import json
import argparse
import os
import glob
from functools import partial

import s3fs

from dask.distributed import Client
from datetime import datetime

from fireatlas.FireMain import Fire_Forward
from fireatlas.FireTypes import Region, TimeStep
from fireatlas.utils import timed
from fireatlas.postprocess import (
    all_dir,
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


def concurrent_copy_from_local_to_s3(local_filepath: str):
    copy_from_local_to_s3(local_filepath, fs)


def job_fire_forward(region: Region, tst: TimeStep, ted: TimeStep):
    logger.info(f"Running FireForward code for {region[0]} from {tst} to {ted} with source {settings.FIRE_SOURCE}")

    allfires, allpixels = Fire_Forward(tst=tst, ted=ted, region=region, restart=False)
    
    allpixels_filepath = allpixels_filepath(tst, ted, region, location="local")
    allfires_filepath = allfires_filepath(tst, ted, region, location="local")
    
    copy_from_local_to_s3(allpixels_filepath, fs)
    copy_from_local_to_s3(allfires_filepath, fs)

    save_snapshots(allfires.gdf, region, tst, ted)

    large_fires = find_largefires(allfires.gdf)
    save_large_fires_nplist(allpixels, region, large_fires, tst)
    save_large_fires_layers(allfires.gdf, region, large_fires, tst)


def job_preprocess_region_t(region: Region, t: TimeStep):
    logger.info(f"Running preprocess-region-t code for {region[0]} at {t=} with source {settings.FIRE_SOURCE}")
    filepath = preprocess_region_t(t, region=region)
    copy_from_local_to_s3(filepath, fs)


def job_preprocess_region(region: Region):
    logger.info(f"Running preprocess-region JSON for {region[0]}")
    filepath = preprocess_region(region)
    copy_from_local_to_s3(filepath, fs)


def job_data_update_checker():
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
    client = Client(n_workers=MAX_WORKERS)
    logger.info(f"dask workers = {len(client.cluster.workers)}")

    # run the first two jobs in parallel
    data_update_futures = client.submit(job_data_update_checker)
    region_futures = client.submit(job_preprocess_region, region)
    
    # block until data update is complete
    client.gather(data_update_futures)

    # uploads raw satellite files from `job_data_update_checker` in parallel
    data_upload_futures = client.map(
        concurrent_copy_from_local_to_s3,
        glob.glob(f"{settings.LOCAL_PATH}/{settings.PREPROCESSED_DIR}/*/*.txt")
    )
    # block until half-day timesteps and region are on s3
    client.gather([*data_upload_futures, *region_futures])

    # then run all region-plus-t in parallel
    region_and_t_futures = client.map(
        partial(job_preprocess_region_t, region=region),
        list_of_timesteps
    )
    # block until preprocessing is complete
    client.gather(region_and_t_futures)

    # run fire forward algorithm (which cannot be run in parallel)
    job_fire_forward(region=region, tst=tst, ted=ted)

    # take all fire forward output and upload all snapshots/largefire outputs in parallel
    fgb_upload_futures = client.map(
        concurrent_copy_from_local_to_s3,
        glob.glob(os.path.join(all_dir, "*", "*", "*.fgb"))
    )
    # block until everything is uploaded
    client.gather(fgb_upload_futures)
    client.close()


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

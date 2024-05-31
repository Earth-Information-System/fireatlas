import json
import argparse
import os
import glob
from functools import partial

import s3fs

from dask.distributed import Client
from datetime import datetime, date, timezone

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
from fireatlas.preprocess import (
    check_preprocessed_file,
    preprocessed_filename,
    preprocess_input_file,
    preprocess_region_t,
    preprocess_region,
)

from fireatlas.DataCheckUpdate import update_VNP14IMGTDL, update_VJ114IMGTDL
from fireatlas.FireIO import copy_from_local_to_s3, copy_from_local_to_veda_s3, VNP14IMGML_filepath, VJ114IMGML_filepath
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


def get_timesteps_needing_region_t_processing(
    tst: TimeStep,
    ted: TimeStep,
    region: Region,
    sat = None,
):
    needs_processing = []
    for t in t_generator(tst, ted):
        filepath = preprocessed_filename(t, sat=sat, region=region)
        if not settings.fs.exists(filepath):
            needs_processing.append(t)
    return needs_processing


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
    save_large_fires_layers(allfires.gdf, region, large_fires, tst, ted)


def job_preprocess_region_t(t: TimeStep, region: Region):
    logger.info(f"Running preprocess-region-t code for {region[0]} at {t=} with source {settings.FIRE_SOURCE}")
    filepath = preprocess_region_t(t, region=region)
    copy_from_local_to_s3(filepath, fs)


def job_preprocess_region(region: Region):
    logger.info(f"Running preprocess-region JSON for {region[0]}")
    filepath = preprocess_region(region)
    copy_from_local_to_s3(filepath, fs)


def job_data_update_checker(client: Client, tst: TimeStep, ted: TimeStep):
    source = settings.FIRE_SOURCE

    futures = []
    if source == "VIIRS":
        sats = ["SNPP", "NOAA20"]
    else:
        sats = [source]
    
    for sat in sats:
        if sat == "SNPP":
            monthly_filepath_func = VNP14IMGML_filepath
            NRT_update_func = update_VNP14IMGTDL
        if sat == "NOAA20":
            monthly_filepath_func = VJ114IMGML_filepath
            NRT_update_func = update_VJ114IMGTDL

        # first check if there are any monthly files that need preprocessing
        timesteps = check_preprocessed_file(tst, ted, sat=sat, freq="NRT")
        monthly_timesteps = list(set([(t[0], t[1]) for t in timesteps]))
        monthly_filepaths = [monthly_filepath_func(t) for t in monthly_timesteps]
        
        # narrow down to the monthly filepaths and timesteps that actually exist
        indices = [i for i, f in enumerate(monthly_filepaths) if f is not None]
        monthly_filepaths = [monthly_filepaths[i] for i in indices]
        monthly_timesteps = [monthly_timesteps[i] for i in indices]

        # set up monthly jobs
        futures.extend(client.map(preprocess_input_file, monthly_filepaths))

        # calculate any remaining missing dates
        missing_dates = [date(*t) for t in timesteps if (t[0], t[1]) not in monthly_timesteps]
        logger.info(f"looking for {[str(d) for d in missing_dates]} on {sat}")
        
        # set up NRT jobs
        futures.extend(client.map(NRT_update_func, missing_dates))

    return futures

@timed
def Run(region: Region, tst: TimeStep, ted: TimeStep):

    ctime = datetime.now(tz=timezone.utc)
    if tst in (None, "", []):  # if no start is given, run from beginning of year
        tst = [ctime.year, 1, 1, 'AM']

    if ted in (None, "", []):  # if no end time is given, set it as the most recent time
        if ctime.hour >= 18:
            ampm = 'PM'
        else:
            ampm = 'AM'
        ted = [ctime.year, ctime.month, ctime.day, ampm]

    client = Client(n_workers=MAX_WORKERS)
    logger.info(f"dask workers = {len(client.cluster.workers)}")
 
    # run the first two jobs in parallel
    data_update_futures = job_data_update_checker(client, tst, ted)
    region_future = client.submit(job_preprocess_region, region)
    
    # block until data update is complete
    client.gather(data_update_futures)

    # uploads raw satellite files from `job_data_update_checker` in parallel
    data_upload_futures = client.map(
        partial(copy_from_local_to_s3, fs=fs),
        glob.glob(f"{settings.LOCAL_PATH}/{settings.PREPROCESSED_DIR}/*/*.txt")
    )
    # block until half-day timesteps and region are on s3
    client.gather([*data_upload_futures, region_future])

    # then run all region-plus-t in parallel that need it
    timesteps_needing_processing = get_timesteps_needing_region_t_processing(
        tst, ted, region
    )
    region_and_t_futures = client.map(
        partial(job_preprocess_region_t, region=region),
        timesteps_needing_processing
    )
    # block until preprocessing is complete
    client.gather(region_and_t_futures)

    # run fire forward algorithm (which cannot be run in parallel)
    job_fire_forward(region=region, tst=tst, ted=ted)

    # take all fire forward output and upload all outputs in parallel
    fgb_s3_upload_futures = client.map(
        partial(copy_from_local_to_s3, fs=fs),
        glob.glob(os.path.join(all_dir, "*", "*", "*.fgb"))
    )

    # take latest fire forward output and upload to VEDA S3 in parallel
    fgb_veda_upload_futures = client.map(
        partial(copy_from_local_to_veda_s3, fs=fs, regnm=region[0]),
        glob.glob(os.path.join(all_dir, "*", f"{ted[0]}{ted[1]:02}{ted[2]:02}{ted[3]}", "*.fgb"))
    )
    # block until everything is uploaded
    client.gather([*fgb_s3_upload_futures, *fgb_veda_upload_futures])
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

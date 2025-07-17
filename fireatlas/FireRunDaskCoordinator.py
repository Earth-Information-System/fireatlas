import json
import argparse
import os
import glob
import fsspec
import datetime as dt
from functools import partial

import s3fs

import dask.config
from dask.distributed import Client
from datetime import datetime, date, timezone, timedelta

from fireatlas.FireMain import Fire_Forward
from fireatlas.FireTypes import Region, TimeStep
from fireatlas.utils import timed
from fireatlas.postprocess import (
    all_dir,
    allfires_filepath,
    allpixels_filepath,
    save_snapshots,
    find_largefires,
    save_large_fires_layers,
    save_large_fires_nplist,
    read_allfires_gdf,
    read_allpixels,
    combined_lf_perims_nifc_join
)
from fireatlas.preprocess import (
    check_preprocessed_file,
    preprocessed_filename,
    preprocess_input_file,
    preprocess_region_t,
    preprocess_region,
    preprocessed_region_filename,
)

from fireatlas.DataCheckUpdate import update_FIRMS, update_VJ114IMGTDL, update_VNP14IMGTDL, get_FIRMS_data_availability
from fireatlas.FireIO import (
    copy_from_local_to_s3, 
    copy_from_local_to_veda_s3, 
    VNP14IMGML_filepath, 
    VJ114IMGML_filepath, 
    VJ114IMGTDL_filepath, 
    VNP14IMGTDL_filepath, 
    FIRMS_VIIRS_SNPP_SP_filepath, 
    FIRMS_VIIRS_SNPP_NRT_filepath, 
    FIRMS_VIIRS_NOAA20_SP_filepath, 
    FIRMS_VIIRS_NOAA20_NRT_filepath, 
    FIRMS_VIIRS_NOAA21_NRT_filepath
)
from fireatlas.FireTime import t_generator, t_nb, t2dt, dt2t, d2t
from fireatlas.FireLog import logger
from fireatlas import settings

dask.config.set({'logging.distributed': 'error'})


# NOTE: this expects credentials to be resolvable globally
# via boto3/botocore common resolution paths
fs = s3fs.S3FileSystem(config_kwargs={"max_pool_connections": 10})

logger.info(settings.model_dump())

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
    force = False,
):
    needs_processing = []
    for t in t_generator(tst, ted):
        filepath = preprocessed_filename(t, sat=sat, region=region)
        if not settings.fs.exists(filepath):
            needs_processing.append(t)

    if force:
        # for NRT make sure the current day window
        # is constantly being refreshed to incorporate batch updates
        # and ignore if a extra timesteps end up duplicated in processing
        dtst = date(*ted[:-1]) - timedelta(days=1)
        dted = date(*ted[:-1])
        if dtst < dted:
            needs_processing.extend([t for t in t_generator(d2t(dtst.year, dtst.month, dtst.day, 'AM'), ted)
                                     if t not in needs_processing])
    return needs_processing


def job_fire_forward(client: Client, region: Region, tst: TimeStep, ted: TimeStep):
    logger.info(f"Running FireForward code for {region[0]} from {tst} to {ted} with source {settings.FIRE_SOURCE}")

    try:
        allfires, allpixels, t_saved = Fire_Forward(tst=tst, ted=ted, region=region, restart=False)
        copy_from_local_to_s3(allpixels_filepath(tst, ted, region, location="local"), fs)
        copy_from_local_to_s3(allfires_filepath(tst, ted, region, location="local"), fs)
        allfires_gdf = allfires.gdf
        if t_saved is None:
            # NOTE: this happens if we're running a region full-on
            # from start to finish that has never been run before
            # and therefore no existin allpixels/allfires save has been found
            t_saved = tst
    except KeyError as e:
        logger.warning(f"Fire forward has already run. {e}")
        allpixels = read_allpixels(tst, ted, region)
        allfires_gdf = read_allfires_gdf(tst, ted, region)
        # NOTE: this means we've already found an
        # allfires and allpixels save for this ted timestep
        t_saved = ted

    snapshot_futures = save_snapshots(allfires_gdf, region, t_saved, ted, client=client)

    large_fires = find_largefires(allfires_gdf)
    save_large_fires_nplist(allpixels, region, large_fires, tst)
    save_large_fires_layers(allfires_gdf, region, large_fires, tst, ted, client=client)
    
    client.gather(snapshot_futures)


def job_preprocess_region_t(t: TimeStep, region: Region):
    logger.info(f"Running preprocess-region-t code for {region[0]} at {t=} with source {settings.FIRE_SOURCE}")
    filepath = preprocess_region_t(t, region=region)
    copy_from_local_to_s3(filepath, fs)


def job_preprocess_region(region: Region):
    output_filepath = preprocessed_region_filename(region)
    if settings.fs.exists(output_filepath):
        logger.info(f"Preprocessed region is already on {settings.READ_LOCATION}.")
        return
    
    logger.info(f"Running preprocess-region JSON for {region[0]}")
    filepath = preprocess_region(region)
    copy_from_local_to_s3(filepath, fs)


def job_nrt_current_day_updates(client: Client):
    """hourly update the NRT files and prep
    Updates files for today and the two previous days. 
    """
    futures, source, now = [], settings.FIRE_SOURCE, datetime.now()

    if source == "VIIRS":
        sats = ["SNPP", "NOAA20", "NOAA21"]
    else: 
        sats=[source]

    for sat in sats:
        futures.extend(client.map(update_FIRMS, *[
            (now.date(), (now-timedelta(days=1)).date(), (now-timedelta(days=2)).date()), 
            (sat, sat, sat),
            ("NRT", "NRT", "NRT")
        ]))

    return futures

def job_data_update_checker(client: Client, tst: TimeStep, ted: TimeStep):
    """
    Checks to see if any input data within the time range needs to be downloaded or preprocessed. 

    Prefers SP (standard product) data over NRT data. If SP is available, it will be used
    instead of the NRT data for the same day. 

    Returns: 
    --------
    futures : list[dask.distributed.client.Future]
        List of download and preprocessing jobs to execute

    """
    source = settings.FIRE_SOURCE
    location = settings.READ_LOCATION

    fs = fsspec.filesystem(location, use_listings_cache=False)

    futures = []
    if source == "VIIRS": 
        sats = ["SNPP", "NOAA20", "NOAA21"]
    else: 
        sats = [source]

    for sat in sats: 

        sp_start, sp_end, nrt_start, nrt_end = get_FIRMS_data_availability(sat)

        if sat == "SNPP":
            monthly_filepath_func = VNP14IMGML_filepath
            sp_filepath_func = FIRMS_VIIRS_SNPP_SP_filepath
            nrt_filepath_func = FIRMS_VIIRS_SNPP_NRT_filepath
        elif sat == "NOAA20":
            monthly_filepath_func = VJ114IMGML_filepath
            sp_filepath_func = FIRMS_VIIRS_NOAA20_SP_filepath
            nrt_filepath_func = FIRMS_VIIRS_NOAA20_NRT_filepath
        elif sat == "NOAA21":
            sp_filepath_func = None 
            nrt_filepath_func = FIRMS_VIIRS_NOAA21_NRT_filepath

        # gives list of timesteps for which there is no preprocessed file available
        timesteps = check_preprocessed_file(tst, ted, sat=sat, freq="NRT")

        if len(timesteps) < 1: # no processing needed
            return futures
        
        # there are no monthly arachive files for NOAA21 yet, so only check for SNPP and NOAA20
        if sat in ["SNPP", "NOAA20"]:
            monthly_timesteps = list(set([(t[0], t[1]) for t in timesteps]))
            monthly_filepaths = [monthly_filepath_func(t) for t in monthly_timesteps]

            # narrow down to the monthly filepaths and timesteps that actually exist 
            indices = [i for i, f in enumerate(monthly_filepaths) if f is not None]
            monthly_timesteps = [monthly_timesteps[i] for i in indices]
            monthly_filepaths = [monthly_filepaths[i] for i in indices]

            # set up jobs to preprocess existing monthly files that need it
            futures.extend(client.map(preprocess_input_file, monthly_filepaths))

            # calculate any remaining missing timesteps not in monthly files
            timesteps = [
                t for t in timesteps if (t[0], t[1]) not in monthly_timesteps
            ]

        for t in timesteps:
            d = dt.datetime(t[0], t[1], t[2])

            if sp_start and (sp_start <= d <= sp_end):  
                # SP is now available for this day; do we already have it? 
                # If not, download and overwrite old preprocessed file.
                fp = sp_filepath_func(t)
                if fs.exists(fp):

                    pfp = preprocessed_filename(t, sat, location=location)
                    if not fs.exists(pfp):
                        futures.append(client.submit(preprocess_input_file, fp))
                else:
                    futures.append(client.submit(update_FIRMS, d, sat, "SP"))
            elif nrt_start <= d <= nrt_end:
                fp = nrt_filepath_func(t)
                if fs.exists(fp):
                    pfp = preprocessed_filename(t, sat, location=location)
                    if not fs.exists(pfp):
                        futures.append(client.submit(preprocess_input_file, fp)) 
                else: 
                    futures.append(client.submit(update_FIRMS, d, sat, "NRT"))
            else:
                # Neither NRT nor SP available
                logger.warning(f"No data available for {sat} on {t[0]}-{t[1]}-{t[2]}.")

    return futures

@timed
def Run_local(region: Region, tst: TimeStep, ted: TimeStep, copy_to_veda: bool=False):
    """
    Coordinates all parts of a run: region preprocessing, downloading and preprocessing
    input fire detection data if needed, running FireForward, and saving snapshot layers.
    Similar to Run, but does not attempt to read from or write to s3 at all. Like Run, 
    uses a Dask client to parallelize some computations on the local machine, making 
    use of multiple CPU cores when available. 
    Remember to set settings.READ_LOCATION to "local"!
    """

    client = Client(n_workers=settings.N_DASK_WORKERS)
    region_future = client.submit(preprocess_region, region)
    logger.info(f"Running preprocess-region JSON for {region[0]}")
    data_update_futures = job_data_update_checker(client, tst, ted)

    client.gather(data_update_futures)
    client.gather(region_future)

    logger.info("------------- Done with preprocessing t -------------")

    # then run all region-plus-t in parallel that need it
    timesteps_needing_processing = get_timesteps_needing_region_t_processing(
        tst, ted, region, force=True
    )
    region_and_t_futures = client.map(
        partial(preprocess_region_t, region=region, force=True),
        timesteps_needing_processing
    )
    # block until preprocessing is complete
    client.gather(region_and_t_futures)

    logger.info("------------- Done with preprocessing region + t -------------")

    # run fire forward algorithm (which cannot be run in parallel)

    logger.info(f"Running FireForward code for {region[0]} from {tst} to {ted} with source {settings.FIRE_SOURCE}")

    try:
        allfires, allpixels, t_saved = Fire_Forward(tst=tst, ted=ted, region=region, restart=False)
        allfires_gdf = allfires.gdf
        if t_saved is None:
            # NOTE: this happens if we're running a region full-on
            # from start to finish that has never been run before
            # and therefore no existin allpixels/allfires save has been found
            t_saved = tst
    except KeyError as e:
        logger.warning(f"Fire forward has already run. {e}")
        allpixels = read_allpixels(tst, ted, region)
        allfires_gdf = read_allfires_gdf(tst, ted, region)
        # NOTE: this means we've already found an
        # allfires and allpixels save for this ted timestep
        t_saved = ted

    snapshot_futures = save_snapshots(allfires_gdf, region, t_saved, ted, client=client)

    large_fires = find_largefires(allfires_gdf)
    save_large_fires_nplist(allpixels, region, large_fires, tst)
    save_large_fires_layers(allfires_gdf, region, large_fires, tst, ted, client=client)

    client.gather(snapshot_futures)

    # If flag matching flat set, add overlaps with this year's NIFC incidents to 
    # CombinedLargefire/lf_perimeter.fgb for ted only. 
    if settings.DO_NIFC_MATCHING:
        logger.info("Started NIFC matching")
        combined_lf_perims_nifc_join(
            tst, 
            ted, 
            region, 
            active_only=settings.NIFC_MATCHING_ACTIVE_ONLY, 
            time_filter=None
        )
        logger.info("Finished NIFC matching")

    logger.info("------------- Done -------------")

    client.close()

@timed
def Run(region: Region, tst: TimeStep, ted: TimeStep, copy_to_veda: bool):

    ctime = datetime.now(tz=timezone.utc)
    if tst in (None, "", []):  # if no start is given, run from beginning of year
        tst = [ctime.year, 1, 1, 'AM']

    if ted in (None, "", []):  # if no end time is given, set it as the most recent time
        if ctime.hour >= 18:
            ampm = 'PM'
        else:
            ampm = 'AM'
        ted = [ctime.year, ctime.month, ctime.day, ampm]
    
    logger.info(f"------------- Starting full run from {tst=} to {ted=} -------------")

    client = Client(n_workers=settings.N_DASK_WORKERS)
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
    
    # run fire forward algorithm (which cannot be run in parallel)
    job_fire_forward(region=region, tst=tst, ted=ted, client=client)

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

    logger.info("------------- Done -------------")

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
    parser.add_argument('--no-veda-copy', dest='copy_to_veda', action='store_false', default=True,
                        help="defaults to True but if passed will stop a copy to VEDA s3 bucket")
    args = parser.parse_args()
    Run([args.regnm, args.bbox], args.tst, args.ted, args.copy_to_veda)

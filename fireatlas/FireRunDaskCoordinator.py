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
    preprocess_region_t,
    preprocess_region,
    preprocessed_region_filename,
    preprocess_monthly_file,
    preprocess_daily_file
)

from fireatlas.DataCheckUpdate import update_FIRMS, get_FIRMS_data_availability
from fireatlas.FireIO import (
    copy_from_local_to_s3,
    copy_from_local_to_veda_s3,
    VNP14IMGML_filepath,
    VJ114IMGML_filepath,
    FIRMS_VIIRS_SNPP_SP_filepath,
    FIRMS_VIIRS_SNPP_NRT_filepath,
    FIRMS_VIIRS_NOAA20_SP_filepath,
    FIRMS_VIIRS_NOAA20_NRT_filepath,
    FIRMS_VIIRS_NOAA21_NRT_filepath
)
from fireatlas.FireTime import t_generator, d2t, t_nm, t_nd
from fireatlas.FireLog import logger
from fireatlas import settings
import geopandas as gpd

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
        sats = [source]

    for sat in sats:
        futures.extend(client.map(update_FIRMS, *[
            (now.date(), (now-timedelta(days=1)).date(), (now-timedelta(days=2)).date()),
            (sat, sat, sat),
            ("NRT", "NRT", "NRT")
        ]))

    return futures


def job_data_update_checker(client: Client, tst: TimeStep, ted: TimeStep):
    """
    Checks to see if any input data within the time range needs to be preprocessed.

    If settings.FIRE_NRT == False, only tries to preprocess already existing monthly
    input files (VNP14IMGML and VJ114IMGML).
    If settings.FIRE_NRT == True, tries to download any missing NRT input data from FIRMS,
    then preprocess any unprocessed NRT data (FIRMS_VIIRS_SNPP_NRT, FIRMS_VIIRS_NOAA20_NRT,
    FIRMS_VIIRS_NOAA21_NRT). Does not try to use monthly files.

    NOTE: If settings.FIRE_NRT and any input files are needed,
    blocks for downloads inside this function and returns only preprocessing futures.

    NOTE: Does not automatically reprocess a timestep that was previously preprocessed
    from NRT data when the standard data product becomes available.

    Generally, assumes that if a preprocessed file for a date already exists, it does
    not need to be reprocessed unless it is from the most recent two days of the
    NRT record.

    Returns:
    --------
    futures : list[dask.distributed.client.Future]
        List of preprocessing jobs to execute
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
        if not settings.FIRE_NRT:

            # look for already-downloaded monthly files
            if sat == "SNPP":
                monthly_filepath_func = VNP14IMGML_filepath
            elif sat == "NOAA20":
                monthly_filepath_func = VJ114IMGML_filepath
            elif sat == "NOAA21":
                logger.warning("No standard products available for NOAA21. "
                               "Did you mean to set FireConsts.FIRE_NRT = True?")
                continue

            # gives list of timesteps for which there is no preprocessed file available
            timesteps = check_preprocessed_file(tst, ted, sat=sat, freq="monthly")

            if len(timesteps) < 1:  # no processing needed for this sat
                continue

            monthly_filepaths = [monthly_filepath_func(t) for t in timesteps]

            # we don't need to preprocess outside of time range, but we do need
            # the previous and next months to preprocess the first and last days
            prev_month, next_month = t_nm(tst, "previous"), t_nm(ted, "next")
            for m in [prev_month, next_month]:
                if not monthly_filepath_func(m):
                    logger.warning(f"No monthly input file found for {m} for {sat}")

            indices = [i for i, f in enumerate(monthly_filepaths) if f is not None]
            missing_indices = [i for i, f in enumerate(monthly_filepaths) if f is None]

            for i in missing_indices:
                logger.warning(f"No monthly input file found for {timesteps[i]} for {sat}")

            existing_timesteps = [timesteps[i] for i in indices]

            futures.extend(client.map(partial(preprocess_monthly_file, sat=sat), existing_timesteps))

        elif settings.FIRE_NRT:

            if sat == "SNPP":
                nrt_filepath_func = FIRMS_VIIRS_SNPP_NRT_filepath
                sp_filepath_func = FIRMS_VIIRS_SNPP_SP_filepath
            elif sat == "NOAA20":
                nrt_filepath_func = FIRMS_VIIRS_NOAA20_NRT_filepath
                sp_filepath_func = FIRMS_VIIRS_NOAA20_SP_filepath
            elif sat == "NOAA21":
                nrt_filepath_func = FIRMS_VIIRS_NOAA21_NRT_filepath
                sp_filepath_func = None

            # gives list of timesteps for which there is no preprocessed file available
            timesteps = check_preprocessed_file(tst, ted, sat=sat, freq="NRT")

            if len(timesteps) < 1:  # no processing needed for this sat
                continue

            # check firms data availability
            sp_start, sp_end, nrt_start, nrt_end = get_FIRMS_data_availability(sat)

            # use these to ensure all downloads are done before any preprocessing starts
            download_futures = {}  # (t, sat) -> dask future
            preprocess_tasks = {}  # (t, sat) -> filepath

            for t in timesteps:
                d = dt.datetime(t[0], t[1], t[2])

                if d > nrt_end:
                    logger.warning(f"No data available for {sat} on {t[0]}-{t[1]}-{t[2]}: date out of range.")
                    continue
                elif d >= nrt_start:  # in NRT availability range
                    fp = nrt_filepath_func(t)

                    if fs.exists(fp):
                        preprocess_tasks[(t, sat)] = fp
                    # if we don't already have this input file, try to download from FIRMS
                    else:
                        download_futures[(t, sat)] = client.submit(update_FIRMS, d, sat, "NRT")
                elif sp_start and d >= sp_start:  # check if sp_start because NOAA21 does not have yet
                    # in standard product availability range
                    fp = sp_filepath_func(t)
                    if fs.exists(fp):
                        preprocess_tasks[(t, sat)] = fp
                    else:
                        download_futures[(t, sat)] = client.submit(update_FIRMS, d, sat, "SP")
                else:
                    # either before sp_start, or this is NOAA21 (so, no sp_start) and it is before
                    # nrt start. either way, warn but allow
                    logger.warning(f"No data available for {sat} on {t[0]}-{t[1]}-{t[2]}. "
                                   "Date may be out of range.")

            # need to have these available to preprocess tst and ted, if possible
            prev_day = t_nd(tst, "previous")
            next_day = t_nd(ted, "next")

            for t in [prev_day, next_day]:
                d = dt.datetime(t[0], t[1], t[2])

                if d > nrt_end:
                    logger.warning(f"No data available for {sat} on {t[0]}-{t[1]}-{t[2]}. Date out of range.")
                elif d >= nrt_start:
                    fp = nrt_filepath_func(t)
                    if not fs.exists(fp):
                        update_FIRMS(d, sat, "NRT")
                elif d >= sp_start:
                    fp = sp_filepath_func(t)
                    if not fs.exists(fp):
                        update_FIRMS(d, sat, "SP")
                else:
                    logger.warning(f"No data available for {sat} on {t[0]}-{t[1]}-{t[2]}. "
                                   "Date may be out of range.")

            if len(download_futures) > 0:
                # block to finish downloads before starting any preprocessing
                downloaded_paths = client.gather(download_futures)
                preprocess_tasks.update(downloaded_paths)

            # schedule preprocessing
            for (tk, satk), fp in preprocess_tasks.items():
                tk = list(tk)
                futures.append(client.submit(preprocess_daily_file, fp, tk, satk))

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

    gpd.show_versions()

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

    CLI script for coordinating full FEDS runs.

    This script is the main entry point for orchestrating DPS runs
    and can also be used locally.

    Parameters
    ----------
    --regnm : str
        Name of the region to run FEDS over, e.g. "CONUS", "Central_Asia"
        If using a predefined shapefile or looking to pick up settings from
        a .env file in FEDSpreprocessed, this must match the name of the
        region shapefile and folder exactly.

    --bbox : str (JSON list)
        Rectangular lat/lon bounding box. FEDS will only run over this region.
        Only active fire detections within this region will be used.
        If reg_shp is provided, it will override bbox and this can be left empty.
        If the region has already been run before and has an existing geometry
        in FEDSpreprocessed/REGION/REGION.json, that will override the bbox.


    --tst : str (JSON list)
        Time start. FEDS will start running at this timestep.
        Provided as a JSON list- be careful to escape the quotes around AM/PM.
        If tst is "" or "[]", the first day of the current year will be used.
        Example: "[2023,6,1,\"AM\"]"

    --ted : str (JSON list)
        Time end. FEDS will run through this timestep (inclusive).
        If no ted is given, the most recent timestep will be used.
        Leave empty ("" or "[]") for NRT runs.


    --no-veda-copy : flag, optional
        If set, disables copying output to the VEDA S3 bucket to be ingested
        and served via OGC API. Use for all testing and local runs.


    Usage example:
    python3 FireRunDaskCoordinator.py --regnm="example_CONUS" \\
        --bbox="[-126,24,-61,49]" \\
        --reg_shp="" \\
        --tst="[2023,6,1,\\"AM\\"]" \\
        --ted="[2023,9,1,\\"AM\\"]" \\
        --no-veda-copy
    """

    parser = argparse.ArgumentParser()
    parser.add_argument("--regnm", type=str)
    parser.add_argument("--bbox", type=validate_json)
    parser.add_argument("--tst", type=validate_json)
    parser.add_argument("--ted", type=validate_json)
    parser.add_argument('--no-veda-copy', dest='copy_to_veda', action='store_false', default=True,
                        help="defaults to True but if passed will stop a copy to VEDA s3 bucket")
    args = parser.parse_args()

    Run([args.regnm, args.bbox], args.tst, args.ted, args.copy_to_veda)

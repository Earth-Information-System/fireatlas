import json
import argparse
import os
import glob
from typing import Tuple
from dask.distributed import Client
from dask import delayed
from dask.delayed import Delayed

from FireLog import logger
from FireTypes import Region, TimeStep
from utils import timed
import datetime

MAX_WORKERS = 3


def validate_json(s):
    try:
        return json.loads(s)
    except ValueError:
        raise argparse.ArgumentTypeError("Not a valid JSON string")


def job_fire_forward(eventual_results: Tuple[Delayed], region: Region, tst: TimeStep, ted: TimeStep):
    import FireIO, FireConsts, FireMain, postprocess
    from FireLog import logger

    ctime = datetime.now()

    if tst in (None, ""):  # if no start is given, run from beginning of year
        tst = [ctime.year, 1, 1, 'AM']

    if ted in (None, ""):  # if no end time is given, set it as the most recent time
        if ctime.hour >= 18:
            ampm = 'PM'
        else:
            ampm = 'AM'
        ted = [ctime.year, ctime.month, ctime.day, ampm]

    logger.info(f"Running code for {region[0]} from {tst} to {ted} with source {FireConsts.firesrc}")

    allfires, allpixels = FireMain.Fire_Forward(tst=tst, ted=ted, restart=False, region=region)
    allpixels_filepath = postprocess.save_allpixels(allpixels, tst, ted, region)
    allfires_filepath = postprocess.save_allfires_gdf(allfires.gdf, tst, ted, region)

    FireIO.copy_from_local_to_s3(allpixels_filepath)
    FireIO.copy_from_local_to_s3(allfires_filepath)

    postprocess.save_snapshots(allfires.gdf, region, tst, ted)

    large_fires = postprocess.find_largefires(allfires.gdf)
    postprocess.save_large_fires_nplist(allpixels, region, large_fires, tst)
    postprocess.save_large_fires_layers(allfires.gdf, region, large_fires, tst)

    for filepath in glob.glob(
            os.path.join(FireConsts.get_diroutdata(location="local"), region[0], str(tst[0]), "Snapshot", "*",
                         "*.fgb")):
        FireIO.copy_from_local_to_s3(filepath)

    for filepath in glob.glob(
            os.path.join(FireConsts.get_diroutdata(location="local"), region[0], str(tst[0]), "Largefire", "*",
                         "*.fgb")):
        FireIO.copy_from_local_to_s3(filepath)


def job_preprocess_region_t(
        eventual_result1: Delayed,
        eventual_result2: Delayed, region: Region, t: TimeStep):
    import FireIO, FireConsts, preprocess
    from FireLog import logger

    logger.info(f"Running preprocessing code for {region[0]} at {t=} with source {FireConsts.firesrc}")
    output_filepath = preprocess.preprocess_region_t(t, sensor=FireConsts.firesrc, region=region)
    FireIO.copy_from_local_to_s3(output_filepath)


def job_preprocess_region(region: Region):
    import FireIO, preprocess

    filepath = preprocess.preprocess_region(region)
    FireIO.copy_from_local_to_s3(filepath)


def job_data_update_checker():
    """"""
    import DataCheckUpdate, FireConsts, FireIO

    try:
        # Download SUOMI-NPP
        DataCheckUpdate.update_VNP14IMGTDL()
        # Download NOAA-20
        DataCheckUpdate.update_VJ114IMGTDL()
    except Exception as exc:
        logger.exception(exc)
    finally:
        for filepath in glob.glob(os.path.join(FireConsts.get_dirprpdata(location="local"), "*", "*.txt")):
            FireIO.copy_from_local_to_s3(filepath)


@timed
def Run(region: Region, tst: TimeStep, ted: TimeStep):
    import FireTime

    list_of_timesteps = list(FireTime.t_generator(tst, ted))
    dask_client = Client(n_workers=MAX_WORKERS)
    logger.info(f"dask workers = {len(dask_client.cluster.workers)}")

    data_input_results = delayed(job_data_update_checker)()
    region_results = delayed(job_preprocess_region)(region)
    region_and_t_results = [
        delayed(job_preprocess_region_t)(data_input_results, region_results, region, t)
        for t in list_of_timesteps
    ]
    fire_forward_results = delayed(job_fire_forward)(region_and_t_results, region, tst, ted)
    fire_forward_results.compute()
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

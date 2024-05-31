import json
import argparse
import fireatlas
from datetime import datetime, UTC
from functools import partial

from dask.distributed import Client
from fireatlas.utils import timed
from fireatlas import FireRunDaskCoordinator
from fireatlas import settings
from fireatlas.FireLog import logger


def validate_json(s):
    try:
        return json.loads(s)
    except ValueError:
        raise argparse.ArgumentTypeError("Not a valid JSON string")


@timed
def Run(region: fireatlas.FireTypes.Region, tst: fireatlas.FireTypes.TimeStep, ted: fireatlas.FireTypes.TimeStep):
    logger.info(f"Running preprocess region-at-t code for {region[0]} at {tst=} with source {settings.FIRE_SOURCE}")
    
    ctime = datetime.now(tz=UTC)
    if tst in (None, "", []):  # if no start is given, run from beginning of year
        tst = [ctime.year, 1, 1, 'AM']

    if ted in (None, "", []):  # if no end time is given, set it as the most recent time
        if ctime.hour >= 18:
            ampm = 'PM'
        else:
            ampm = 'AM'
        ted = [ctime.year, ctime.month, ctime.day, ampm]

    client = Client(n_workers=FireRunDaskCoordinator.MAX_WORKERS)
    logger.info(f"dask workers = {len(client.cluster.workers)}")

    # then run all region-plus-t in parallel that need it
    timesteps_needing_processing = FireRunDaskCoordinator.get_timesteps_needing_region_t_processing(
        tst, ted, region
    )
    region_and_t_futures = client.map(
        partial(FireRunDaskCoordinator.job_preprocess_region_t, region=region),
        timesteps_needing_processing
    )
    client.gather(region_and_t_futures)
    client.close()


if __name__ == "__main__":
    """ The main code to run preprocessing for a region and time period. It writes to a dedicated directory on s3.

    Example:
    python3 FireRunByRegionAndT.py --regnm="WesternUS" --tst="[2023,6,1,\"AM\"]"
    """

    parser = argparse.ArgumentParser()
    parser.add_argument("--regnm", type=str)
    parser.add_argument("--tst", type=validate_json)
    parser.add_argument("--ted", type=validate_json)
    args = parser.parse_args()
    Run([args.regnm, None], args.tst, args.ted)

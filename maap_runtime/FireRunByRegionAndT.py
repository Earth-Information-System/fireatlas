import json
import argparse
import fireatlas
from functools import partial

from dask.distributed import Client
from fireatlas.utils import timed
from fireatlas import FireRunDaskCoordinator
from fireatlas import settings
from fireatlas.FireLog import logger
from fireatlas.FireTime import t_generator


def validate_json(s):
    try:
        return json.loads(s)
    except ValueError:
        raise argparse.ArgumentTypeError("Not a valid JSON string")


@timed
def Run(region: fireatlas.FireTypes.Region, tst: fireatlas.FireTypes.TimeStep, ted: fireatlas.FireTypes.TimeStep):
    logger.info(f"Running preprocess region-at-t code for {region[0]} at {tst=} with source {settings.FIRE_SOURCE}")
    list_of_timesteps = list(t_generator(tst, ted))
    client = Client(n_workers=FireRunDaskCoordinator.MAX_WORKERS)
    logger.info(f"dask workers = {len(client.cluster.workers)}")
    region_and_t_futures = client.map(
        partial(FireRunDaskCoordinator.job_preprocess_region_t, region=region),
        list_of_timesteps
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

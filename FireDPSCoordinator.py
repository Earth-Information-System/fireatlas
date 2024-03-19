import json
import argparse
import os

from FireLog import logger
import FireConsts
from utils import timed
from FireTypes import TimeStep, Region
import preprocess


def validate_json(s):
    try:
        return json.loads(s)
    except ValueError:
        raise argparse.ArgumentTypeError("Not a valid JSON string")


@timed
def Run(region: Region, t: TimeStep):
    import FireIO

    # check to see if region has been processed.
    # maybe process region

    # for each t between tst and ted
        # for each satellite specified
            # check to see if half day files have been created
            # maybe kick off job to create

        # check to see if the the region_t files have been created
        # maybe kick off job to create


if __name__ == "__main__":
    """ The main code to run preprocessing for a region and time period. It writes to a dedicated directory on s3.
    
    Example:
    python3 FireRunByRegionAndT.py --regnm="WesternUS" --t="[2023,6,1,\"AM\"]"
    """
    
    parser = argparse.ArgumentParser()
    parser = argparse.ArgumentParser()
    parser.add_argument("--regnm", type=str)
    parser.add_argument("--bbox", type=validate_json)

    parser.add_argument("--tst", type=validate_json)
    parser.add_argument("--ted", type=validate_json)
    args = parser.parse_args()

    try:
        Run([args.regnm, None], args.tst)
    except Exception as e:
        logger.exception(e)

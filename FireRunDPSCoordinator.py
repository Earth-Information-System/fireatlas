import json
import argparse

from FireLog import logger
from utils import timed
from FireTypes import TimeStep, Region


def validate_json(s):
    try:
        return json.loads(s)
    except ValueError:
        raise argparse.ArgumentTypeError("Not a valid JSON string")


@timed
def Run(region: Region, tst: TimeStep, ted: TimeStep):
    import FireConsts, FireIO, preprocess

    # check to see if region has been processed.
    # maybe process region

    # for each t between tst and ted
        # for each satellite specified
            # check to see if half day files have been created
            # maybe kick off job to create

        # check to see if the the region_t files have been created
        # maybe kick off job to create

    # run fire forward which goes for whole year and writes all the outputs

if __name__ == "__main__":
    """ The main code to run preprocessing for a region and time period. It writes to a dedicated directory on s3.
    
    Example:
    python3 FireRunDPSCoordinator.py --regnm="WesternUS" --t="[2023,6,1,\"AM\"]"
    """
    
    parser = argparse.ArgumentParser()
    parser = argparse.ArgumentParser()
    parser.add_argument("--regnm", type=str)
    parser.add_argument("--bbox", type=validate_json)

    parser.add_argument("--tst", type=validate_json)
    parser.add_argument("--ted", type=validate_json)
    args = parser.parse_args()

    try:
        Run([args.regnm, args.bbox], args.tst, args.ted)
    except Exception as e:
        logger.exception(e)

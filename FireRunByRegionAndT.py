import json
import argparse
import os

from utils import timed
from FireTypes import TimeStep, Region


def validate_json(s):
    try:
        return json.loads(s)
    except ValueError:
        raise argparse.ArgumentTypeError("Not a valid JSON string")


@timed
def RegionAndTRun(region: Region, t: TimeStep):
    import FireIO, FireConsts, preprocess
    from FireLog import logger

    logger.info(f"Running preprocessing code for {region[0]} at {t=} with source {FireConsts.firesrc}")

    output_filepath = preprocess.preprocess_region_t(t, sensor=FireConsts.firesrc, region=region)
    FireIO.copy_from_local_to_s3(output_filepath)


if __name__ == "__main__":
    """ The main code to run preprocessing for a region and time period. It writes to a dedicated directory on s3.
    
    Example:
    python3 FireRunByRegionAndT.py --regnm="WesternUS" --t="[2023,6,1,\"AM\"]"
    """
    
    parser = argparse.ArgumentParser()
    parser.add_argument("--regnm", type=str)
    parser.add_argument("--t", type=validate_json)
    args = parser.parse_args()
    try:
        RegionAndTRun([args.regnm, None], args.t)
    except Exception as e:
        from FireLog import logger
        logger.exception(e)

import json
import argparse
import os

from FireLog import logger
from utils import timed
from FireTypes import Region


def validate_json(s):
    try:
        return json.loads(s)
    except ValueError:
        raise argparse.ArgumentTypeError("Not a valid JSON string")


@timed
def Run(region: Region):
    import FireIO, preprocess

    filepath = preprocess.preprocess_region(region)

    FireIO.copy_from_local_to_s3(filepath)


if __name__ == "__main__":
    """ The main code to run preprocessing for a region and time period. It writes to a dedicated directory on s3.
    
    Example:
    python3 FireRunByRegionAndT.py --regnm="WesternUS" --t="[2023,6,1,\"AM\"]"
    """
    
    parser = argparse.ArgumentParser()
    parser.add_argument("--regnm", type=str)
    parser.add_argument("--bbox", type=validate_json)
    args = parser.parse_args()

    try:
        Run([args.regnm, args.bbox])
    except Exception as e:
        logger.exception(e)

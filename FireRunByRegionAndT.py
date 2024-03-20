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
def RegionAndTRun(region: Region, t: TimeStep):
    import FireIO
    
    print(f"Running preprocessing code for {region[0]} at {t=} with source {FireConsts.firesrc}")

    output_filepaths = preprocess.preprocess_region_t(t, FireConsts.firesrc, region=region)
    for filepath in output_filepaths:
        dst = os.path.join(FireConsts.get_dirprpdata(location="s3", strip_protocol=True), *filepath.split("/")[-3:])

        FireIO.copy_from_local_to_s3(filepath, dst)


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
        logger.exception(e)

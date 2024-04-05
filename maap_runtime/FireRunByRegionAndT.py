import json
import argparse

import fireatlas
from fireatlas.utils import timed
from fireatlas import FireRunDaskCoordinator
from fireatlas import settings
from fireatlas import FireLog


def validate_json(s):
    try:
        return json.loads(s)
    except ValueError:
        raise argparse.ArgumentTypeError("Not a valid JSON string")


@timed
def Run(region: fireatlas.FireType.Region, tst: fireatlas.FireType.TimeStep):
    FireLog.logger.info(f"Running preprocessing code for {region[0]} at {tst=} with source {settings.FIRE_SOURCE}")
    FireRunDaskCoordinator.job_preprocess_region_t(None, None, region, tst)


if __name__ == "__main__":
    """ The main code to run preprocessing for a region and time period. It writes to a dedicated directory on s3.
    
    Example:
    python3 FireRunByRegionAndT.py --regnm="WesternUS" --tst="[2023,6,1,\"AM\"]"
    """
    
    parser = argparse.ArgumentParser()
    parser.add_argument("--regnm", type=str)
    parser.add_argument("--tst", type=validate_json)
    args = parser.parse_args()
    Run([args.regnm, None], args.tst)

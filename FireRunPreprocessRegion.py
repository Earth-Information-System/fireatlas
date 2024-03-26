import json
import argparse

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
    python3 FireRunPreprocessRegion.py --regnm="WesternUS" --bbox="[-119.5, 36.8, -118.9, 37.7]"
    """
    
    parser = argparse.ArgumentParser()
    parser.add_argument("--regnm", type=str)
    parser.add_argument("--bbox", type=validate_json)
    args = parser.parse_args()
    Run([args.regnm, args.bbox])

import json
import glob
import os
import argparse
import fireatlas

from dask import delayed
from itertools import chain
from fireatlas.utils import timed
from fireatlas import FireRunDaskCoordinator
from fireatlas import settings


def validate_json(s):
    try:
        print(f"${s}$")
        return json.loads(s)
    except ValueError:
        raise argparse.ArgumentTypeError("Not a valid JSON string")


@timed
def Run(region: fireatlas.FireTypes.Region, tst: fireatlas.FireTypes.TimeStep, ted: fireatlas.FireTypes.TimeStep):
    """
    """
    fire_forward_results = FireRunDaskCoordinator.job_fire_forward([None, ], region, tst, ted)

    # take all fire forward output and upload all snapshots/largefire outputs in parallel
    data_dir = os.path.join(settings.LOCAL_PATH, settings.OUTPUT_DIR, region[0], str(tst[0]))
    fgb_upload_results = [
        FireRunDaskCoordinator.concurrent_copy_outputs_from_local_to_s3(fire_forward_results, local_filepath)
        for local_filepath in list(chain(
            glob.glob(os.path.join(data_dir, "Snapshot", "*", "*.fgb")),
            glob.glob(os.path.join(data_dir, "Largefire", "*", "*.fgb"))
        ))
    ]
    dag = delayed(lambda x: x)(fgb_upload_results)
    dag.compute()

if __name__ == "__main__":
    """ The main code to run time forwarding for a time period

    Example:
    python3 FireRunFireForward.py --regnm="CaliTestRun" --tst="[2023,6,1,\"AM\"]" --ted="[2023,9,1,\"AM\"]"
    """

    parser = argparse.ArgumentParser()
    parser.add_argument("--regnm", type=str)
    parser.add_argument("--tst", type=validate_json)
    parser.add_argument("--ted", type=validate_json)
    args = parser.parse_args()
    Run([args.regnm, None], args.tst, args.ted)

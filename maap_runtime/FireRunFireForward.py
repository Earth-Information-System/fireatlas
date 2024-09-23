import json
import glob
import os
import argparse
from functools import partial
import fireatlas

from dask.distributed import Client

from fireatlas import FireRunDaskCoordinator
from fireatlas.FireIO import copy_from_local_to_s3, copy_from_local_to_veda_s3
from fireatlas.utils import timed
from fireatlas.postprocess import all_dir


def validate_json(s):
    try:
        print(f"${s}$")
        return json.loads(s)
    except ValueError:
        raise argparse.ArgumentTypeError("Not a valid JSON string")


@timed
def Run(region: fireatlas.FireTypes.Region, tst: fireatlas.FireTypes.TimeStep, ted: fireatlas.FireTypes.TimeStep, timeout_param: fireatlas.FireTypes.Timeout_param):
    """ run Fire_Forward and then upload outputs to s3 in parallel
    """
    FireRunDaskCoordinator.job_fire_forward(region, tst, ted, timeout_param)

    client = Client(n_workers=FireRunDaskCoordinator.MAX_WORKERS)

    # take all fire forward output and upload all outputs in parallel
    data_dir = all_dir(tst, region, location="local")
    fgb_s3_upload_futures = client.map(
        partial(copy_from_local_to_s3, fs=FireRunDaskCoordinator.fs),
        glob.glob(os.path.join(data_dir, "*", "*", "*.fgb"))
    )

    # take latest fire forward output and upload to VEDA S3 in parallel
    fgb_veda_upload_futures = client.map(
        partial(copy_from_local_to_veda_s3, fs=FireRunDaskCoordinator.fs, regnm=region[0]),
        glob.glob(os.path.join(data_dir, "*", f"{ted[0]}{ted[1]:02}{ted[2]:02}{ted[3]}", "*.fgb"))
    )
    # block until everything is uploaded
    client.gather([*fgb_s3_upload_futures, *fgb_veda_upload_futures])
    client.close()


if __name__ == "__main__":
    """ The main code to run time forwarding for a time period

    Example:
    python3 FireRunFireForward.py --regnm="CaliTestRun" --tst="[2023,6,1,\"AM\"]" --ted="[2023,9,1,\"AM\"]"
    """

    parser = argparse.ArgumentParser()
    parser.add_argument("--regnm", type=str)
    parser.add_argument("--tst", type=validate_json)
    parser.add_argument("--ted", type=validate_json)
    parser.add_argument("--timeout_param", type=validate_json)
    args = parser.parse_args()
    Run([args.regnm, None], args.tst, args.ted, args.timeout_param)

import json
import argparse

from utils import timed
from fireatlas import FireRunDaskCoordinator


def validate_json(s):
    try:
        return json.loads(s)
    except ValueError:
        raise argparse.ArgumentTypeError("Not a valid JSON string")


@timed
def Run():
    """ download data from different satellite sensors at the
    start of every 2nd hour from 1am through 11pm: `0 1-23/2 * * *`

    the jobs are being scheduled here: https://repo.ops.maap-project.org/eorland_gee/fireatlas_nrt/-/pipeline_schedules

    :return: None
    """
    FireRunDaskCoordinator.job_data_update_checker()



if __name__ == "__main__":
    """Downloads the NRT files at a regular interval
    
    Example:
    python3 FireRunDataUpdateChecker.py
    """
    Run()

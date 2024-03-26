import glob
import json
import argparse
import os

from utils import timed


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
    import FireIO
    import FireConsts
    from FireLog import logger
    import DataCheckUpdate

    try:
        # Download SUOMI-NPP
        DataCheckUpdate.update_VNP14IMGTDL()
        # Download NOAA-20
        DataCheckUpdate.update_VJ114IMGTDL()
    except Exception as exc:
        logger.exception(exc)
    finally:
        for filepath in glob.glob(os.path.join(FireConsts.get_dirprpdata(location="local"), "*", "*.txt")):
            FireIO.copy_from_local_to_s3(filepath)


if __name__ == "__main__":
    """Downloads the NRT files at a regular interval
    
    Example:
    python3 FireRunDataUpdateChecker.py
    """
    try:
        Run()
    except Exception as e:
        from FireLog import logger
        logger.exception(e)

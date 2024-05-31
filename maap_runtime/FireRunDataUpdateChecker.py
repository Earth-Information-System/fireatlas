import glob
from functools import partial
from datetime import datetime, UTC

from dask.distributed import Client
from fireatlas.utils import timed
from fireatlas import FireRunDaskCoordinator, settings
from fireatlas.FireIO import copy_from_local_to_s3


@timed
def Run():
    """ download data from different satellite sensors at the
    start of every 2nd hour from 1am through 11pm: `0 1-23/2 * * *`

    the jobs are being scheduled here: https://repo.ops.maap-project.org/eorland_gee/fireatlas_nrt/-/pipeline_schedules

    :return: None
    """
    ctime = datetime.now(UTC) 
    if tst in (None, "", []):  # if no start is given, run from beginning of year
        tst = [ctime.year, 1, 1, 'AM']

    if ted in (None, "", []):  # if no end time is given, set it as the most recent time
        if ctime.hour >= 18:
            ampm = 'PM'
        else:
            ampm = 'AM'
        ted = [ctime.year, ctime.month, ctime.day, ampm]

    client = Client(n_workers=FireRunDaskCoordinator.MAX_WORKERS)

    data_update_futures = FireRunDaskCoordinator.job_data_update_checker(client, tst, ted)
    client.gather(data_update_futures)

    data_upload_futures = [
        partial(copy_from_local_to_s3, fs=FireRunDaskCoordinator.fs),
        glob.glob(f"{settings.LOCAL_PATH}/{settings.PREPROCESSED_DIR}/*/*.txt")
    ]
    client.gather(data_upload_futures)
    client.close()


if __name__ == "__main__":
    """Downloads the NRT files at a regular interval

    Example:
    python3 FireRunDataUpdateChecker.py
    """
    Run()

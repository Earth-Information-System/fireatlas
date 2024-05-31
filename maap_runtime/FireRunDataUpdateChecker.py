import glob

from dask.distributed import Client
from fireatlas.utils import timed
from fireatlas import FireRunDaskCoordinator
from fireatlas import settings



@timed
def Run():
    """ download data from different satellite sensors at the
    start of every 2nd hour from 1am through 11pm: `0 1-23/2 * * *`

    the jobs are being scheduled here: https://repo.ops.maap-project.org/eorland_gee/fireatlas_nrt/-/pipeline_schedules

    :return: None
    """
    client = Client(n_workers=FireRunDaskCoordinator.MAX_WORKERS)

    data_update_futures = client.submit(FireRunDaskCoordinator.job_data_update_checker)
    client.gather(data_update_futures)

    data_upload_futures = [
        FireRunDaskCoordinator.concurrent_copy_from_local_to_s3,
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

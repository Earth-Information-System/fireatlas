import glob

from dask import delayed
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
    data_update_results = FireRunDaskCoordinator.job_data_update_checker()
    data_update_results.compute()

    data_upload_results = [
        FireRunDaskCoordinator.concurrent_copy_inputs_from_local_to_s3([None,], local_filepath)
        for local_filepath in glob.glob(f"{settings.LOCAL_PATH}/{settings.PREPROCESSED_DIR}/*/*.txt")
    ]
    dag = delayed(lambda x: x)(data_upload_results)
    dag.compute()


if __name__ == "__main__":
    """Downloads the NRT files at a regular interval

    Example:
    python3 FireRunDataUpdateChecker.py
    """
    Run()

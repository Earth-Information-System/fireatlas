import glob
import s3fs

from functools import partial
from dask.distributed import Client
from fireatlas.utils import timed
from fireatlas import FireRunDaskCoordinator, settings
from fireatlas.FireIO import copy_from_local_to_s3

# NOTE: this expects credentials to be resolvable globally
# via boto3/botocore common resolution paths
fs = s3fs.S3FileSystem(config_kwargs={"max_pool_connections": 10})

@timed
def Run():
    client = Client(n_workers=settings.N_DASK_WORKERS)

    data_update_futures = FireRunDaskCoordinator.job_nrt_current_day_updates(client)
    client.gather(data_update_futures)

    # uploads raw satellite files from `job_nrt_current_day_updates` in parallel
    data_upload_futures = client.map(
        partial(copy_from_local_to_s3, fs=fs),
        glob.glob(f"{settings.LOCAL_PATH}/{settings.PREPROCESSED_DIR}/*/*.txt")
    )
    # block until half-day timesteps and region are on s3
    timed(client.gather, text=f"Dask upload of {len(data_upload_futures) + 1} files")([*data_upload_futures])
    client.close()


if __name__ == "__main__":
    """Downloads the NRT files at a regular interval

    Example:
    python3 FireRunDataUpdateChecker.py
    """
    Run()

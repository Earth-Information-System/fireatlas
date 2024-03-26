import functools
import json
import argparse
import os
import time
from typing import Tuple
import concurrent
from concurrent.futures import ThreadPoolExecutor
import preprocess
from FireLog import logger
from FireTypes import Region, TimeStep
from utils import timed
from maap.maap import MAAP
from maap.utils import algorithm_utils
from maap.dps.dps_job import DPSJob

class RetryException(Exception):
    pass


class MaxRetryException(Exception):
    pass

class JobSubmissionException(Exception):
    pass


def retry(max_retries, exception_to_check, delay=1):

    def decorator_retry(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            retries = 0
            while True:
                try:
                    return func(*args, **kwargs)
                except exception_to_check as e:
                    retries += 1
                    if retries > max_retries:
                        raise MaxRetryException(e)
                    logger.info(f"Retry {retries} for {func.__name__} due to {e}. Retrying in {delay} seconds...")
                    time.sleep(delay)
        return wrapper
    return decorator_retry


def validate_json(s):
    try:
        return json.loads(s)
    except ValueError:
        raise argparse.ArgumentTypeError("Not a valid JSON string")


def get_algorithm_config_filepath(dir_names):
    current_file_dir = os.path.dirname(os.path.abspath(__file__))
    return [
        os.path.join(current_file_dir, 'maap_runtime', f'{dir_name}', 'algorithm_config.yaml')
        for dir_name in dir_names
    ]


def validate_job_submission(submitted_jobs: Tuple[DPSJob]) -> Tuple[DPSJob]:
    """we don't retry job submissions, they should ideally always work

    validate status of job submission results and return result 'job_id'
    """
    failed_statuses = [result for result in submitted_jobs if result.status == 'failed']
    if any(failed_statuses):
        raise JobSubmissionException(f"[ SUBMISSION FAILED ]: the following jobs failed to submit {failed_statuses}")
    return submitted_jobs


@retry(max_retries=2, exception_to_check=(RetryException,))
def submit_preprocess_region(
        maap_api: MAAP,
        algo_kwargs: dict,
        region: Region,
) -> Tuple[DPSJob]:
    import FireIO

    output_filepath = preprocess.preprocessed_region_filename(region, location="s3")
    if FireIO.os_path_exists(output_filepath):
        logger.info(f"skipping 'preprocess_region' b/c file \
        already exists for region {region[0]}, {output_filepath}")
        return

    submitted_jobs = []
    result = maap_api.submitJob(
        **algo_kwargs,
        **{
            "regnm": region[0],
            "bbox": json.dumps(region[1])
        }
    )
    submitted_jobs.append(result)
    failed_jobs = track_submitted_jobs(submitted_jobs)
    if failed_jobs:
        logger.warning(f"'submit_preprocess_region' attempt retry for job_ids={[job.id for job in failed_jobs]}")
        raise RetryException()


@retry(max_retries=2, exception_to_check=(RetryException,))
def submit_preprocess_per_t(
        maap_api: MAAP,
        algo_kwargs: dict,
        region: Region,
        list_of_time_steps: Tuple[TimeStep],
) -> Tuple[DPSJob]:
    import FireConsts
    import FireIO

    submitted_jobs = []
    for t in list_of_time_steps:
        output_filepath = preprocess.preprocessed_filename(t, sat=FireConsts.firesrc, region=region, location="s3")
        if FireIO.os_path_exists(output_filepath):
            logger.info(f"skipping 'preprocess_region_t' b/c file \
            already exists for region {region[0]}, {output_filepath}")
            continue

        result = maap_api.submitJob(
            **algo_kwargs,
            **{
                "regnm": region[0],
                "t": json.dumps(t)
            }
        )
        submitted_jobs.append(result)

    failed_jobs = track_submitted_jobs(submitted_jobs)
    if failed_jobs:
        logger.warning(f"'submit_preprocess_per_t' attempt retry for job_ids={[job.id for job in failed_jobs]}")
        raise RetryException()


@retry(max_retries=2, exception_to_check=(RetryException,))
def submit_update_checker(maap_api: MAAP, algo_kwargs: dict) -> Tuple[DPSJob]:
    # TODO: maybe check in the future
    submitted_jobs = []
    result = maap_api.submitJob(**algo_kwargs)
    submitted_jobs.append(result)

    failed_jobs = track_submitted_jobs(submitted_jobs)
    if failed_jobs:
        logger.warning(f"'submit_update_checker' attempt retry for job_ids={[job.id for job in failed_jobs]}")
        raise RetryException()


@retry(max_retries=2, exception_to_check=(RetryException,))
def submit_fire_forward(
        maap_api: MAAP,
        algo_kwargs: dict,
        region: Region,
        tst: TimeStep,
        ted: TimeStep,
) -> Tuple[DPSJob]:
    submitted_jobs = []
    result = maap_api.submitJob(
        **algo_kwargs,
        **{
            "regnm": region[0],
            "tst": json.dumps(tst),
            "ted": json.dumps(ted)
        }
    )
    submitted_jobs.append(result)
    failed_jobs = track_submitted_jobs(submitted_jobs)
    if failed_jobs:
        logger.warning(f"'submit_fire_forward' attempt retry for job_ids={[job.id for job in failed_jobs]}")
        raise RetryException()


def wait_for_job(dps_job: DPSJob) -> DPSJob:
    """this internal DPSJob function will block until job completes and use exponential backoff
    https://github.com/MAAP-Project/maap-py/blob/master/maap/dps/dps_job.py#L80C9-L80C28

    it seems the statuses.lower() are: ['failed', 'succeeded', 'accepted', 'running']
    https://github.com/MAAP-Project/maap-py/blob/master/maap/dps/dps_job.py
    """
    return dps_job.wait_for_completion()


def poll_on_job_status(jobs: Tuple[DPSJob]) -> Tuple[DPSJob]:
    failed_jobs = []
    # don't want to overwhelm the MAAP api so keeping max_workers relatively small
    with ThreadPoolExecutor(max_workers=5) as executor:
        dps_job_futures = [executor.submit(wait_for_job, dps_job) for dps_job in jobs]
        for dps_job in concurrent.futures.as_completed(dps_job_futures):
            try:
                if dps_job.result().retrieve_status().lower() != 'succeeded':
                    failed_jobs.append(dps_job)
            except Exception as e:
                logger.exception(f"'poll_on_jobs_status' failed with {e}")
    return failed_jobs


def track_submitted_jobs(submitted_jobs: Tuple[DPSJob]) -> Tuple[DPSJob]:
    queued_jobs = validate_job_submission(submitted_jobs)
    failed_jobs = poll_on_job_status(queued_jobs)
    return failed_jobs


@timed
def Run(region: Region, tst: TimeStep = None, ted: TimeStep = None):
    import FireTime

    dps_child_jobs = (
        'data_update_checker',
        'preprocess_region',
        'preprocess_region_t',
        'fire_forward'
    )
    list_of_time_steps = list(FireTime.t_generator(tst, ted))
    jobs_to_run = zip(dps_child_jobs, get_algorithm_config_filepath(dps_child_jobs))

    try:
        for dps_job_key, algo_config in jobs_to_run.items():
            maap_api = MAAP(maap_host='api.maap-project.org')
            algo_config = algorithm_utils.read_yaml_file(algo_config)
            submit_job_kwargs = {
                "identifier": f"job-{algo_config['algorithm_name']}:{algo_config['algorithm_version']}",
                "algo_id": algo_config["algorithm_name"],
                "version": algo_config["algorithm_version"],
                "username": "gcorradini",
                "queue": algo_config["queue"],
            }

            if dps_job_key == 'data_update_checker':
                submit_update_checker(maap_api, submit_job_kwargs)
            elif dps_job_key == 'preprocess_region':
                submit_preprocess_region(maap_api, submit_job_kwargs, region)
            elif dps_job_key == 'preprocess_region_and_t':
                submit_preprocess_per_t(maap_api, submit_job_kwargs, region, list_of_time_steps)
            elif dps_job_key == 'fire_forward':
                submit_fire_forward(maap_api, submit_job_kwargs, region, tst, ted)
            else:
                raise ValueError(f"[ JOB MATCH ]: dps_job_key='{dps_job_key}'")
    except (MaxRetryException, JobSubmissionException) as exc:
        logger.exception(exc)


if __name__ == "__main__":
    """Coordinating DPS jobs for all others
    
    Example:
    python3 FireRunDPSCoordinator.py --regnm="CaliTestRun" --bbox="[-125,36,-117,42]" --tst="[2023,6,1,\"AM\"]" --ted="[2023,9,1,\"AM\"]"
    """
    
    parser = argparse.ArgumentParser()
    parser = argparse.ArgumentParser()
    parser.add_argument("--regnm", type=str)
    parser.add_argument("--bbox", type=validate_json)
    parser.add_argument("--tst", type=validate_json)
    parser.add_argument("--ted", type=validate_json)
    args = parser.parse_args()
    Run([args.regnm, args.bbox], args.tst, args.ted)

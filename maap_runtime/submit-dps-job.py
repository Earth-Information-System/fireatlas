# NOTE: this script should only be run in gitlab CI
import uuid
import argparse
import json
import sys
from maap.maap import MAAP

parser = argparse.ArgumentParser()
parser.add_argument("algo", type=str, help="the name of the registered algorithm in MAAP")
parser.add_argument("github_ref", type=str, help="the git branch or tag to run against")
parser.add_argument("username", type=str, help="the username who submits the job")
parser.add_argument("queue", type=str, help="the queue where we want the job to run")
parser.add_argument("maap_environment", type=str, help="the MAAP image environment to run the job inside of")
parser.add_argument("--params", type=str, help="an optional string of json serialized additional params to "
                                               "pass to the job",
                    default="{}")


def submit_job(algo_id, github_ref, username, queue, maap_environment, params=None):
    maap = MAAP(maap_host='api.maap-project.org')
    if params is None:
        params = {}
    return maap.submitJob(
        identifier=f"job-{algo_id}:{github_ref}",
        algo_id=f"{algo_id}",  # MAAP seems to expect `submitJob` to identify things this way
        version=github_ref,
        username=username,
        queue=queue,
        **params
    )


if __name__ == '__main__':
    args = parser.parse_args()
    deserialized_params = json.loads(args.params)
    pargs = [args.algo, args.github_ref, args.username, args.queue, args.maap_environment]
    print(f"[ ARGS ]: {pargs}")
    print(f"[ KWARGS ]: {deserialized_params}")
    job = submit_job(*pargs, params=deserialized_params)
    print(job)
    if job.status == 'failed':
        msg = 'Job submission failed, please checkout logs for a message. If message is "Not Authorized" your PGT token might be old'
        print(msg)
        sys.exit(1)  # exit with a non-zero status to bubble up failure into GH actions



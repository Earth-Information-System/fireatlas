# NOTE: this script should only be run in gitlab CI
import uuid
import argparse
import json
from maap.maap import MAAP

parser = argparse.ArgumentParser()
parser.add_argument("algo", type=str, help="the name of the registered algorithm in MAAP")
parser.add_argument("version", type=str, help="the git branch to run against")
parser.add_argument("username", type=str, help="the username who submits the job")
parser.add_argument("queue", type=str, help="the queue where we want the job to run")
parser.add_argument("maap_environment", type=str, help="the MAAP image environment to run the job inside of")
parser.add_argument("--params", type=str, help="an optional string of json serialized additional params to "
                                               "pass to the job",
                    default="{}")


def submit_job(algo_id, version, username, queue, maap_environment, params=None):
    maap = MAAP(maap_host='api.maap-project.org')
    if params is None:
        params = {}
    return maap.submitJob(
        identifier=f"job-{algo_id}:{version}",
        algo_id=f"{algo_id}",  # MAAP seems to expect `submitJob` to identify things this way
        version=version,
        username=username,
        queue=queue,
        **params
    )


if __name__ == '__main__':
    args = parser.parse_args()
    deserialized_params = json.loads(args.params)
    pargs = [args.algo, args.version, args.username, args.queue, args.maap_environment]
    print(f"[ ARGS ]: {pargs}")
    print(f"[ KWARGS ]: {deserialized_params}")
    response = submit_job(*pargs, params=deserialized_params)
    print(response)



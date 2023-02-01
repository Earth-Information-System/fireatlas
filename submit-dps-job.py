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
parser.add_argument("--params", type=str, help="an optional json serialized string of additonal params to "
                                               "pass to the job",
                    default="{}")


def submit_job(algo_id, version, username, queue, params={}):
    maap = MAAP(maap_host='api.ops.maap-project.org')
    return maap.submitJob(
        identifier=f"job-{algo_id}_ubuntu:{version}",
        algo_id=algo_id,
        version=version,
        username=username,
        queue=queue,
        **params
    )


if __name__ == '__main__':
    args = parser.parse_args()
    deserialized_params = json.loads(args.params)
    pargs = [args.algo, args.version, args.username, args.queue]
    print(f"[ ARGS ]: {pargs}")
    response = submit_job(*pargs, params=deserialized_params)
    print(response)



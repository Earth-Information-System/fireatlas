import json
import argparse
import fireatlas
from fireatlas.utils import timed
from fireatlas import FireRunDaskCoordinator


def validate_json(s):
    try:
        return json.loads(s)
    except ValueError:
        raise argparse.ArgumentTypeError("Not a valid JSON string")


@timed
def Run(region: fireatlas.FireTypes.Region, tst: fireatlas.FireTypes.TimeStep, ted:  fireatlas.FireTypes.TimeStep):
    # NOTE: this set up has to happen before `import FireConsts`
    # or any other modules also import from FireConsts
    # so set os environ variables that will override
    # `FireConsts.settings` for all Python interpreters
    # (meaning even those spawned during fork in multiple processes)
    # os.environ['EPSG_CODE'] = FireEnums.EPSG.HI_LAT
    # os.environ['FTYP_OPT'] = 2
    # os.environ['CONT_OPT'] = 2
    # import FireConsts
    FireRunDaskCoordinator.job_fire_forward([None], region, tst, ted)


if __name__ == "__main__":
    """ The main code to run time forwarding for a time period
    
    Example:
    python3 FireRunFireForward.py --regnm="CaliTestRun" --tst="[2023,6,1,\"AM\"]" --ted="[2023,9,1,\"AM\"]"
    """
    
    parser = argparse.ArgumentParser()
    parser.add_argument("--regnm", type=str)
    parser.add_argument("--tst", type=validate_json)
    parser.add_argument("--ted", type=validate_json)
    args = parser.parse_args()
    Run([args.regnm, None], args.tst, args.ted)

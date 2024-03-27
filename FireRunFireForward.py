import os
import glob
import json
from datetime import datetime

import argparse
from FireTypes import Region
from utils import timed


def validate_json(s):
    try:
        return json.loads(s)
    except ValueError:
        raise argparse.ArgumentTypeError("Not a valid JSON string")


@timed
def Run(region: Region, tst=None, ted=None):
    # NOTE: this set up has to happen before `import FireConsts`
    # or any other modules also import from FireConsts
    # so set os environ variables that will override
    # `FireConsts.settings` for all Python interpreters
    # (meaning even those spawned during fork in multiple processes)
    # os.environ['EPSG_CODE'] = FireEnums.EPSG.HI_LAT
    # os.environ['FTYP_OPT'] = 2
    # os.environ['CONT_OPT'] = 2
    # import FireConsts
    import FireIO, FireConsts, FireMain, postprocess
    from FireLog import logger

    ctime = datetime.now()
    
    if tst in (None, ""): # if no start is given, run from beginning of year
        tst = [ctime.year, 1, 1, 'AM']

    if ted in (None, ""): # if no end time is given, set it as the most recent time
        
        if ctime.hour >= 18:
            ampm = 'PM'
        else:
            ampm = 'AM'

        ted = [ctime.year, ctime.month, ctime.day, ampm]
    
    logger.info(f"Running code for {region[0]} from {tst} to {ted} with source {FireConsts.firesrc}")

    allfires, allpixels = FireMain.Fire_Forward(tst=tst, ted=ted, restart=False, region=region)
    allpixels_filepath = postprocess.save_allpixels(allpixels, tst, ted, region)
    allfires_filepath = postprocess.save_allfires_gdf(allfires.gdf, tst, ted, region)

    FireIO.copy_from_local_to_s3(allpixels_filepath)
    FireIO.copy_from_local_to_s3(allfires_filepath)

    postprocess.save_snapshots(allfires.gdf, region, tst, ted)

    large_fires = postprocess.find_largefires(allfires.gdf)
    postprocess.save_large_fires_nplist(allpixels, region, large_fires, tst)
    postprocess.save_large_fires_layers(allfires.gdf, region, large_fires, tst)

    for filepath in glob.glob(os.path.join(FireConsts.get_diroutdata(location="local"), region[0], str(tst[0]), "Snapshot", "*", "*.fgb")):
        FireIO.copy_from_local_to_s3(filepath)

    for filepath in glob.glob(os.path.join(FireConsts.get_diroutdata(location="local"), region[0], str(tst[0]), "Largefire", "*", "*.fgb")):
        FireIO.copy_from_local_to_s3(filepath)


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

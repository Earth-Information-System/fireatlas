import sys
import os
import json
import time
import argparse
from FireLog import logger
from datetime import datetime
import FireConsts

def validate_json(s):
    try:
        return json.loads(s)
    except ValueError:
        raise argparse.ArgumentTypeError("Not a valid JSON string")


def RegionRun(regnm, bbox, tst=None, ted=None):
    # NOTE: this set up has to happen before `import FireConsts`
    # or any other modules also import from FireConsts
    # so set os environ variables that will override
    # `FireConsts.settings` for all Python interpreters
    # (meaning even those spawned during fork in multiple processes)
    # os.environ['EPSG_CODE'] = FireEnums.EPSG.HI_LAT
    # os.environ['FTYP_OPT'] = 2
    # os.environ['CONT_OPT'] = 2
    # import FireConsts
    import FireIO, FireMain, FireGpkg, FireGpkg_sfs

    ctime = datetime.now()

    region = (regnm,bbox)
    
    if tst in (None, ""): # if no start is given, parse from the most recent run
        
        lts = FireIO.get_lts_serialization(regnm=region[0])
        if lts == None: # if no prior runs, start anew
            tst = [ctime.year, 1, 1, 'AM']
        else:
            tst = lts

    if ted in (None, ""): # if no end time is given, set it as the most recent time
        
        if ctime.hour >= 18:
            ampm = 'PM'
        else:
            ampm = 'AM'

        ted = [ctime.year, ctime.month, ctime.day, ampm]
    
    print(f"Running code for {region[0]} from {tst} to {ted} with source {FireConsts.firesrc}")

    FireMain.Fire_Forward(tst=tst, ted=ted, restart=False, region=region)
    FireGpkg.save_gdf_trng(tst=tst, ted=ted, regnm=region[0])
    FireGpkg_sfs.save_sfts_trng(tst, ted, regnm=region[0])


if __name__ == "__main__":
    """ The main code to run time forwarding for a time period
    """
    from dask.distributed import performance_report
    
    parser = argparse.ArgumentParser()
    parser.add_argument("--regnm", type=str)
    parser.add_argument("--bbox", type=validate_json)
    parser.add_argument("--tst", type=validate_json)
    parser.add_argument("--ted", type=validate_json)
    parser.add_argument("--export", action="store_true", help="export to veda or not")
    args = parser.parse_args()
    t1 = time.time()
    try:
        FireConsts.export_to_veda = args.export
        RegionRun(args.regnm, args.bbox, args.tst, args.ted)
    except Exception as e:
        logger.exception(e)
    t2 = time.time()
    print(f"{(t2-t1)/60.} minutes used to run the whole code")


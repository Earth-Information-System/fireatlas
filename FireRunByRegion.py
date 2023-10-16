""" FireRun
Module to control different runs
"""
import sys
import time
import argparse
from FireLog import logger


def RegionRun(regnm, bbox, tst=None, ted=None):
    
    import FireIO, FireMain, FireGpkg, FireGpkg_sfs, FireObj
    import FireConsts
    from datetime import datetime
    import os
        
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
    parser.add_argument("regnm")
    parser.add_argument("bbox")
    parser.add_argument("tst")
    parser.add_argument("ted")
    args = parser.parse_args()
    t1 = time.time()
    try:
        RegionRun(args.regnm, args.bbox, args.tst, args.ted)
    except Exception as e:
        logger.exception(e)
    t2 = time.time()
    print(f"{(t2-t1)/60.} minutes used to run the whole code")


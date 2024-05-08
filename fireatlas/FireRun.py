""" FireRun
Module to control different runs
"""
import time
import argparse
from datetime import datetime
import geopandas as gpd
from pandas import Series
import os
from fireatlas.FireLog import logger
from fireatlas.utils import timed
from fireatlas import settings


def Yearbatchrun(year, tst=None, ted=None, restart=False):
    """ Run the code for each single year
    """
    from fireatlas import FireMain, FireSummary, FireGdf_merge, FireGdf_sfs_merge, FireGdf_ign, FireGdf_final

    t1 = time.time()
    # set the start and end time
    if tst is None:
        tst = (year, 6, 1, "AM")
    if ted is None:
        ted = (year, 8, 31, "PM")
    # if year == 2012: tst = (year,1,20,'AM')

    # Run the time forward and record daily fire objects .pkl data and fire attributes .GeoJSON data
    FireMain.Fire_Forward(tst=tst, ted=ted, restart=restart, region="AK")
    t2 = time.time()
    print(f"{(t2-t1)/60.} minutes to run algorithm")

    # Run to save geojson files for each time step
    FireGdf_merge.save_gdf_trng(tst=tst, ted=ted, fperim=True)
    # FireGdf_merge.save_gdf_trng(tst=tst,ted=ted,fall=True)
    FireGdf_merge.save_gdf_trng(tst=tst, ted=ted, NFP_txt=True)
    t3 = time.time()
    print(f"{(t3-t2)/60.} minutes to save gpkg files")

    # Run to save ignition point layer for each time step
    FireGdf_ign.save_gdf_trng(tst, ted)
    FireGdf_final.save_gdf_trng(tst, ted)
    t31 = time.time()
    print(f"{(t31-t3)/60.} minutes to save ognitions and final perimeters")

    # Run to save large fire geosjon files for each time step
    # FireGdf_sfs.save_gdf_trng(tst=tst,ted=ted,fperim=True)
    FireGdf_sfs_merge.save_gdf_trng(ted=ted, fperim=True)
    t4 = time.time()
    print(f"{(t4-t31)/60.} minutes to save large fires")

    # Run to save year end summary and file lists
    FireSummary.save_sum_1d(tst, ted)
    FireSummary.add_heritage(ted)
    FireSummary.add_largefirelist(ted)
    t5 = time.time()
    print(f"{(t5-t4)/60.} minutes to generate summary")

    # Run to clean up large fire geojson files not needed
    # FireGdf_sfs.yrend_clean(ted)

    t2 = time.time()
    print(f"{(t2-t1)/60.} minutes to run code")


def CreekSamplerun(firesrc='SNPP'):
    """
    - Run Creek fire tracking
    - firesrc options: 'SNPP', 'NOAA20', 'VIIRS'
    - before running, need to set corresponding firesrc in FireConsts.py
    """
    from fireatlas import FireMain, FireGpkg, FireGpkg_sfs

    tst = (2020, 9, 5, "AM")
    ted = (2020, 9, 8, "PM")
    region = ("CreekMAAPDEMO"+firesrc, [-119.5, 36.8, -118.9, 37.7])

    settings.FIRE_SOURCE = firesrc

    # # do fire tracking
    FireMain.Fire_Forward(tst=tst, ted=ted, restart=True, region=region)
    #
    # # calculate and save snapshot files
    FireGpkg.save_gdf_trng(tst=tst, ted=ted, regnm=region[0])
    #
    # # calculate and save single fire files
    FireGpkg_sfs.save_sfts_trng(tst, ted, regnm=region[0])

    #FireGpkg_sfs.convert_sfts(region[0],2020,[0])

def DixieSamplerun(firesrc='SNPP'):
    """
    - Run Creek fire tracking
    - firesrc options: 'SNPP', 'NOAA20', 'VIIRS'
    - before running, need to set corresponding firesrc in FireConsts.py
    """
    from fireatlas import FireGpkg_sfs

    tst = (2021, 7, 13, "AM")
    ted = (2021, 9, 16, "PM")
    region = ("Dixie"+firesrc, [-121.6, 39.8, -120.1, 40.8])

    settings.FIRE_SOURCE = firesrc

    # # do fire tracking
    # FireMain.Fire_Forward(tst=tst, ted=ted, restart=True, region=region)
    #
    # # calculate and save snapshot files
    # FireGpkg.save_gdf_trng(tst=tst, ted=ted, regnm=region[0])
    #
    # # calculate and save single fire files
    # FireGpkg_sfs.save_sfts_trng(tst, ted, regnm=region[0])

    FireGpkg_sfs.convert_sfts(region[0],2021,[5])
# def CreekSamplerunNOAA20():
#     """
#     need to change firesrc in FireConsts.py 'NOAA20'
#     """
#     import FireMain, FireGpkg, FireGpkg_sfs
#
#     tst = (2020, 9, 5, "AM")
#     ted = (2020, 11, 5, "PM")
#     region = ("CreekNOAA20", [-119.5, 36.8, -118.9, 37.7])
#
#     # do fire tracking
#     FireMain.Fire_Forward(tst=tst, ted=ted, restart=True, region=region)
#
#     # calculate and save snapshot files
#     FireGpkg.save_gdf_trng(tst=tst, ted=ted, regnm=region[0])
#
#     # calculate and save single fire files
#     FireGpkg_sfs.save_sfts_trng(tst, ted, regnm=region[0])
#
#
# def CreekSamplerunVIIRS():
#     """
#     need to change firesrc in FireConsts.py 'VIIRS'
#     """
#     import FireMain, FireGpkg, FireGpkg_sfs
#
#     tst = (2020, 9, 5, "AM")
#     ted = (2020, 11, 5, "PM")
#     region = ("CreekVIIRS", [-119.5, 36.8, -118.9, 37.7])
#
#     # do fire tracking
#     FireMain.Fire_Forward(tst=tst, ted=ted, restart=True, region=region)
#
#     # calculate and save snapshot files
#     FireGpkg.save_gdf_trng(tst=tst, ted=ted, regnm=region[0])
#
#     # calculate and save single fire files
#     FireGpkg_sfs.save_sfts_trng(tst, ted, regnm=region[0])


def CreekRegionSamplerun():
    from fireatlas import FireMain, FireGpkg, FireGpkg_sfs

    tst = (2020, 9, 5, "AM")
    ted = (2020, 9, 19, "AM")
    region = ("CreekEliTwoWeeksSNPP", [-120, 36, -118, 38])
    logger.info(f'STARTING RUN FOR {region[0]}')
    tstart = time.time()

    
    # do fire tracking
    FireMain.Fire_Forward(tst=tst, ted=ted, restart=True, region=region)

    # calculate and save snapshot files
    FireGpkg.save_gdf_trng(tst=tst, ted=ted, regnm=region[0])

    # calculate and save single fire files
    FireGpkg_sfs.save_sfts_trng(tst, ted, regnm=region[0])
    tend = time.time()
    
    logger.info(f"{(tend-tstart)/60.} minutes used for CreekRegionSamplerun with dask.")
    

def ChileSampleRun():
    from fireatlas import FireMain, FireGpkg, FireGpkg_sfs

    tst = (2023, 2, 13, 'AM')
    ted = (2023, 4, 1, "AM")
    region = ("ChileGlobalNRT", [-82.39770817214531, -55.54947848623975, 
                                 -65.87427067214531, -14.895459243377251])
    logger.info(f'STARTING RUN FOR {region[0]}')
    tstart = time.time()

    
    # do fire tracking
    FireMain.Fire_Forward(tst=tst, ted=ted, restart=True, region=region)

    # calculate and save snapshot files
    FireGpkg.save_gdf_trng(tst=tst, ted=ted, regnm=region[0])

    # calculate and save single fire files
    FireGpkg_sfs.save_sfts_trng(tst, ted, regnm=region[0])
    tend = time.time()
    
    logger.info(f"{(tend-tstart)/60.} minutes used for ChileSampleRun with dask.")


def BorealNA():
    # NOTE: this set up has to happen before `from fireatlas import settings`
    # so set os environ variables that will override
    # `FireConsts.Settings` for all Python interpreters
    # (meaning even those spawned during fork in multiple processes)
    os.environ['EPSG_CODE'] = 3571
    os.environ['FTYP_OPT'] = "global"
    os.environ['CONT_OPT'] = "global"
    
    from fireatlas import FireIO, FireMain, FireGpkg, FireGpkg_sfs

    ctime = datetime.now()
    region = ("BOREAL_NRT_3571", [-169, 44, -48, 75])
    
    logger.info(f'STARTING RUN FOR {region[0]}')
    tstart = time.time()
    
    lts = FireIO.get_lts_serialization(regnm=region[0])
    if lts == None:
        tst = [ctime.year, 1, 1, 'AM']
    else:
        #tst = FireObj.t_nb(lts, nb="previous") <-- this returns an error
        tst = lts

    if ctime.hour >= 18:
        ampm = 'PM'
    else:
        ampm = 'AM'
    #tst = [2023,7,22,'AM']
    ted = [ctime.year, ctime.month, ctime.day, ampm]
    #ted = [2023,6,28,'AM']
    print(f"Running code from {tst} to {ted}.")
    
    # do fire tracking
    FireMain.Fire_Forward(tst=tst, ted=ted, restart=False, region=region)

    # calculate and save snapshot files
    FireGpkg.save_gdf_trng(tst=tst, ted=ted, regnm=region[0])

    # calculate and save single fire files
    FireGpkg_sfs.save_sfts_trng(tst, ted, regnm=region[0])
    
    tend = time.time()
    
    logger.info(f"{(tend-tstart)/60.} minutes used for BorealNA for {tst} to {ted}.")
    

def ItalyGreeceNRT():
    from fireatlas import FireIO, FireMain, FireGpkg, FireGpkg_sfs

    ctime = datetime.now()
    region = ("ItalyGreeceNRT_DPS", [11, 36, 28, 42])
    
    logger.info(f'STARTING RUN FOR {region[0]}')
    tstart = time.time()
    
    lts = FireIO.get_lts_serialization(regnm=region[0])
    if lts == None:
        tst = [ctime.year, 1, 1, 'AM']
    else:
        #tst = FireObj.t_nb(lts, nb="previous") <-- this returns an error
        tst = lts

    if ctime.hour >= 18:
        ampm = 'PM'
    else:
        ampm = 'AM'
    #tst = [2023,6,1,'AM']
    ted = [ctime.year, ctime.month, ctime.day, ampm]
    #ted = [2023,6,28,'AM']
    print(f"Running code from {tst} to {ted}.")
    
    # do fire tracking
    FireMain.Fire_Forward(tst=tst, ted=ted, restart=False, region=region)

    # calculate and save snapshot files
    FireGpkg.save_gdf_trng(tst=tst, ted=ted, regnm=region[0])

    # calculate and save single fire files
    FireGpkg_sfs.save_sfts_trng(tst, ted, regnm=region[0])
    
    tend = time.time()
    
    logger.info(f"{(tend-tstart)/60.} minutes used for Italy and Greece for {tst} to {ted}.")    



def QuebecSampleRun():
    from fireatlas import FireIO, FireMain, FireGpkg, FireGpkg_sfs
    ctime = datetime.now()
    #tst = (2023, 1, 1, 'AM')
    #ted = (2023, 6, 7, "AM")
    region = ("QuebecGlobalNRT_ELI", [-83.69877641421793, 44.25483911637959, 
                                      -48.45463578921794, 62.94135765648493])
    
    logger.info(f'STARTING RUN FOR {region[0]}')
    tstart = time.time()
    
    lts = FireIO.get_lts_serialization(regnm=region[0])
    if lts == None:
        tst = [ctime.year, 1, 1, 'AM']
    else:
        #tst = FireObj.t_nb(lts, nb="previous") <-- this returns an error
        tst = lts

    if ctime.hour >= 18:
        ampm = 'PM'
    else:
        ampm = 'AM'
    #tst = [2023,6,30,'AM']
    ted = [ctime.year, ctime.month, ctime.day, ampm]
    #ted = [2023,6,28,'AM']
    print(f"Running code from {tst} to {ted}.")
    
    # do fire tracking
    FireMain.Fire_Forward(tst=tst, ted=ted, restart=False, region=region)

    # calculate and save snapshot files
    FireGpkg.save_gdf_trng(tst=tst, ted=ted, regnm=region[0])

    # calculate and save single fire files
    FireGpkg_sfs.save_sfts_trng(tst, ted, regnm=region[0])
    
    tend = time.time()
    
    logger.info(f"{(tend-tstart)/60.} minutes used for QuebecSampleRun with dask.")

def CArun():
    from fireatlas import FireIO, FireMain, FireGpkg, FireGpkg_sfs

    CAshp = FireIO.get_Cal_shp()
    region = ("California", CAshp)

    tst = (2019, 6, 1, "AM")
    ted = (2019, 11, 30, "PM")

    # do fire tracking
    FireMain.Fire_Forward(tst=tst, ted=ted, restart=True, region=region)

    # calculate and save snapshot files
    FireGpkg.save_gdf_trng(tst=tst, ted=ted, regnm=region[0])

    # calculate and save single fire files
    FireGpkg_sfs.save_sfts_trng(tst, ted, regnm=region[0])

def CArun2136():
    from fireatlas import FireIO, FireMain, FireGpkg, FireGpkg_sfs

    CAshp = FireIO.get_Cal_shp()
    region = ("California2136", CAshp)

    tst = (2020, 9, 23, "AM")
    ted = (2020, 9, 30, "PM")

    # do fire tracking
    FireMain.Fire_Forward(tst=tst, ted=ted, restart=True, region=region)

    # calculate and save snapshot files
    FireGpkg.save_gdf_trng(tst=tst, ted=ted, regnm=region[0])

    # calculate and save single fire files
    FireGpkg_sfs.save_sfts_trng(tst, ted, regnm=region[0])

def CArun32610():
    from fireatlas import FireIO, FireMain, FireGpkg, FireGpkg_sfs

    CAshp = FireIO.get_Cal_shp()
    region = ("California32610", CAshp)

    tst = (2012, 1, 20, "AM")  # temporarily set to 'pm' to avoid reading from previous year's data
    # tst = (2012, 1, 1, "AM")

    ted = (2012, 12, 31, "PM")

    # do fire tracking
    # FireMain.Fire_Forward(tst=tst, ted=ted, restart=True, region=region)

    # calculate and save snapshot files
    FireGpkg.save_gdf_trng(tst=tst, ted=ted, regnm=region[0])

#     calculate and save single fire files
    FireGpkg_sfs.save_sfts_trng(tst, ted, regnm=region[0], layers=["perimeter", "fireline", "newfirepix","nfplist"])

    FireGpkg_sfs.combine_sfts(region[0],tst[0],addFRAP=True)

def CArunNRT():
    
    from fireatlas import FireIO, FireMain, FireGpkg, FireGpkg_sfs, FireConsts
    ctime = datetime.now()

    CAshp = FireIO.get_Cal_shp()
    region = ("CA", CAshp)

    lts = FireIO.get_lts_serialization(regnm=region[0])
    if lts == None:
        tst = [ctime.year, 1, 1, 'AM']
    else:
        #tst = FireObj.t_nb(lts, nb="previous") <-- this returns an error
        tst = lts

    if ctime.hour >= 18:
        ampm = 'PM'
    else:
        ampm = 'AM'
    #tst = [2022,1,1,'AM']
    ted = [ctime.year, ctime.month, ctime.day, ampm]
    #ted = [2022,1,10,'AM']
    print(f"Running code from {tst} to {ted}.")

    #FireMain.Fire_Forward(tst=tst, ted=ted, restart=False, region=region)
    #FireGpkg.save_gdf_trng(tst=tst, ted=ted, regnm=region[0])
    FireGpkg_sfs.save_sfts_trng(tst, ted, regnm=region[0])
    

@timed
def CONUSrunNRT():
    from fireatlas import preprocess
    from fireatlas import FireConsts, FireMain, FireGpkg, FireGpkg_sfs

    region = ["CONUS",]
    tst = [2023, 8, 28, 'AM']
    ted = [2023, 9, 6, 'AM']

    logger.info(f'STARTING RUN FOR {region[0]}')    
    
    return FireMain.Fire_Forward(tst=tst, ted=ted, restart=False, region=region)


def WesternUSrunNRT():
    
    from fireatlas import FireIO, FireMain, FireGpkg, FireGpkg_sfs
    if settings.FIRE_NRT != True:
        print('Please set FIRE_NRT to True')
        return
    
    ctime = datetime.now()

    region = ('WesternUSNRT',[-125.698046875,31.176476158707615,
                              -101.00078125,49.51429477264348])

    lts = FireIO.get_lts_serialization(regnm=region[0])
    if lts == None:
        tst = [ctime.year, 1, 1, 'AM']
    else:
        #tst = FireObj.t_nb(lts, nb="previous") <-- this returns an error
        tst = lts

    if ctime.hour >= 18:
        ampm = 'PM'
    else:
        ampm = 'AM'
    
    ted = [ctime.year, ctime.month, ctime.day, ampm]
    #ted = [2022,1,10,'AM']
    print(f"Running code from {tst} to {ted} with source {settings.FIRE_SOURCE}")

    FireMain.Fire_Forward(tst=tst, ted=ted, restart=False, region=region)
    FireGpkg.save_gdf_trng(tst=tst, ted=ted, regnm=region[0])
    FireGpkg_sfs.save_sfts_trng(tst, ted, regnm=region[0])
    
    
def SouthEastUSrunNRT():
    from fireatlas import FireIO, FireMain, FireGpkg, FireGpkg_sfs

    if settings.FIRE_NRT != True:
        print('Please set FIRE_NRT to True')
        return
    
    ctime = datetime.now()

    region = ('SouthEastUSNRT_DPS',[-106.79802059770478,24.457626666909054,
                                    -72.87223934770478,37.309430118635944])
    
    logger.info(f'STARTING RUN FOR {region[0]}')
    
    lts = FireIO.get_lts_serialization(regnm=region[0])
    if lts == None:
        tst = [ctime.year, 1, 1, 'AM']
    else:
        #tst = FireObj.t_nb(lts, nb="previous") <-- this returns an error
        tst = lts

    if ctime.hour >= 18:
        ampm = 'PM'
    else:
        ampm = 'AM'
    
    ted = [ctime.year, ctime.month, ctime.day, ampm]
    #ted = [2022,1,10,'AM']
    print(f"Running code from {tst} to {ted} with source {settings.FIRE_SOURCE}")

    FireMain.Fire_Forward(tst=tst, ted=ted, restart=False, region=region)
    FireGpkg.save_gdf_trng(tst=tst, ted=ted, regnm=region[0])
    FireGpkg_sfs.save_sfts_trng(tst, ted, regnm=region[0])

def SouthEastUS_LF_ONLY():
    from fireatlas import FireIO, FireMain, FireGpkg, FireGpkg_sfs
    if settings.FIRE_NRT != True:
        print('Please set FIRE_NRT to True')
        return
    
    ctime = datetime.now()

    region = ('SouthEastUSNRT_DPS',[-106.79802059770478,24.457626666909054,
                                    -72.87223934770478,37.309430118635944])
    
    logger.info(f'STARTING RUN FOR {region[0]}')
    
    lts = FireIO.get_lts_serialization(regnm=region[0])
    if lts == None:
        tst = [ctime.year, 1, 1, 'AM']
    else:
        #tst = FireObj.t_nb(lts, nb="previous") <-- this returns an error
        tst = lts

    if ctime.hour >= 18:
        ampm = 'PM'
    else:
        ampm = 'AM'
    tst = [2023,2,22,'PM']
    ted = [2023,2,22,'PM']
    #ted = [ctime.year, ctime.month, ctime.day, ampm]
    #ted = [2022,1,10,'AM']
    print(f"Running code from {tst} to {ted} with source {settings.FIRE_SOURCE}")

    #FireMain.Fire_Forward(tst=tst, ted=ted, restart=False, region=region)
    #FireGpkg.save_gdf_trng(tst=tst, ted=ted, regnm=region[0])
    FireGpkg_sfs.save_sfts_trng(tst, ted, regnm=region[0])

    
def NorthEastUSrunNRT():
    from fireatlas import FireIO, FireMain, FireGpkg, FireGpkg_sfs

    if settings.FIRE_NRT != True:
        print('Please set FIRE_NRT to True')
        return
    
    ctime = datetime.now()

    region = ('NorthEastUSNRT_DPS',[-106.79802059770478,35.590087054959234,
                                    -66.11829033856952,49.628319367544776])

    lts = FireIO.get_lts_serialization(regnm=region[0])
    if lts == None:
        tst = [ctime.year, 1, 1, 'AM']
    else:
        #tst = FireObj.t_nb(lts, nb="previous") <-- this returns an error
        tst = lts

    if ctime.hour >= 18:
        ampm = 'PM'
    else:
        ampm = 'AM'
    
    #tst = [ctime.year, 1, 1, 'AM']    
    ted = [ctime.year, ctime.month, ctime.day, ampm]
    #ted = [2022,1,10,'AM']
    print(f"Running code from {tst} to {ted} with source {settings.FIRE_SOURCE}")

    FireMain.Fire_Forward(tst=tst, ted=ted, restart=False, region=region)
    FireGpkg.save_gdf_trng(tst=tst, ted=ted, regnm=region[0])
    FireGpkg_sfs.save_sfts_trng(tst, ted, regnm=region[0])


def WesternUSYrRun(year):
    from fireatlas import FireIO, FireMain, FireGpkg, FireGpkg_sfs, FireConsts
    region = ('WesternUS',[-125.698046875,31.676476158707615,
                           -101.00078125,49.51429477264348])

    tst = (year, 1, 1, "AM")
    ted = (year, 12, 31, "PM")
    
    print(f"Running code from {tst} to {ted} with source {settings.FIRE_SOURCE}")
    
    FireMain.Fire_Forward(tst=tst, ted=ted, restart=False, region=region)
    FireGpkg.save_gdf_trng(tst=tst, ted=ted, regnm=region[0])
    FireGpkg_sfs.save_sfts_trng(tst, ted, regnm=region[0])


def constrainByShape_Fire_Forward(tst, ted, perimeter_gdf, sat='SNPP'):
    """
    Derives region from the perimeter and updates tst/ted based on the fire start/end dates.
    preprocesses the monthly data for the region and runs fire_forward.
    Saves the outpts too--allpixels, allfires, and individual fires.
    """
    from fireatlas import FireTime, preprocess, FireMain, postprocess
    from fireatlas import settings

    settings.FIRE_SOURCE = sat

    # convert gdf to a series if it is not already
    if not isinstance(perimeter_gdf, Series):
        perimeter = perimeter_gdf.iloc[0]
    else:
        perimeter = perimeter_gdf

    # update the time window based on the fire start/end dates
    tst, ted = FireTime.update_tst_ted(perimeter, tst, ted)

    # define region based on the perimeter
    region_name = f"{perimeter['FIRE_NAME']}_{perimeter['FIRE_ID']}"
    region = (region_name, perimeter.geometry)

    logger.info(f"=============== Running: {region_name} ===============")
    # preprocess the monthly data for this region
    # do once per region--removes static flare sources
    preprocess.preprocess_region(region)

    # get lists of times to run fire_forward
    list_of_ts = list(FireTime.t_generator(tst, ted))
    unique_ym = preprocess.check_preprocessed_file(tst, ted, sat)

    # preprocess the monthly files--will only do for those not already processed
    for ym in unique_ym:
        preprocess.preprocess_monthly_file(ym, sat)

    # filter VIIRS to the perimeter for each time step
    region = preprocess.read_region(region)
    for t in list_of_ts:
        preprocess.preprocess_region_t(t, region=region)

    # finally run fire_forward
    allfires, allpixels = FireMain.Fire_Forward(tst=tst, ted=ted, restart=False, region=region)

    # Save outputs
    postprocess.save_allpixels(allpixels, tst, ted, region)
    postprocess.save_allfires_gdf(allfires.gdf, tst, ted, region)
    postprocess.save_individual_fire(allfires.gdf, tst, ted, region)


def constrainByShape_Run(perimeter_gdf_path, tst=None, ted=None, sat = 'SNPP', data_source='FRAP', id_col='INC_NUM'):
    """
    Runs fire_forward within the specified time range and constrains viirs pixels to the supplied fire perimeters.
    Accomodates multiple perimeters (will loop through each one) or a single perimeter.
    Time window further filtered based on the fire start/end dates.

    Parameters
    ----------
    tst : list
        Min start time (year, month, day, "AM" or "PM").
    ted : list
        Max end time (year, month, day, "AM" or "PM").
    perimeter_gdf_path : path to geopandas.GeoDataFrame of fire perimeters
    sat: 'SNPP'--currently the only option for montly VIIRS data

    Returns
    -------
    allpixels csv file, allfires geoparquet, fgb of individual fire with perimeters merged by day
    """
    from fireatlas import FireIO, preprocess

    # preprocessing steps to do at the start
    preprocess.preprocess_landcover()

    # read perimeter data + preprocess
    perimeter_gdf = gpd.read_file(perimeter_gdf_path)
    perimeter_gdf = FireIO.preprocess_polygon(perimeter_gdf, data_source=data_source, id_col=id_col)

    # if there's more than one perimeter, loop through each one
    # wrap fire_forward in a try/except block to catch any errors and continue processing
    def Fire_Forward_try(tst, ted, perimeter_series, sat=sat):
        try:
            constrainByShape_Fire_Forward(tst, ted, perimeter_series, sat)
        except Exception as e:
            logger.error(e)
            pass
    if len(perimeter_gdf) > 1:
        perimeter_gdf.apply(lambda x: Fire_Forward_try(tst, ted, x), axis=1)
    else:
        constrainByShape_Fire_Forward(tst, ted, perimeter_gdf, sat)


if __name__ == "__main__":
    """ The main code to run time forwarding for a time period
    """
    
    from dask.distributed import performance_report
    
    parser = argparse.ArgumentParser(description="registered MAAP.DPS jobs call ./run_dps.sh which delegates to this run function")
    parser.add_argument("run_function_name", help="The name of the function in ./FireRun.py to call")
    args = parser.parse_args()

    t1 = time.time()

    try:
        run_func = globals()[args.run_function_name]
        logger.info(f"[ RUNNING ]: {run_func}")
        #with performance_report(filename="dask-report.html"):
        run_func()
    except Exception as e:
        logger.exception(e)

    t2 = time.time()
    print(f"{(t2-t1)/60.} minutes used to run the whole code")


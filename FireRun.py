""" FireRun
Module to control different runs
"""


def Yearbatchrun(year, tst=None, ted=None, restart=False):
    """ Run the code for each single year
    """
    import FireMain, FireSummary, FireGdf_merge, FireGdf_sfs_merge, FireGdf_ign, FireGdf_final

    import time

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


def CreekSamplerun():
    """
    need to change firesrc in FireConsts.py 'SNPP'
    """
    import FireMain, FireGpkg, FireGpkg_sfs

    tst = (2020, 9, 5, "AM")
    ted = (2020, 11, 5, "PM")
    region = ("Creek", [-119.5, 36.8, -118.9, 37.7])

    # do fire tracking
    FireMain.Fire_Forward(tst=tst, ted=ted, restart=True, region=region)

    # calculate and save snapshot files
    FireGpkg.save_gdf_trng(tst=tst, ted=ted, regnm=region[0])

    # calculate and save single fire files
    FireGpkg_sfs.save_sfts_trng(tst, ted, regnm=region[0])


def CreekSamplerunNOAA20():
    """
    need to change firesrc in FireConsts.py 'NOAA20'
    """
    import FireMain, FireGpkg, FireGpkg_sfs

    tst = (2020, 9, 5, "AM")
    ted = (2020, 11, 5, "PM")
    region = ("CreekNOAA20", [-119.5, 36.8, -118.9, 37.7])

    # do fire tracking
    FireMain.Fire_Forward(tst=tst, ted=ted, restart=True, region=region)

    # calculate and save snapshot files
    FireGpkg.save_gdf_trng(tst=tst, ted=ted, regnm=region[0])

    # calculate and save single fire files
    FireGpkg_sfs.save_sfts_trng(tst, ted, regnm=region[0])


def CreekSamplerunVIIRS():
    """
    need to change firesrc in FireConsts.py 'VIIRS'
    """
    import FireMain, FireGpkg, FireGpkg_sfs

    tst = (2020, 9, 5, "AM")
    ted = (2020, 11, 5, "PM")
    region = ("CreekVIIRS", [-119.5, 36.8, -118.9, 37.7])

    # do fire tracking
    FireMain.Fire_Forward(tst=tst, ted=ted, restart=True, region=region)

    # calculate and save snapshot files
    FireGpkg.save_gdf_trng(tst=tst, ted=ted, regnm=region[0])

    # calculate and save single fire files
    FireGpkg_sfs.save_sfts_trng(tst, ted, regnm=region[0])


def CreekRegionSamplerun():
    import FireMain, FireGpkg, FireGpkg_sfs

    tst = (2020, 9, 1, "AM")
    ted = (2020, 9, 30, "PM")
    region = ("CreekRegion", [-120, 36, -118, 38])

    # do fire tracking
    FireMain.Fire_Forward(tst=tst, ted=ted, restart=True, region=region)

    # calculate and save snapshot files
    FireGpkg.save_gdf_trng(tst=tst, ted=ted, regnm=region[0])

    # calculate and save single fire files
    FireGpkg_sfs.save_sfts_trng(tst, ted, regnm=region[0])


def CArun():
    import FireIO, FireMain, FireGpkg, FireGpkg_sfs

    CAshp = FireIO.get_Cal_shp()
    region = ("California", CAshp)

    tst = (2020, 9, 1, "AM")
    ted = (2020, 9, 30, "PM")

    # do fire tracking
    FireMain.Fire_Forward(tst=tst, ted=ted, restart=True, region=region)

    # calculate and save snapshot files
    FireGpkg.save_gdf_trng(tst=tst, ted=ted, regnm=region[0])

    # calculate and save single fire files
    FireGpkg_sfs.save_sfts_trng(tst, ted, regnm=region[0])

def CArunNRT():
    import FireIO, FireMain, FireGpkg, FireGpkg_sfs, FireObj
    import FireConsts
    import DataCheckUpdate
    from datetime import datetime
    import os
    ctime = datetime.now()

    CAshp = FireIO.get_Cal_shp()
    region = ("CA", CAshp)

    lts = FireIO.get_lts_serialization(regnm=region[0])
    if lts == None:
        tst = [ctime.year, 1, 1, 'AM']
    else:
        tst = FireObj.t_nb(lts, nb="previous")

    if ctime.hour >= 18:
        ampm = 'PM'
    else:
        ampm = 'AM'
    ted = [ctime.year, ctime.month, ctime.day, ampm]
    print(f"Running code from {tst} to {ted}.")

    # Download data
    # Download Suomi-NPP
    DataCheckUpdate.update_VNP14IMGTDL()
    # Download NOAA-20
    DataCheckUpdate.update_VJ114IMGTDL()
    # Download GridMET
    # TODO ...
    
    FireMain.Fire_Forward(tst=tst, ted=ted, restart=False, region=region)
    FireGpkg.save_gdf_trng(tst=tst, ted=ted, regnm=region[0])
    FireGpkg_sfs.save_sfts_trng(tst, ted, regnm=region[0])

if __name__ == "__main__":
    """ The main code to run time forwarding for a time period
    """
    import sys

    sys.path.insert(
        1,
        "/Users/yangchen/GoogleDrive/My/My.Research/UCI/ProjectsFull/California.fire/Code/fireatlas",
    )

    # Yearbatchrun(2020)
    import time

    t1 = time.time()
    # CreekSamplerun()
    # CreekSamplerunNOAA20()
    CreekSamplerunVIIRS()

    # CreekRegionSamplerun()
    # CArun()
    t2 = time.time()
    print(f"{(t2-t1)/60.} minutes used to run the whole code code")

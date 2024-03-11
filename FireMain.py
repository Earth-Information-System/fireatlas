""" FireMain
Main module for running the fire object tracking along time

List of functions
-----------------
* Fobj_init: Initialize the fire object for a given time
* Fire_expand: Use daily new AF pixels to create new Fobj or combine with existing Fobj
* Fire_merge: For newly formed/expanded fire objects close to existing active fires, merge them
* Fire_Forward: The wrapper function to progressively track all fire events for a time period

Modules required
----------------
* Firelog
* FireObj
* FireIO
* FireClustering
* FireVector
* FireConsts
"""

# Use a logger to record console output
from FireLog import logger

from FireTypes import Region, TimeStep
from utils import timed


# Functions
def correct_nested_ids(mergetuple):
    """ correct the target fids after nested merging
    this is done before the merging happens in cases when several fires merge in one time step
    and also in the last time step to correct the heritage when several fires merge in different time steps

    Parameters
    ----------
    mergetuple: a list of tuples
        a list containing source and target ids for merging

    Returns
    -------
    mergetuple: a list of tuples
        a list containing source and corrected target ids for merging
    """
    import collections

    # 1)check if all keys (source fires) are unique
    src, tgt = zip(*mergetuple)
    tgt = list(tgt)
    count_keys = collections.Counter(src)
    not_unique = [key for key in count_keys if count_keys[key] > 1]
    # if not: replace with smallest tgt number
    if len(not_unique) > 0:
        for key in not_unique:
            tgt1 = min(
                [tgt[ind] for ind in range(len(tgt)) if src[ind] == key]
            )  # tgt with smallest fid
            indices = [ind for ind in range(len(tgt)) if src[ind] == key]
            for ind in indices:
                tgt[ind] = tgt1

    mergetuple = list(zip(src, tgt))

    # 2)correct nested ids
    mergedict = dict(mergetuple)
    src, tgt = zip(*mergetuple)
    while any(item in src for item in set(tgt)):
        for i, tgt0 in enumerate(tgt):
            if tgt0 in src:
                mergetuple[i] = (mergetuple[i][0], mergedict[tgt0])
        src, tgt = zip(*mergetuple)
    mergetuple = list(set(mergetuple))

    return mergetuple


def set_eafirerngs(allfires, fids):
    """ Return a list of fire connecting ranges from a list of fire ids
        fire connecting range is the hull plus a buffer (connectivity_fire)
    Parameters
    ----------
    allfires : Allfires object
        the input allfires object
    fids : list
        the list of fire ids

    Returns
    -------
    eafirerngs : list
        the list of fire connecting ranges corresponding to the sequence of fids
    """
    import FireFuncs

    # extract existing active fire data (use extending ranges)
    firerngs = []
    for fid in fids:
        f = allfires.fires[fid]  # fire
        CONNECTIVITY_FIRE_KM = FireFuncs.get_CONNECTIVITY_FIRE(f)
        rng = f.hull.buffer(CONNECTIVITY_FIRE_KM * 1000)
        firerngs.append(rng)
    return firerngs


def set_sleeperrngs(allfires, fids):
    """ Return a list of fire connecting ranges from a list of fire ids
        fire connecting range is the hull plus a buffer (connectivity_sleeper)
    Parameters
    ----------
    allfires : Allfires object
        the input allfires object
    fids : list
        the list of fire ids

    Returns
    -------
    eafirerngs : list
        the list of fire connecting ranges corresponding to the sequence of fids
    """
    import FireConsts

    # extract existing active fire data (use extending ranges)
    sleeperrngs = []
    for fid in fids:
        f = allfires.fires[fid]  # fire
        rng = f.hull.buffer(FireConsts.CONNECTIVITY_SLEEPER_KM * 1000)
        sleeperrngs.append(rng)
    return sleeperrngs


def maybe_remove_static_sources(region: Region, input_data_dir: str) -> Region:
    """ Modify region to exclude static sources

    Parameters
    ----------
    region : the run region. 

    Returns
    -------
    region : region obj
        A region that is the difference between the user-supplied region and the points identified as static flaring/gas according to the source. Creates a "swiss cheese"- like region, with negative space where there were points, with a buffer around points determined by "remove_static_sources_buffer". 
    """
    import os
    from FireIO import get_reg_shp
    from FireConsts import epsg, remove_static_sources_bool, remove_static_sources_sourcefile, remove_static_sources_buffer

    import geopandas as gpd
    import pandas as pd
    
    if not remove_static_sources_bool:
        return region
    
    # get source data geometry
    global_flaring = pd.read_csv(os.path.join(input_data_dir, 'static_sources', remove_static_sources_sourcefile))
    global_flaring = global_flaring.drop_duplicates()
    global_flaring = global_flaring[0:(len(global_flaring.id_key_2017) - 1)]

    global_flaring = gpd.GeoDataFrame(global_flaring, geometry=gpd.points_from_xy(global_flaring.Longitude, global_flaring.Latitude)) # Convert to point geometries
    global_flaring["buffer_geometry"] = global_flaring.buffer(remove_static_sources_buffer)
    global_flaring = global_flaring.set_geometry(col = "buffer_geometry")
    
    # get region geometry
    reg = get_reg_shp(region[1])
    reg_df = gpd.GeoDataFrame.from_dict({"name":[region[0]], "geometry":[reg]}) # Put geometry into dataframe for join
    
    # ensure everything is in the same projection
    global_flaring = global_flaring.set_crs("EPSG:" + str(epsg)) ## Translate to the user-input coordinate system
    reg_df = reg_df.set_crs("EPSG:" + str(epsg))
    
    # Take the difference of points and region
    diff = gpd.tools.overlay(reg_df, global_flaring, how='difference')
    
    region = (diff.name[0], diff.geometry[0])
    return region
    

@timed
def Fire_expand_rtree(allfires, allpixels, tpixels, fids_ea):
    """ Use daily new AF pixels to create new Fobj or combine with existing Fobj

    Parameters
    ----------
    allfires : Allfires obj
        the existing Allfires object for the time step
    afp : pandas.DataFrame
        preprocessed dataframe of new active fire pixels
    fids_ea : list
        fire ids of existing active fires at previous time step

    Returns
    -------
    allfires : Allfires obj
        updated Allfires object for the day with new formed/expanded fire objects
    """
    import pandas as pd
    import FireObj, FireClustering, FireTime
    from FireConsts import expand_only, firessr
    from FireVector import cal_hull

    # initializations
    idmax = max([*allfires.fires.keys(), 0])  # maximum id of existing fires
    fids_expanded = []  # a list of fire ids that is expanded at t
    fids_new = []  # a list of fire ids that is created at t

    # derive fire connecting ranges of existing active fires (fids_ea)
    eafirerngs = set_eafirerngs(allfires, fids_ea)

    # create a spatial index based on geometry bounds of fire connecting ranges
    ea_idx = FireClustering.build_rtree(eafirerngs)

    # loop over all new clusters (0:cid-1) and determine its fate
    FP2expand = {}  # a dict to record {fid : Firepixel objects} pairs
    for ic, pixels in tpixels.groupby("initial_cid"):
        hull = cal_hull(pixels[["x", "y"]].values, sensor=firessr)

        # if the cluster is close enough to an existing active fire object
        #   record all pixels to be added to the existing object (no actual changes on existing fire objects)
        
        # find potential neighbors using spatial index
        id_cfs = FireClustering.idx_intersection(ea_idx, hull.bounds)  
        clusterdone = False

         # loop over all potential eafires
        for id_cf in id_cfs:
            # each cluster can only be appended to one existing object
            if clusterdone == False:
                # determine if cluster touch fire connecting range
                if eafirerngs[id_cf].intersects(hull):
                    # record existing target fire id in fid_expand list
                    # this is the fire id of the existing active fire
                    fmid = fids_ea[id_cf]
                    # record pixels from target cluster (locs and time) along with the existing active fire object id
                    # single existing object, can have multiple new clusters to append
                    if fmid in FP2expand.keys():
                        FP2expand[fmid] = pd.concat([FP2expand[fmid], pixels])  # new pixels
                    else:
                        FP2expand[fmid] = pixels  # new pixels
                    
                    # record fmid to fid_expanded ? is this same as list(FP2expand.keys)?
                    fids_expanded.append(fmid)  
                    
                    # mark the cluster as done (no need to create new Fobj)
                    clusterdone = True

        # if this cluster can't be appended to any existing Fobj:
        #     create a new fire object using the new cluster
        # ignore creating new fires if expand_only is set to True
        if not expand_only:  
            # if the cluster is not added to existing active fires
            if clusterdone is False:  
                # create a new fire id and add it to the fid_new list
                id_newfire = idmax + 1
                fids_new.append(id_newfire)  # record id_newfire to fid_new

                # use the fire id and new fire pixels to create a new Fire object
                newfire = FireObj.Fire(id_newfire, allfires.t, allpixels, sensor=firessr)
                newfire.t_st = newfire.t
                newfire.pixels = pixels
                newfire.hull = hull
                newfire.updatefline()
                newfire.updateftype()  # update the fire type

                # add the new fire object to the fires list in the Allfires object
                allfires.fires[id_newfire] = newfire

                # increase the maximum id
                idmax += 1

    # update the expanded fire object (do the actual pixel appending and attributes changes)
    # fire attributes need to be manually changed:
    #  - end time; - pixels; - newpixels, - hull
    if len(FP2expand) > 0:
        for fmid, newpixels in FP2expand.items():
            # the target existing fire object
            f = allfires.fires[fmid]

            # update current time, end time
            f.t = allfires.t
            f.t_ed = allfires.t

            # extend pixels with newpixels
            f.pixels = pd.concat([f.pixels, newpixels])

            f.updatefhull()
            f.updatefline()

            # update the fire type
            f.updateftype()

    # remove duplicates and sort the fid_expanded
    fids_expanded = sorted(set(fids_expanded))

    # record fid change for expanded and new
    allfires.record_fids_change(fids_expanded=fids_expanded, fids_new=fids_new)

    return allfires


@timed
def Fire_merge_rtree(allfires, fids_ne, fids_ea, fids_sleep):
    """ For newly formed/expanded fires close to existing active fires or sleepers, merge them

    Parameters
    ----------
    allfires : Allfires obj
        the existing Allfires object for the time step
    fids_ne : list
        ids of newly formed/expanded fires
    fids_ea : list
        ids of existing active fire objects (including newly formed/expanded fires)

    Returns
    -------
    allfires : Allfires obj
        Allfires obj after fire merging
    """
    import pandas as pd
    import FireClustering

    # extract existing active fire data (use extending ranges)
    eafirerngs = set_eafirerngs(allfires, fids_ea)

    # create a spatial index based on geometry bounds of fire connecting ranges
    ea_idx = FireClustering.build_rtree(eafirerngs)

    # extract new and recently expanded fire data (use hulls without buffer)
    nefires = [allfires.fires[fid] for fid in fids_ne]
    nefirehulls = [f.hull for f in nefires]

    # loop over all fire objects that have newly expanded or formed, record merging fire id pairs
    fids_merge = []  # initialize the merged fire id pairs (source id:target id)
    # flag to mark an newly expanded fire obj that has been invalidated
    firedone = {i: False for i in fids_ne}  
    for id_ne in range(len(nefires)):
        fid_ne = fids_ne[id_ne]  # newly formed/expanded fire id
        if (
            firedone[fid_ne] == False
        ):  # skip objects that have been merged to others in earlier loop
            # potential neighbors
            id_cfs = FireClustering.idx_intersection(ea_idx, nefirehulls[id_ne].bounds)
            # loop over all potential neighbor fobj candidiates
            for id_ea in id_cfs:
                fid_ea = fids_ea[id_ea]  # fire id of existing active fire
                # if fid_ne == fid_ea, skip;
                # if the expanded fire has been merged to a existing active fire, skip the rest loops
                if fid_ne != fid_ea:
                    # if fire fmid is within distance of fire fid, two objects will merge
                    if nefirehulls[id_ne].intersects(eafirerngs[id_ea]):
                        # the fire id of neighboring active Fobj
                        # depending on which fid is smaller, merge the two fire objects in different directions
                        if fid_ea > fid_ne:  # merge fid_ea to fid_ne
                            fids_merge.append((fid_ea, fid_ne))
                            if fid_ea in firedone.keys():
                                firedone[
                                    fid_ea
                                ] = True  # remove fid_ea from the newly expanded fire list (since it has been invalidated)
                        else:  # merge fid_ne to fid_ea
                            fids_merge.append((fid_ne, fid_ea))
                            # fid_ne is merged to others, so stop it and check the next id_ne
                            ## technically the eafirerngs and nefirehulls have to be updated before the next loop happens
                            ## else it can happen that intersections are not being detected
                            ## if more than two fires grow together!!
                            # break

    # now check if any of the sleeper fires may have reactivated by new/expanded fires
    if len(fids_sleep) > 0:  # check if there are potential sleepers
        # extract existing sleeping fires and their firelines
        sleepfires = [allfires.fires[fid] for fid in fids_sleep]

        # sleepflines = [f.fline for f in sleepfires]; if no fline, use fline_prior
        sleepflines = [
            f.fline if f.fline is not None else f.fline_prior for f in sleepfires
        ]

        # extract ne fires sleeper range
        nefiresleeperrangs = set_sleeperrngs(allfires, fids_ne)

        # create a spatial index based on geometry bounds of ne fire sleeper ranges
        ne_idx = FireClustering.build_rtree(nefiresleeperrangs)

        # do the check analoguous to above; loop over each sleeper fire
        firedone = {
            i: False for i in fids_sleep
        }  # flag to mark an sleeper fire obj that has been invalidated
        for id_sleep in range(len(sleepfires)):
            fid_sleep = fids_sleep[id_sleep]  # sleeper fire id
            if (
                sleepflines[id_sleep] == None
            ):  # if there is no fire line (last active detection within), skip
                continue
            if (
                firedone[fid_sleep] == False
            ):  # skip objects that have been merged to others in earlier loop
                # potential neighbors
                id_cfs = FireClustering.idx_intersection(
                    ne_idx, sleepflines[id_sleep].bounds
                )
                # loop over all potential neighbour fobj candidates
                for id_ne in id_cfs:
                    fid_ne = fids_ne[id_ne]
                    if nefiresleeperrangs[id_ne].intersects(sleepflines[id_sleep]):
                        # depending on which fid is smaller, merge the two fire objects in different directions
                        if (
                            fid_ne > fid_sleep
                        ):  # merge new fire to sleeper, reactivate sleeper
                            fids_merge.append((fid_ne, fid_sleep))
                            if fid_ne in firedone.keys():
                                firedone[
                                    fid_ne
                                ] = True  # remove fid_ne from the newly expanded fire list (since it has been invalidated)
                        else:  # merge sleeper to new or expanded fire
                            fids_merge.append((fid_sleep, fid_ne))

    # loop over each pair in the fids_merge, and do modifications for both target and source objects
    #  - target: t_ed; pixels, newpixels, hull
    #  - source: invalidated
    if len(fids_merge) > 0:
        # fids_merge needs to be corrected if several fires merge at once!
        # i.e. if fire 2 merges into fire 1 and fire 3 merges into fire 2
        # in this case not correcting fids_merge will lead to invalidation of fire 3!!!
        fids_merge = correct_nested_ids(fids_merge)

        for fid1, fid2 in fids_merge:
            # update source and target objects
            f_source = allfires.fires[fid1]
            f_target = allfires.fires[fid2]

            # - target fire t_ed set to current time
            f_target.t = allfires.t
            f_target.t_ed = allfires.t

            # just in case: set target to valid (is this needed?)
            f_target.invalid = False

            # - target fire add source pixels to pixels and newpixels
            f_target.pixels = pd.concat([f_target.pixels, f_source.pixels])

            # - update the hull using previous hull and new pixels
            f_target.updatefhull()
            f_target.updatefline()

            # invalidate and deactivate source object
            f_source.invalid = True
            f_source.mergeid = f_target.mergeid

            # update target fire ftype
            f_target.updateftype()

            # record the heritages
            allfires.heritages.append((fid1, fid2))

        # remove duplicates and record fid change for merged and invalidated
        fids_invalid, fids_merged = zip(*fids_merge)
        fids_merged = sorted(set(fids_merged))
        fids_invalid = sorted(set(fids_invalid))
        allfires.record_fids_change(fids_merged=fids_merged, fids_invalid=fids_invalid)

    return allfires

@timed
def Fire_Forward_one_step(allfires, allpixels, t, region):
    import FireConsts, FireTime
    
    logger.info("--------------------")
    logger.info(f"Fire tracking at {t}")

    if FireTime.isyearst(t):
        allfires.newyear_reset(region[0])

    # 1. record existing active fire ids (before fire tracking at t)
    fids_ea = allfires.fids_active

    # 2. update t of allfires, clean up allfires and fire object
    allfires.cleanup(t)

    tpixels = allpixels[allpixels["t"] == FireTime.t2dt(t)]

    # 4.5. if active fire pixels are detected, do fire expansion/merging
    if len(tpixels) > 0:
        # 4. do fire expansion/creation using allpixels
        allfires = Fire_expand_rtree(allfires, allpixels, tpixels, fids_ea)

        # 5. do fire merging using updated fids_ne, fid_ea, fid_sleep
        fids_ne = allfires.fids_ne  # new or expanded fires id
        fids_ea = sorted(set(fids_ea + allfires.fids_new))  # existing active fires (new fires included)
        fids_sleep = allfires.fids_sleeper
        if len(fids_ne) > 0:
            allfires = Fire_merge_rtree(allfires, fids_ne, fids_ea, fids_sleep)

    # 7. manualy invalidate static fires (with exceptionally large fire density)
    if FireConsts.opt_rmstatfire:
        allfires.invalidate_statfires()

    # 8. log changes
    #  - record fid_updated (the fid of fires that change in the time step) to allfires object and logger
    logger.info(f"fids_expand: {len(allfires.fids_expanded)}")
    logger.info(f"fids_new: {len(allfires.fids_new)}")
    logger.info(f"fids_merged: {len(allfires.fids_merged)}")
    logger.info(f"fids_invalid: {len(allfires.fids_invalid)}")

    # 9. correct heritages at each time step?
    if len(allfires.heritages) > 0:
        allfires.heritages = correct_nested_ids(allfires.heritages)

    # 10. update allfires gdf
    allfires.update_gdf()

    return allfires

@timed
def Fire_Forward(tst: TimeStep, ted: TimeStep, sat=None, restart=False, region=None):
    """ The wrapper function to progressively track all fire events for a time period

    Parameters
    ----------
    tst : tuple, (int,int,int,str)
        the year, month, day and 'AM'|'PM' at start time
    ted : tuple, (int,int,int,str)
    the year, month, day and 'AM'|'PM' at end time
    sat : str, 'SNPP', 'NOAA20', 'VIIRS', 'BAMOD
        if set, overrides `FireConsts.firesrc`
    restart : bool
        if set to true, force to initiate an object

    Returns
    -------
    allfires : FireObj allfires object
        the allfires object at end date
    """
    import FireTime, FireObj, FireConsts
    import preprocess
    import pandas as pd

    if sat is None:
        sat = FireConsts.firesrc

    # initialize allfires object
    allfires = FireObj.Allfires(tst)

    # read in all preprocessed pixel data
    allpixels = pd.concat([
        preprocess.read_preprocessed(t, sat=sat, region=region)
        for t in FireTime.t_generator(tst, ted)
    ])
    allpixels["fid"] = -1
    allpixels["in_fline"] = None    

    # loop over every t during the period and mutate allfires, allpixels
    for t in FireTime.t_generator(tst, ted):
        allfires = Fire_Forward_one_step(allfires, allpixels, t, region)

    return allfires, allpixels


if __name__ == "__main__":
    """ The main code to run time forwarding for a time period
    """

    import time
    import FireGpkg
    import FireGpkg_sfs

    t1 = time.time()
    print(t1)

    # # set the start and end time
    # tst=(2021,7,13,'AM')
    # ted=(2021,9,15,'PM')
    # region = ('Dixie',[-121.6,39.8,-120.5,40.6])

    tst = (2020, 9, 25, 'PM')
    ted = (2020, 12, 31, "PM")
    #print(tst)
    #region = ("Creek", [-119.5, 36.8, -118.9, 37.7])
    #region = ('CA',[-124.409591, 32.534155999999996, -114.131211, 42.009518])
    #region = ('Thomas',[-119.79311723013566,34.162521752180936,-118.87850541372941,34.791948775281746])
    #region = ('Meyers',[-113.77550545807375,45.84172304592036,-113.426689540105,46.13941078007829])
    #region = ('Bighorn',[-111.2483396974153,32.107921771038576,-110.2980223146028,32.73852603996812])
    #region = ('Frye',[-110.02849137185659,32.568462386228475,-109.69752823709096,32.84117803581184])
    #region = ('Fish',[-117.98761271684656,34.14492745779149,-117.9021253511239,34.21536501245526])
    #region = ('GrizzlyCreek',[-107.32359115976944,39.51527120096794,-107.04481308359756,39.698839413262284])
    region = ('WesternUS_REDO',[-125.698046875,31.676476158707615,-101.00078125,49.51429477264348])
    #region = ('HermitsPeakREDO',[-105.62745646198083,35.373429737675505,-105.18251017291833,36.26028722617026])
    # region = ('CONUS',[-126.401171875,-61.36210937500001,24.071240929282325,49.40003415463647])
    #region = ('Caldor',[-120.69305873258455,38.52600288552201,-119.90341639860017,38.916006495378696])
    # Run the time forward and record daily fire objects .pkl data and fire attributes .GeoJSON data
    print("----------------------------------------")
    print("Running Fire_Forward")
    print("----------------------------------------")
    #Fire_Forward(tst=tst, ted=ted, restart=True, region=region)

    # calculate and save snapshot files
    print("----------------------------------------")
    print("Running save_gdf_trng")
    print("----------------------------------------")
    #FireGpkg.save_gdf_trng(tst=tst, ted=ted, regnm=region[0])

    # calculate and save single fire files
    print("----------------------------------------")
    print("Running save_sfts_trng")
    print("----------------------------------------")
    FireGpkg_sfs.save_sfts_trng(tst, ted, regnm=region[0])

    t2 = time.time()
    print(f"{(t2-t1)/60.} minutes used to run code")

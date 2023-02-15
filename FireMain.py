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
    import FireFuncs, FireVector

    # extract existing active fire data (use extending ranges)
    firerngs = []
    for fid in fids:
        f = allfires.fires[fid]  # fire
        CONNECTIVITY_FIRE_KM = FireFuncs.get_CONNECTIVITY_FIRE(f)
        rng = FireVector.addbuffer(f.hull, CONNECTIVITY_FIRE_KM * 1000)
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
    import FireFuncs, FireVector

    # extract existing active fire data (use extending ranges)
    sleeperrngs = []
    for fid in fids:
        f = allfires.fires[fid]  # fire
        CONNECTIVITY_SLEEPER_KM = FireFuncs.get_CONNECTIVITY_SLEEPER()
        rng = FireVector.addbuffer(f.hull, CONNECTIVITY_SLEEPER_KM * 1000)
        sleeperrngs.append(rng)
    return sleeperrngs


def Fobj_init(tst, regnm, restart=False):
    """ Initialize the fire object for a given time. This can be from the object
    saved at previous time, or can be initialized using Allfires().

    Parameters
    ----------
    tst : tuple, (int,int,int,str)
        the year, month, day and 'AM'|'PM' during the intialization
    restart : bool
        if set to true, force to initiate an object

    Returns
    -------
    allfires : Allfires obj
        the fire object for the previous time step
    """
    import FireObj, FireIO, FireTime

    pst = FireTime.t_nb(tst, nb="previous")  # previous time step
    if FireIO.check_fobj(pst, regnm, activeonly=False) & (restart == False):
        allfires = FireIO.load_fobj(pst, regnm, activeonly=False) # load all fires (including dead)
        allfires.cleanup(tst)  # update time and reset lists
        # # if it's the first time step of a calendar year, reset all fires id
        # if (tst[1]==1 & tst[2]==1 & tst[3]=='AM'):
        #     allfires.newyear_reset()
    else:  # if no pkl file at previous time step or restart is set to True
        allfires = FireObj.Allfires(tst)

    return allfires


def remove_static_sources(region, source):
    """ Modify region to exclude static sources

    Parameters
    ----------
    region : the run region. 
    source : str
        Name of the file in dirextdata that containts statics sources with a Longitude and Latitude column.  

    Returns
    -------
    region : region obj
        A region that is the diference between the user-supplied region and the points identified as static flaring/gas according to the source. Creates a "swiss cheese"- like region, with negative space where there were points, with a buffer around points determined by "remove_static_sources_buffer". 
    """
    # import
    import os
    from FireIO import get_reg_shp
    from FireConsts import dirextdata, epsg, remove_static_sources_buffer
    import shapely
    import shapely.geometry
    from shapely.geometry import Point, Polygon
    import geopandas as gpd
    
    # get source data geometry
    global_flaring = gpd.read_file(os.path.join(dirextdata,'static_sources', source))
    global_flaring = global_flaring.drop_duplicates()
    global_flaring = global_flaring[0:(len(global_flaring.id_key_2017) - 1)]

    global_flaring = gpd.GeoDataFrame(global_flaring, geometry=gpd.points_from_xy(global_flaring.Longitude, global_flaring.Latitude)) # Convert to point geometries
    global_flaring["buffer_geometry"] = global_flaring.buffer(remove_static_sources_buffer)
    global_flaring = global_flaring.set_geometry(col = "buffer_geometry")
    
    # get region geometry
    reg = get_reg_shp(region[1])
    reg_df = gpd.GeoDataFrame.from_dict({"name":[region[0]], "geometry":[reg]}) # Put geometry into dataframe for join
    
    # ensure everything is in the same projection
    global_flaring.set_crs("EPSG:4346") ## Set a lat lon system
    global_flaring.to_crs("EPSG:" + epsg) ## Translate to the user-input coordinate system
    reg_df.set_crs("EPSG:" + epsg)
    
    # Take the difference of points and region
    diff = gpd.tools.overlay(reg_df, global_flaring, how='difference')
    
    region = (diff.name[0], diff.geometry[0])
    return(region)
    

def Fire_expand_rtree(allfires, afp, fids_ea, log=True):
    """ Use daily new AF pixels to create new Fobj or combine with existing Fobj

    Parameters
    ----------
    allfires : Allfires obj
        the existing Allfires object for the time step
    afp : 5-element list
        (lat, lon, line, sample, FRP) of new active fire pixels
    fids_ea : list
        fire ids of existing active fires at previous time step

    Returns
    -------
    allfires : Allfires obj
        updated Allfires object for the day with new formed/expanded fire objects
    """
    # import time
    import FireObj, FireClustering, FireVector, FireFuncs
    from FireConsts import expand_only, firessr

    # initializations
    idmax = (
        allfires.number_of_fires - 1
    )  # maximum id of existing fires (max(allfires.fires.keys())?)
    fids_expanded = []  # a list of fire ids that is expanded at t
    fids_new = []  # a list of fire ids that is created at t

    # derive fire connecting ranges of existing active fires (fids_ea)
    eafirerngs = set_eafirerngs(allfires, fids_ea)

    # create a spatial index based on geometry bounds of fire connecting ranges
    ea_idx = FireClustering.build_rtree(eafirerngs)

    # do preliminary clustering using new active fire locations (assign cid to each pixel)
    afp_loc = list(zip(afp.x, afp.y))
    CONNECTIVITY_CLUSTER = FireFuncs.get_CONNECTIVITY_CLUSTER()
    cid = FireClustering.do_clustering(
        afp_loc, CONNECTIVITY_CLUSTER
    )  # cluster id for each afp_loc
    if log:
        logger.info(f"New fire clusters of {max(cid)} at this time step")

    # loop over all new clusters (0:cid-1) and determine its fate
    FP2expand = {}  # a diction to record {fid : Firepixel objects} pairs
    for ic in range(max(cid) + 1):
        # create cluster object using all newly detected active fires within a cluster
        #   (-1 means no source fireid)
        pixels = [
            FireObj.FirePixel(
                afp.iloc[i].x,
                afp.iloc[i].y,
                afp.iloc[i].Lon,
                afp.iloc[i].Lat,
                afp.iloc[i].FRP,
                afp.iloc[i].DS,
                afp.iloc[i].DT,
                afp.iloc[i].ampm,
                afp.iloc[i].YYYYMMDD_HHMM,
                afp.iloc[i].Sat,
                -1,
            )
            for i, v in enumerate(cid)
            if v == ic
        ]  # pixels
        cluster = FireObj.Cluster(
            ic, pixels, allfires.t, sensor=firessr
        )  # form cluster
        hull = cluster.hull  # save hull to reduce computational cost

        # if the cluster is close enough to an existing active fire object
        #   record all pixels to be added to the existing object (no actuall changes on existing fire objects)
        id_cfs = FireClustering.idx_intersection(
            ea_idx, cluster.b_box
        )  # potential neighbours using spatial index
        clusterdone = False
        for id_cf in id_cfs:  # loop over all potential eafires
            if (
                clusterdone == False
            ):  # one cluster can only be appended to one existing object
                if eafirerngs[id_cf].intersects(
                    hull
                ):  # determine if cluster touch fire connecting range
                    # record existing target fire id in fid_expand list
                    fmid = fids_ea[
                        id_cf
                    ]  # this is the fire id of the existing active fire
                    # record pixels from target cluster (locs and time) along with the existing active fire object id
                    # newFPs = [FireObj.FirePixel(p.x,p.y,p.lon,p.lat,p.frp,p.DS,p.DT,p.ampm,p.datetime,p.sat,fmid) for p in pixels] # new FirePixels from the cluster
                    if (
                        fmid in FP2expand.keys()
                    ):  # single existing object, can have multiple new clusters to append
                        FP2expand[fmid] = FP2expand[fmid] + pixels  # newFPs
                    else:
                        FP2expand[fmid] = pixels  # newFPs
                    fids_expanded.append(
                        fmid
                    )  # record fmid to fid_expanded ? is this same as list(FP2expand.keys)?
                    clusterdone = (
                        True  # mark the cluster as done (no need to create new Fobj)
                    )

        # if this cluster can't be appended to any existing Fobj, create a new fire object using the new cluster
        if not expand_only:  # ignore creating new fires if expand_only is set to True
            if (
                clusterdone is False
            ):  # if the cluster is not added to existing active fires
                # create a new fire id and add it to the fid_new list
                id_newfire = idmax + 1
                fids_new.append(id_newfire)  # record id_newfire to fid_new

                # use the fire id and new fire pixels to create a new Fire object
                newfire = FireObj.Fire(id_newfire, allfires.t, pixels, sensor=firessr)
                newfire.updateftype()  # update the fire type

                # add the new fire object to the fires list in the Allfires object
                allfires.fires[id_newfire] = newfire

                # increase the maximum id
                idmax += 1

    # update the expanded fire object (do the actual pixel appending and attributes changes)
    # fire attributes need to be manualy changed:
    #  - end time; - pixels; - newpixels, - hull, - extpixels
    if len(FP2expand) > 0:
        for fmid, newFPs in FP2expand.items():
            # the target existing fire object
            f = allfires.fires[fmid]

            # update end time
            f.t_ed = allfires.t

            # update pixels
            f.pixels = f.pixels + newFPs
            f.newpixels = newFPs
            # if len(newFPs) > 0:
            #     f.actpixels = newFPs

            # update the hull using previous hull and previous exterior pixels
            # phull = f.hull   # previous hull
            pextlocs = [p.loc for p in f.extpixels]  # previous external pixels
            newlocs = [p.loc for p in newFPs]  # new added pixels
            # f.hull = FireVector.update_hull(phull,pextlocs+newlocs)  # use update_hull function to save time
            f.updatefhull(pextlocs + newlocs)

            # update exterior pixels
            # f.updateextpixels(f.extpixels+newFPs)
            f.updateextpixels(newFPs)
            # f.extpixels = FireVector.cal_extpixels(f.extpixels+newFPs,f.hull)

            f.updateftype()  # update the fire type
            # t1 = time.time()
            # logger.info(f'Update external pixels: {t1-t2}')

    # remove duplicates and sort the fid_expanded
    fids_expanded = sorted(set(fids_expanded))

    # record fid change for expanded and new
    allfires.record_fids_change(fids_expanded=fids_expanded, fids_new=fids_new)

    return allfires


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

    import FireClustering, FireVector, FireFuncs
    from FireConsts import firessr  # CONNECTIVITY_THRESHOLD_KM

    # extract existing active fire data (use extending ranges)
    eafirerngs = set_eafirerngs(allfires, fids_ea)

    # create a spatial index based on geometry bounds of fire connecting ranges
    ea_idx = FireClustering.build_rtree(eafirerngs)

    # extract new and recently expanded fire data (use hulls without buffer)
    nefires = [allfires.fires[fid] for fid in fids_ne]
    nefirehulls = [f.hull for f in nefires]

    # loop over all fire objects that have newly expanded or formed, record merging fire id pairs
    fids_merge = []  # initialize the merged fire id pairs (source id:target id)
    firedone = {
        i: False for i in fids_ne
    }  # flag to mark an newly expanded fire obj that has been invalidated
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

        # nefirebuf  = [FireVector.addbuffer(hull,sleeperthresh*1000) for hull in nefirehulls]
        # ne_idx = FireClustering.build_rtree(nefirebuf)

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
    #  - target: t_ed; pixels, newpixels, hull, extpixels
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
            f_target.t_ed = allfires.t

            # just in case: set target to valid (is this needed?)
            f_target.invalid = False

            # - target fire add source pixels to pixels and newpixels
            f_target.pixels = f_target.pixels + f_source.pixels
            f_target.newpixels = f_target.newpixels + f_source.newpixels

            # - update the hull using previous hull and previous exterior pixels
            phull = f_target.hull
            pextlocs = [p.loc for p in f_target.extpixels]
            newlocs = [p.loc for p in f_source.pixels]
            # f_target.hull = FireVector.update_hull(phull,pextlocs+newlocs, sensor=firessr)
            f_target.updatefhull(pextlocs + newlocs)

            # - use the updated hull to update exterior pixels
            f_target.extpixels = FireVector.cal_extpixels(
                f_target.extpixels + f_source.pixels, f_target.hull
            )

            # invalidate and deactivate source object
            f_source.invalid = True
            f_source.mergeid = f_target.mergeid

            # update target fire ftype
            f_target.updateftype()

            # record the heritages
            # allfires.heritages.append((fid1,fid2,allfires.t))
            allfires.heritages.append((fid1, fid2))

        # remove duplicates and record fid change for merged and invalidated
        fids_invalid, fids_merged = zip(*fids_merge)
        fids_merged = sorted(set(fids_merged))
        fids_invalid = sorted(set(fids_invalid))
        allfires.record_fids_change(fids_merged=fids_merged, fids_invalid=fids_invalid)

    return allfires


def Fire_Forward(tst, ted, restart=False, region=None):
    """ The wrapper function to progressively track all fire events for a time period
           and save fire object to pkl file and gpd to geojson files

    Parameters
    ----------
    tst : tuple, (int,int,int,str)
        the year, month, day and 'AM'|'PM' at start time
    ted : tuple, (int,int,int,str)
        the year, month, day and 'AM'|'PM' at end time
    restart : bool
        if set to true, force to initiate an object

    Returns
    -------
    allfires : FireObj allfires object
        the allfires object at end date
    """

    # import libraries
    import FireObj, FireIO, FireTime
    from FireConsts import firesrc, firenrt, opt_rmstatfire, remove_static_sources_bool, remove_static_sources_sourcefile

    import os
    import glob

    # used to record time of script running
    import time

    t1 = time.time()
    t0 = t1

    # initialize allfires object
    allfires = Fobj_init(tst, region[0], restart=restart)
    
    # remove static sources
    if remove_static_sources_bool: 
        region = remove_static_sources(region, remove_static_sources_sourcefile)

    # loop over all days during the period
    endloop = False  # flag to control the ending of the loop
    t = list(tst)  # t is the time (year,month,day,ampm) for each step
    while endloop == False:
        logger.info("")
        logger.info(t)
        print("Fire tracking at", t)

        if FireTime.isyearst(t):
            allfires.newyear_reset(region[0])

        # 1. record existing active fire ids (before fire tracking at t)
        fids_ea = allfires.fids_active

        # 2. update t of allfires, clean up allfires and fire object
        allfires.cleanup(t)

        # 3. read active fire pixels from VIIRS dataset
        i = 0
        while i < 5:
            try: 
                afp = FireIO.read_AFP(t, src=firesrc, nrt=firenrt, region=region)
                break
            except Exception as e:
                print(f"Attempt {i}/5 failed.")
                print(e)
                i += 1
                if not i < 5:
                    raise e
        
        # 4.5. if active fire pixels are detected, do fire expansion/merging
        if len(afp) > 0:
            # 4. do fire expansion/creation using afp
            t_expand = time.time()
            allfires = Fire_expand_rtree(allfires, afp, fids_ea)
            t_expand2 = time.time()
            logger.info(f"expanding fires {(t_expand2-t_expand)}")

            # 5. do fire merging using updated fids_ne, fid_ea, fid_sleep
            fids_ne = allfires.fids_ne  # new or expanded fires id
            fids_ea = sorted(
                set(fids_ea + allfires.fids_new)
            )  # existing active fires (new fires included)
            fids_sleep = allfires.fids_sleeper
            t_merge = time.time()
            if len(fids_ne) > 0:
                allfires = Fire_merge_rtree(allfires, fids_ne, fids_ea, fids_sleep)
            t_merge2 = time.time()
            logger.info(f"merging fires {(t_merge2-t_merge)}")

        # # 6. determine or update fire type
        # allfires.updateftypes()

        # 7. manualy invalidate static fires (with exceptionally large fire density)
        if opt_rmstatfire:
            allfires.invalidate_statfires()

        # 8. log and save
        #  - record fid_updated (the fid of fires that change in the time step) to allfires object and logger
        logger.info(f"fids_expand: {allfires.fids_expanded}")
        logger.info(f"fids_new: {allfires.fids_new}")
        logger.info(f"fids_merged: {allfires.fids_merged}")
        logger.info(f"fids_invalid: {allfires.fids_invalid}")

        # correct heritages at each time step?
        if len(allfires.heritages) > 0:
            allfires.heritages = correct_nested_ids(allfires.heritages)

        # 9. loop control
        #  - if t reaches ted, set endloop to True to stop the next loop
        if FireTime.t_dif(t, ted) == 0:
            endloop = True
            # # correct fire heritage of final time step
            # if len(allfires.heritages) > 0:
            #     allfires.heritages = correct_nested_ids(allfires.heritages)

        # 10. fire object save
        FireIO.save_fobj(allfires, t, region[0], activeonly=False)
        FireIO.save_fobj(allfires, t, region[0], activeonly=True)
        # if FireTime.t_dif(t,ted)==0:
        #     FireIO.save_fobj(allfires,t,region[0],activeonly=False)

        # 11. update t with the next time stamp
        #  - record running times for the loop
        t2 = time.time()
        logger.info(f"{(t2-t1)/60.} minutes used to run alg {t}")
        t1 = t2

        # - move to next time step
        t = FireTime.t_nb(t, nb="next")

    # record total running time
    t3 = time.time()
    logger.info(f"This running takes {(t3-t0)/60.} minutes")

    return allfires


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

    tst = (2016, 1, 1, 'AM')
    ted = (2016, 12, 31, "PM")
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
    #region = ('HermitsPeak',[-105.62745646198083,35.373429737675505,-105.18251017291833,36.26028722617026])
    # region = ('CONUS',[-126.401171875,-61.36210937500001,24.071240929282325,49.40003415463647])
    #region = ('Caldor',[-120.69305873258455,38.52600288552201,-119.90341639860017,38.916006495378696])
    # Run the time forward and record daily fire objects .pkl data and fire attributes .GeoJSON data
    print("----------------------------------------")
    print("Running Fire_Forward")
    print("----------------------------------------")
    Fire_Forward(tst=tst, ted=ted, restart=True, region=region)

    # calculate and save snapshot files
    print("----------------------------------------")
    print("Running save_gdf_trng")
    print("----------------------------------------")
    FireGpkg.save_gdf_trng(tst=tst, ted=ted, regnm=region[0])

    # calculate and save single fire files
    print("----------------------------------------")
    print("Running save_sfts_trng")
    print("----------------------------------------")
    FireGpkg_sfs.save_sfts_trng(tst, ted, regnm=region[0])

    t2 = time.time()
    print(f"{(t2-t1)/60.} minutes used to run code")

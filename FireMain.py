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
def get_CONNECTIVITY_THRESHOLD(fire):
    ''' get the CONNECTIVITY_THRESHOLD_KM value for a given fire
    '''
    from FireConsts import CONT_opt

    if CONT_opt == 0:      # use preset values
        from FireConsts import CONT_preset
        return CONT_preset
    elif CONT_opt == 1:    # use values from Fire Atlas data
        from FireConsts import CONT_CA
        return CONT_CA[fire.ftype]
    elif CONT_opt == 2:   # global universal run
        from FireConsts import FTYP_Glb
        ftype = FTYP_Glb[fire.ftype]
        fnpix = fire.n_pixels
        if  ftype == 'Temp Forest':
            v = fnpix * 2 /25 + 0.8    # 0.8 (0) - 2.8 (25)
        elif ftype == 'Trop Forest':
            v = fnpix * 0.7 /25 + 0.7    # 0.7 (0) - 1.4 (25)
        elif ftype == 'Bore Forest':
            v = fnpix * 3.2 /25 + 1.0    # 1.0 (0) - 4.2 (25)
        elif ftype == 'Savana':
            v = fnpix * 2.5 /25 + 0.9    # 0.9 (0) - 3.4 (25)
        else:
            v = 0.7
        return v


def correct_nested_ids(mergetuple):
    ''' correct the target fids after nested merging
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
    '''
    import collections

    # 1)check if all keys are unique
    src, tgt = zip(*mergetuple)
    tgt = list(tgt)
    count_keys = collections.Counter(src)
    not_unique = [key for key in count_keys if count_keys[key] > 1]
    # if not: replace with smallest tgt number
    if len(not_unique)>0:
        for key in not_unique:
            tgt1 = min([tgt[ind] for ind in range(len(tgt)) if src[ind] == key])
            indices = [ind for ind in range(len(tgt)) if src[ind] == key]
            for ind in indices:
                tgt[ind] = tgt1
    mergetuple = list(zip(src,tgt))

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


def Fobj_init(tst,regnm,restart=False):
    ''' Initialize the fire object for a given time. This can be from the object
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
    '''
    import FireObj, FireIO

    # if option set to restart, do Allfires objection initialiation
    if restart == True:
        allfires = FireObj.Allfires(tst)
    else:
        # if restart==False, and previous time step data are available, read from file
        pst = FireObj.t_nb(tst,nb='previous') # previous time step
        if FireIO.check_fobj(pst,regnm):
            allfires = FireIO.load_fobj(pst,regnm,activeonly=False)
            # if it's the first time step of a calendar year, reset all fires id
            if (tst[1]==1 & tst[2]==1 & tst[3]=='AM'):
                # allfires = NewYearReset(allfires)
                allfires.newyear_reset()
        else:
            allfires = FireObj.Allfires(tst)

    # # Initialize allfires using previous time Fobj value in a file
    # if (FireIO.check_fobj(pst) & (restart==False)):
    #     allfires = FireIO.load_fobj(pst)
    #
    # # If previous time value is unavailable, initialize an empty Fobj
    # else:
    #     allfires = FireObj.Allfires(tst)

    return allfires

def Fire_expand_rtree(allfires,afp,fids_ea,sensor='viirs',log=True):
    ''' Use daily new AF pixels to create new Fobj or combine with existing Fobj

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
    '''
    # import time
    import FireObj,FireClustering,FireVector
    from FireConsts import SPATIAL_THRESHOLD_KM, expand_only #,CONNECTIVITY_THRESHOLD_KM

    # t0 = time.time()
    # record current time for later use (t in allfires has been updated in the Fire_Forward function)
    t = allfires.t

    # do some initializations
    idmax = allfires.number_of_fires-1  # maximum id of existing fires
    fids_expanded = []      # a list recording Fobj ids which has been expanded
    fids_new = []           # a list recording ids of new Fobjs

    # expanding ranges of existing active fires (extracted using fids_ea)
    eafires     = [allfires.fires[fid]  for fid in fids_ea]

    eafirerngs = []
    for f in eafires:
        CONNECTIVITY_THRESHOLD_KM = get_CONNECTIVITY_THRESHOLD(f)
        eafirerngs.append(FireVector.addbuffer(f.hull,CONNECTIVITY_THRESHOLD_KM*1000))
    # eafirerngs  = [FireVector.addbuffer(f.hull,CONNECTIVITY_THRESHOLD_KM[f.ftype]*1000) for f in eafires]

    # inserting boxes into a spatial index
    ea_idx = FireClustering.build_rtree(eafirerngs)
    # t1 = time.time()
    # logger.info(f'Time to build rtree: {t1-t0}')

    # do preliminary clustering using new active fire locations (assign cid to each pixel)
    # afp_loc = [(x,y) for x,y,a,b,c in afp]
    afp_loc = list(zip(afp.Lat, afp.Lon))

    cid = FireClustering.do_clustering(afp_loc,SPATIAL_THRESHOLD_KM)  # this is the cluster id (starting from 0)
    if log:
        logger.info(f'New fire clusters of {max(cid)} at this time step')
    # t1 = time.time()
    # logger.info(f'Time for initial clustering: {t1-t0}')

    # loop over each of the new clusters (0:cid-1) and determine its fate
    FP2expand = {}  # the diction used to record fire pixels assigned to existing active fires {fid:Firepixels}
    for ic in range(max(cid)+1):
        # create cluster object using all newly detected active fires within a cluster
        # pixels = [afp[i] for i, v in enumerate(cid) if v==ic]
        pixels = [FireObj.FirePixel(afp.iloc[i].Lat,afp.iloc[i].Lon,afp.iloc[i].FRP,t,-1) for i, v in enumerate(cid) if v==ic]
        cluster = FireObj.Cluster(ic,pixels,t,sensor=sensor)  # create a Cluster object using the pixel locations
        hull = cluster.hull  # the hull of the cluster

        # extract potential neighbours using spatial index (used for prefilter)
        id_cfs = FireClustering.idx_intersection(ea_idx, cluster.b_box)

        # now check if the cluster is truely close to an existing active fire object
        # if yes, record all pixels to be added to the existing object
        clusterdone = False
        for id_cf in id_cfs:  # loop over all potential eafires
            if clusterdone == False:  # one cluster can only be appended to one existing object
                # if cluster touch the extending range of an existing fire
                if eafirerngs[id_cf].intersects(hull):
                    # record existing target fire id in fid_expand list
                    fmid = fids_ea[id_cf]  # this is the fire id of the existing active fire
                    # record pixels from target cluster (locs and time) along with the existing active fire object id
                    newFPs = [FireObj.FirePixel(p.lat,p.lon,p.frp,t,fmid) for p in pixels] # new FirePixels from the cluster
                    if fmid in FP2expand.keys():   # for a single existing object, there can be multiple new clusters to append
                        FP2expand[fmid] = FP2expand[fmid] + newFPs
                    else:
                        FP2expand[fmid] = newFPs

                    # logger.info(f'Fire {fmid} expanded with pixels from new cluster {ic}')

                    fids_expanded.append(fmid) # record fmid to fid_expanded

                    clusterdone = True   # mark the cluster as done (no need to create new Fobj)


        # if this cluster can't be appended to any existing Fobj, create a new fire object using the new cluster
        if not expand_only:
            #print('creating new fires')
            if clusterdone is False:
                # create a new fire id and add it to the fid_new list
                id_newfire = idmax + 1
                # logger.info(f'Fire {id_newfire} created with pixels from new cluster {ic}')
                fids_new.append(id_newfire)  # record id_newfire to fid_new

                # use the fire id and new fire pixels to create a new Fire object
                newfire = FireObj.Fire(id_newfire,t,pixels,sensor=sensor)

                # set the fire type for the newly created fire
                newfire.set_ftype()

                # add the new fire object to the fires list in the Allfires object
                allfires.fires.append(newfire)

                # increase the maximum id
                idmax += 1
    # t1 = time.time()
    # logger.info(f'Time for intersections: {t1-t0}')
    # update the expanded fire object (do the actual pixel appending)
    #  - fire attributes to change: end time; pixels; newpixels, hull, extpixels
    if len(FP2expand) > 0:
        for fmid, newFPs in FP2expand.items():

            # the target existing fire object
            f = allfires.fires[fmid]

            # update end time
            f.t_ed = t

            # update pixels
            f.pixels = f.pixels + newFPs
            f.newpixels = newFPs
            if len(newFPs) > 0:
                f.actpixels = newFPs

            # update the hull using previous hull and previous exterior pixels
            phull = f.hull
            pextlocs = [p.loc for p in f.extpixels]
            newlocs = [p.loc for p in newFPs]
            # t2 = time.time()
            f.hull = FireVector.update_hull(phull,pextlocs+newlocs)  # use update_hull function to save time
            # t1 = time.time()
            # logger.info(f'Updated hull: {t1-t2}')
            # use the updated hull to update exterior pixels
            f.extpixels = FireVector.cal_extpixels(f.extpixels+newFPs,f.hull)
            # t1 = time.time()
            # logger.info(f'Update external pixels: {t1-t2}')

    # remove duplicates and sort the fid_expanded
    fids_expanded = sorted(set(fids_expanded))

    # record fid change for expanded and new
    allfires.record_fids_change(fids_expanded=fids_expanded, fids_new=fids_new)

    # logger.info(f'In total, {len(fids_expanded)} fires expanded, and {len(fids_new)} fires created')

    return allfires

def Fire_merge_rtree(allfires,fids_ne,fids_ea,fids_sleep,sensor ='viirs'):
    ''' For newly formed/expanded fire objects close to existing active fires, merge them

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
    '''

    import FireClustering,FireVector
    from FireConsts import sleeperthresh  # CONNECTIVITY_THRESHOLD_KM

    # extract existing active fire data (use extending ranges)
    eafires     = [allfires.fires[fid]  for fid in fids_ea]
    # eafirerngs  = [FireVector.addbuffer(f.hull,CONNECTIVITY_THRESHOLD_KM[f.ftype]*1000) for f in eafires]
    eafirerngs = []
    for f in eafires:
        CONNECTIVITY_THRESHOLD_KM = get_CONNECTIVITY_THRESHOLD(f)
        eafirerngs.append(FireVector.addbuffer(f.hull,CONNECTIVITY_THRESHOLD_KM*1000))

    # extract existing active fire data (use hulls to avoid duplicate buffers)
    nefires     = [allfires.fires[fid]  for fid in fids_ne]
    nefirehulls = [f.hull for f in nefires]

    # inserting boxes into a spatial index
    ea_idx = FireClustering.build_rtree(eafirerngs)

    # loop over all fire objects that have newly expanded or formed, record merging fire id pairs
    fids_merge = []  # initialize the merged fire id pairs {source id:target id}
    firedone = {i:False for i in fids_ne}  # flag to mark an newly expanded fire obj that has been invalidated
    for id_ne in range(len(nefires)):
        fid_ne = fids_ne[id_ne]    # newly formed/expanded fire id
        if firedone[fid_ne] == False: # skip objects that have been merged to others in earlier loop
            # potential neighbors
            id_cfs = FireClustering.idx_intersection(ea_idx, nefirehulls[id_ne].bounds)
            # loop over all potential neighbor fobj candidiates
            for id_ea in id_cfs:
                fid_ea = fids_ea[id_ea]  # fire id of existing active fire
                # if fid_ne == fid_ea, skip;
                # if the expanded fire has been merged to a existing active fire, skip the rest loops
                if (fid_ne != fid_ea):
                    # if fire fmid is within distance of fire fid, two objects will merge
                    if nefirehulls[id_ne].intersects(eafirerngs[id_ea]):
                        # the fire id of neighboring active Fobj
                        # depending on which fid is smaller, merge the two fire objects in different directions
                        if fid_ea > fid_ne:  # merge fid_ea to fid_ne
                            fids_merge.append((fid_ea,fid_ne))
                            if fid_ea in firedone.keys():
                                firedone[fid_ea] = True  # remove fid_ea from the newly expanded fire list (since it has been invalidated)
                        else:            # merge fid_ne to fid_ea
                            fids_merge.append((fid_ne,fid_ea))
                            # fid_ne is merged to others, so stop it and check the next id_ne
                            ## technically the eafirerngs and nefirehulls have to be updated before the next loop happens
                            ## else it can happen that intersections are not being detected
                            ## if more than two fires grow together!!
                            #break

    # now check if any of the sleeper fires may have reactivated based on its fire line
    if len(fids_sleep) > 0: # check if there are potential sleepers
        # extract existing sleeping fires and their firelines
        sleepfires = [allfires.fires[fid]  for fid in fids_sleep]
        sleepflines = [f.fline for f in sleepfires]

        # inserting boxes into a spatial index
        nefirebuf  = [FireVector.addbuffer(hull,sleeperthresh*1000) for hull in nefirehulls]
        ne_idx = FireClustering.build_rtree(nefirebuf)

        # do the check analoguous to above
        firedone = {i:False for i in fids_sleep}  # flag to mark an sleeper fire obj that has been invalidated
        for id_sleep in range(len(sleepfires)):
            fid_sleep = fids_sleep[id_sleep]    # sleeper fire id
            if sleepflines[id_sleep] == None: # if there is no fire line (last active detection within), skip
                continue
            if firedone[fid_sleep] == False: # skip objects that have been merged to others in earlier loop
                # potential neighbors
                id_cfs = FireClustering.idx_intersection(ne_idx, sleepflines[id_sleep].bounds)
                # loop over all potential neighbour fobj candidates
                for id_ne in id_cfs:
                    fid_ne = fids_ne[id_ne]
                    if nefirebuf[id_ne].intersects(sleepflines[id_sleep]):
                        # depending on which fid is smaller, merge the two fire objects in different directions
                        if fid_ne > fid_sleep:  # merge new fire to sleeper, reactivate sleeper
                            fids_merge.append((fid_ne,fid_sleep))
                            if fid_ne in firedone.keys():
                                firedone[fid_ne] = True  # remove fid_ne from the newly expanded fire list (since it has been invalidated)
                        else:            # merge sleeper to new or expanded fire
                            fids_merge.append((fid_sleep,fid_ne))

    # loop over each pair in the fids_merge, and do modifications for both target and source objects
    #  - target: t_ed; pixels, newpixels, hull, extpixels
    #  - source: invalidated
    if len(fids_merge) > 0:
        # fids_merge needs to be corrected if several fires merge at once!
        # i.e. if fire 2 merges into fire 1 and fire 3 merges into fire 2
        # in this case not correcting fids_merge will lead to invalidation of fire 3!!!
        fids_merge =  correct_nested_ids(fids_merge)
        # logger.info(f'IDs to merge: {fids_merge}')

        for fid1,fid2 in fids_merge:
            #logger.info(f'Fire {fid1} was merged to Fire {fid2}')

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
            f_target.hull = FireVector.update_hull(phull,pextlocs+newlocs, sensor=sensor)

            # - use the updated hull to update exterior pixels
            f_target.extpixels = FireVector.cal_extpixels(f_target.extpixels+f_source.pixels,f_target.hull)

            # invalidate and deactivate source object
            f_source.invalid = True
            f_source.mergeid = f_target.mergeid
            #logger.info(f'Source id overwritten: {f_source.id} = {f_target.id}')

            # record the heritages
            allfires.heritages.append((fid1,fid2))

        # remove duplicates and record fid change for merged and invalidated
        fids_invalid,fids_merged = zip(*fids_merge)
        fids_merged = sorted(set(fids_merged))
        fids_invalid = sorted(set(fids_invalid))
        allfires.record_fids_change(fids_merged = fids_merged, fids_invalid = fids_invalid)

    return allfires

def Fire_Forward(tst,ted,restart=False,region=None,src='SNPP',nrt=False):
    ''' The wrapper function to progressively track all fire events for a time period
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
    '''

    # import libraries
    import FireObj
    import FireIO
    import os
    import glob
    import copy

    if src in ['SNPP', 'NOAA20', 'VIIRS', 'BAMOD']:
        sensor = 'viirs'
    elif src in ['BAMOD']:
        sensor = 'mcd64'
    else:
        print('Please set correct src value')

    # used to record time of script running
    import time
    t1 = time.time()
    t0 = t1

    # initialize allfires object (using previous day data or an empty fire object)
    # if FireObj.t_dif(tst,(tst[0],1,1,'AM'))==0:  # force restart at the start of a year
    #     restart = True
    allfires = Fobj_init(tst,region[0],restart=restart)

    # # clean temp folder (used for faster loading of active fire data)
    # temp_files = os.path.join(dirpjdata,'temp') + '/*'
    # files = glob.glob(temp_files)
    # for f in files:
    #     os.remove(f)

    # loop over all days during the period
    endloop = False  # flag to control the ending of the loop
    t = list(tst)    # t is the time (year,month,day,ampm) for each step
    while endloop == False:
        logger.info('')
        logger.info(t)

        # 1. record existing active fire ids (for the previous time step)
        fids_ea = allfires.fids_active

        # 2. update allfires and fire object changes due to temporal progression
        # all fires
        allfires.update_t(t)           # update t for allfires
        allfires.reset_fids_updated()  # reset fids_expanded, fids_new, fids_merged, fids_invalid
        # for each fire
        allfires.update_t_allfires(t)  # update t
        allfires.reset_newpixels()     # reset newpixels

        # 3. read active fire pixels from VIIRS dataset
        t_read = time.time()

        afp = FireIO.read_AFP(t,src=src,nrt=nrt,region=region)
        t_read2 = time.time()
        logger.info(f'reading file {(t_read2-t_read)}')
        if len(afp) > 0:
            t_expand = time.time()
            # 4. do fire expansion/creation using afp
            allfires = Fire_expand_rtree(allfires,afp,fids_ea)
            t_expand2 = time.time()
            logger.info(f'expanding fires {(t_expand2-t_expand)}')
            # 5. do fire merging using updated fids_ne and fid_ea
            fids_ne = allfires.fids_ne                         # new or expanded fires id
            fids_ea = sorted(set(fids_ea+allfires.fids_new))   # existing active fires (new fires included)
            fids_sleep = allfires.fids_sleeper
            t_merge = time.time()
            if len(fids_ne) > 0:
                allfires = Fire_merge_rtree(allfires,fids_ne,fids_ea,fids_sleep)
            t_merge2 = time.time()
            logger.info(f'merging fires {(t_merge2-t_merge)}')

        # 6. determine or update fire type
        allfires.updateftypes()

        # 7. manualy invalidate static fires (with exceptionally large fire density)
        allfires.invalidate_statfires()

        # 8. log and save
        #  - record fid_updated (the fid of fires that change in the time step) to allfires object and logger
        logger.info(f'fids_expand: {allfires.fids_expanded}')
        logger.info(f'fids_new: {allfires.fids_new}')
        logger.info(f'fids_merged: {allfires.fids_merged}')
        logger.info(f'fids_invalid: {allfires.fids_invalid}')


        # 9. loop control
        #  - if t reaches ted, set endloop to True to stop the next loop
        if FireObj.t_dif(t,ted)==0:
            endloop = True
            # correct fire heritage of last time step
            if len(allfires.heritages) > 0:
                allfires.heritages = correct_nested_ids(allfires.heritages)

        # 10. fire object save
        FireIO.save_fobj(allfires,t,region[0],activeonly=True)
        if FireObj.t_dif(t,ted)==0:
            FireIO.save_fobj(allfires,t,region[0],activeonly=False)

        # 11. update t with the next time stamp
        #  - record running times for the loop
        t2 = time.time()
        logger.info(f'{(t2-t1)/60.} minutes used to run alg {t}')
        t1 = t2

        # - move to next time step
        t = FireObj.t_nb(t,nb='next')

    # record total running time
    t3 = time.time()
    logger.info(f'This running takes {(t3-t0)/60.} minutes')

    return allfires

if __name__ == "__main__":
    ''' The main code to run time forwarding for a time period
    '''

    import time
    t1 = time.time()

    # set the start and end time
    tst=(2021,7,13,'AM')
    ted=(2021,9,15,'PM')

    region = ('Dixie',[-121.6,39.8,-120.5,40.6])
    # Run the time forward and record daily fire objects .pkl data and fire attributes .GeoJSON data
    Fire_Forward(tst=tst,ted=ted,restart=True,region=region,src='VIIRS',nrt=False)

    t2 = time.time()
    print(f'{(t2-t1)/60.} minutes used to run code')

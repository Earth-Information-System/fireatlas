# -*- coding: utf-8 -*-
"""
Benchmarking differnt versions of the functions that take longest

Created on Wed Nov  3 12:26:11 2021

@author: rebec
"""

#%% next benchmark: implementation of rtree

# instead of FireClustering.filter_centroid: FireClustering.filter_nearest
def build_rtree(geoms, fids=False):
    '''Builds Rtree from a shapely multipolygon shape
    and optionally uses list of fids as identifier'''
    import rtree
    
    idx = rtree.index.Index() # create new index
    for ind, geom in enumerate(geoms):
        if fids:
            idx.insert(ind, geom.bounds, fids[ind])
        else:
            idx.insert(ind, geom.bounds, ind)
    
    return idx

def idx_intersection(idx, geom):
    ''' 
    Finds all objects in an index whcih bounding boxes intersect with a geometry's bounding box
    '''
    bbox = geom.bounds
    intersections = list(idx.intersection(bbox, objects = True))
    fids, bbox = zip(*[(item.object, item.bbox) for item in intersections])
    return fids

def Fire_expand_rTree(allfires,afp,fids_ea):
    ''' Use daily new AF pixels to create new Fobj or combine with existing Fobj

    Parameters
    ----------
    allfires : Allfires obj
        the existing Allfires object for the time step
    afp : 3-element list
        (lat, lon, FRP) of new active fire pixels
    fids_ea : list
        fire ids of existing active fires at previous time step

    Returns
    -------
    allfires : Allfires obj
        updated Allfires object for the day with new formed/expanded fire objects
    '''

    import FireObj,FireClustering,FireVector
    from FireConsts import SPATIAL_THRESHOLD_KM,CONNECTIVITY_THRESHOLD_KM

    # record current time for later use (t in allfires has been updated in the Fire_Forward function)
    t = allfires.t

    # do some initializations
    idmax = allfires.number_of_fires-1  # maximum id of existing fires
    fids_expanded = []      # a list recording Fobj ids which has been expanded
    fids_new = []           # a list recording ids of new Fobjs

    # expanding ranges of existing active fires (extracted using fids_ea)
    eafires     = [allfires.fires[fid]  for fid in fids_ea]
    eafirerngs  = [FireVector.addbuffer(f.hull,CONNECTIVITY_THRESHOLD_KM[f.ftype]*1000) for f in eafires]
    
    # inserting boxes into a spatial index
    ea_idx = build_rtree(eafirerngs)

    # do preliminary clustering using new active fire locations (assign cid to each pixel)
    afp_loc = [(x,y) for x,y,z in afp]
    cid = FireClustering.do_clustering(afp_loc,SPATIAL_THRESHOLD_KM)  # this is the cluster id (starting from 0)
    logger.info(f'New fire clusters of {max(cid)} at this time step')

    # loop over each of the new clusters (0:cid-1) and determine its fate
    FP2expand = {}  # the diction used to record fire pixels assigned to existing active fires {fid:Firepixels}
    for ic in range(max(cid)+1):
        # create cluster object using all newly detected active fires within a cluster
        pixels = [afp[i] for i, v in enumerate(cid) if v==ic]
        cluster = FireObj.Cluster(ic,pixels,t)  # create a Cluster object using the pixel locations
        hull = cluster.hull  # the hull of the cluster

        # extract potential neighbors using centroid distance (used for prefilter)
        # id_cfs is the indices in eafirecents list, not the fire id
        # id_cfs = FireClustering.filter_centroid(cluster.centroid,eafirecents,MAX_THRESH_CENT_KM)
        id_cfs = idx_intersection(ea_idx, hull) # extract potential neighbours using spatial index

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
                    newFPs = [FireObj.FirePixel((p[0],p[1]),p[2],t,fmid) for p in pixels] # new FirePixels from the cluster
                    if fmid in FP2expand.keys():   # for a single existing object, there can be multiple new clusters to append
                        FP2expand[fmid] = FP2expand[fmid] + newFPs
                    else:
                        FP2expand[fmid] = newFPs

                    logger.info(f'Fire {fmid} expanded with pixels from new cluster {ic}')

                    fids_expanded.append(fmid) # record fmid to fid_expanded

                    clusterdone = True   # mark the cluster as done (no need to create new Fobj)

        # if this cluster can't be appended to any existing Fobj, create a new fire object using the new cluster
        if clusterdone is False:
            # create a new fire id and add it to the fid_new list
            id_newfire = idmax + 1
            logger.info(f'Fire {id_newfire} created with pixels from new cluster {ic}')
            fids_new.append(id_newfire)  # record id_newfire to fid_new

            # use the fire id and new fire pixels to create a new Fire object
            newfire = FireObj.Fire(id_newfire,t,pixels)

            # add the new fire object to the fires list in the Allfires object
            allfires.fires.append(newfire)

            # increase the maximum id
            idmax += 1

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

            # update the hull using previous hull and previous exterior pixels
            phull = f.hull
            pextlocs = [p.loc for p in f.extpixels]
            newlocs = [p.loc for p in newFPs]
            f.hull = FireVector.update_hull(phull,pextlocs+newlocs)  # use update_hull function to save time

            # use the updated hull to update exterior pixels
            f.extpixels = FireVector.cal_extpixels(f.extpixels+newFPs,f.hull)

    # remove duplicates and sort the fid_expanded
    fids_expanded = sorted(set(fids_expanded))

    # record fid change for expanded and new
    allfires.record_fids_change(fids_expanded=fids_expanded, fids_new=fids_new)

    # logger.info(f'In total, {len(fids_expanded)} fires expanded, and {len(fids_new)} fires created')

    return allfires

#%% VIIRS read function bottleneck
def read_VNP14ML04_clip2(t, ext = False):
    ''' Read monthly S-NPP VIIRS fire location (C1.04)

    Parameters
    ----------
    t : tuple, (int,int,int,str)
        the year, month, day and 'AM'|'PM' during the intialization

    Returns
    -------
    vlist : list
        (lat,lon) tuple of all daily active fires
        
        
    THIS IS SLOWER THAN clip
    '''
    from FireConsts import dirextdata

    from datetime import date
    import os
    import pandas as pd
    from shapely.geometry import Point
    import geopandas as gpd


    d = date(*t[:-1])
    # set monthly file name
    dirFC = os.path.join(dirextdata,'VNP14IMGML05') + '/'
    fnmFC = os.path.join(dirFC,'VNP14IMGML.'+d.strftime('%Y%m')+'.C1.05.txt')

    # read and extract
    ext = [65.5, -156.0, 66.3, -152.0]
    # shp_name = 'AK.shp'
    # shp = get_any_shp(shp_name) # shapefile of test region
    # ext = list(shp.bounds)
    # if shp_name == 'AK.shp': ext[2] = 140*-1  # Alaska extends over Antimeridian, therefore bounds are screwed up
    usecols = ['YYYYMMDD','HHMM','Lat','Lon','Sample','FRP','Confidence','Type','DNFlag']
    if os.path.exists(fnmFC):
        # create empty list
        df_list = []
        # read and filter dataframe
        with open(fnmFC) as fileobj:
            head = fileobj.readline().replace(" ", "").split(',')
            idx = [head.index(col) for col in usecols]
            for lines in fileobj:
                line = lines.replace(" ", "").split(',')
                cond1 = line[idx[6]] != 'low'
                cond2 = (line[idx[7]] == '0') | (line[idx[7]] == '3')
                cond3 = (float(line[idx[2]]) > ext[1]) & (float(line[idx[2]]) < ext[3])
                cond4 = (float(line[idx[3]]) > ext[0]) & (float(line[idx[3]]) < ext[2])
                if cond4: 
                    if cond3: 
                        if cond2: 
                            if cond1: 
                                df_list.append([line[i] for i in range(len(line)) if i in idx])
                
        newfirepixels = pd.DataFrame(df_list,columns=usecols)

        # in Collection 04, some FRP values are '*******', causing problems
        if newfirepixels.dtypes.FRP.name == 'object':
            newfirepixels.FRP = newfirepixels.FRP.replace('*******',0).astype('float')

        # temporal (daily) filter
        newfirepixels = newfirepixels.loc[(newfirepixels['YYYYMMDD']==d.strftime('%Y-%m-%d'))]

        # overpass time filter (AM or PM)
        # DNFlag - 0: night; 1: day
        # vlh = (df['HHMM']*0.01+df['lon']/15.)  # local time
        # tpm = (vlh > 7) & (vlh < 19)           # if local time in [7,19], set as PM overpass
        tpm = (newfirepixels['DNFlag'] == 1)
        if t[-1] == 'AM':
            newfirepixels = newfirepixels.loc[~tpm]
        elif t[-1] == 'PM':
            newfirepixels = newfirepixels.loc[tpm]

        # spatial filter
        point_data = [Point(xy) for xy in zip(newfirepixels['Lon'], newfirepixels['Lat'])]
        vmon_gpd = gpd.GeoDataFrame(newfirepixels, geometry=point_data)
        # vmon_gpd = vmon_gpd[vmon_gpd['geometry'].within(shp)]

        # convert to list of (lat,lon)
        vlist = vmon_gpd[['Lat','Lon','FRP']].values.tolist()

        return vlist
    else:
        return None

def read_VNP14ML04_clip2(t, ext = False):
    ''' Read monthly S-NPP VIIRS fire location (C1.04)

    Parameters
    ----------
    t : tuple, (int,int,int,str)
        the year, month, day and 'AM'|'PM' during the intialization

    Returns
    -------
    vlist : list
        (lat,lon) tuple of all daily active fires
        
        
    THIS IS SLOWER THAN ORIGINAL
    '''
    from FireConsts import dirextdata

    from datetime import date
    import os
    import pandas as pd
    from shapely.geometry import Point
    import geopandas as gpd


    d = date(*t[:-1])
    # set monthly file name
    dirFC = os.path.join(dirextdata,'VNP14IMGML05') + '/'
    fnmFC = os.path.join(dirFC,'VNP14IMGML.'+d.strftime('%Y%m')+'.C1.05.txt')

    # read and extract
    ext = [65.5, -156.0, 66.3, -152.0]
    # shp_name = 'AK.shp'
    # shp = get_any_shp(shp_name) # shapefile of test region
    # ext = list(shp.bounds)
    # if shp_name == 'AK.shp': ext[2] = 140*-1  # Alaska extends over Antimeridian, therefore bounds are screwed up
    usecols = ['YYYYMMDD','HHMM','Lat','Lon','Sample','FRP','Confidence','Type','DNFlag']
    if os.path.exists(fnmFC):
        # create empty list
        df_list = []
        # read and filter dataframe
        with open(fnmFC) as fileobj:
            head = fileobj.readline().replace(" ", "").split(',')
            idx = [head.index(col) for col in usecols]
            for lines in fileobj:
                line = lines.replace(" ", "").split(',')
                cond1 = line[idx[6]] != 'low'
                cond2 = (line[idx[7]] == '0') | (line[idx[7]] == '3')
                cond3 = (float(line[idx[2]]) > ext[1]) & (float(line[idx[2]]) < ext[3])
                cond4 = (float(line[idx[3]]) > ext[0]) & (float(line[idx[3]]) < ext[2])
                if cond4: 
                    if cond3: 
                        if cond2: 
                            if cond1: 
                                df_list.append([line[i] for i in range(len(line)) if i in idx])
                
        newfirepixels = pd.DataFrame(df_list,columns=usecols)

        # in Collection 04, some FRP values are '*******', causing problems
        if newfirepixels.dtypes.FRP.name == 'object':
            newfirepixels.FRP = newfirepixels.FRP.replace('*******',0).astype('float')

        # temporal (daily) filter
        newfirepixels = newfirepixels.loc[(newfirepixels['YYYYMMDD']==d.strftime('%Y-%m-%d'))]

        # overpass time filter (AM or PM)
        # DNFlag - 0: night; 1: day
        # vlh = (df['HHMM']*0.01+df['lon']/15.)  # local time
        # tpm = (vlh > 7) & (vlh < 19)           # if local time in [7,19], set as PM overpass
        tpm = (newfirepixels['DNFlag'] == 1)
        if t[-1] == 'AM':
            newfirepixels = newfirepixels.loc[~tpm]
        elif t[-1] == 'PM':
            newfirepixels = newfirepixels.loc[tpm]

        # spatial filter
        point_data = [Point(xy) for xy in zip(newfirepixels['Lon'], newfirepixels['Lat'])]
        vmon_gpd = gpd.GeoDataFrame(newfirepixels, geometry=point_data)
        # vmon_gpd = vmon_gpd[vmon_gpd['geometry'].within(shp)]

        # convert to list of (lat,lon)
        vlist = vmon_gpd[['Lat','Lon','FRP']].values.tolist()

        return vlist
    else:
        return None
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 28 10:52:41 2022

@author: rebec
"""

def find_all_end(tst, ted):
    ''' find all final perimeters points in the half-daily snapshots and
    save them to one gdf

    Parameters
    ----------
    tst : Start date of fire season
    ted: End date of fire season

    Returns
    -------
    gdf_all : geopandas DataFrame
        the gdf containing half daily fire basic attributes and final perimeter
    '''
    import FireObj,FireIO
    
    # initialise list of already checked fire ids
    checked_ids = []
    id_ted_dict = {}
    
    endloop = False  # flag to control the ending of the loop
    creategdf = True # needed since ted cannot be used anymore
    t = list(ted)    # t is the time (year,month,day,ampm) for each step
    while endloop == False:
        #print(t)
        # read daily gdf
        if FireIO.check_gdfobj(t,op=''):
            gdf = FireIO.load_gdfobj(t)
            gdf_active = gdf[gdf.isactive == 1]
    
            # append daily row to gdf_all
            if creategdf:
                gdf_all = gdf_active
                creategdf = False
            else:
                # exclude ids that have an older active perimeter
                gdf_active = gdf_active[~gdf_active['mergid'].isin(checked_ids)]
                gdf_all = gdf_all.append(gdf_active)
            
            # add ids that have been written out to checked_ids, these will be skipped next time
            if len(gdf_active) > 0:
                for id in gdf_active['mergid'].tolist():
                    id_ted_dict[id] = tuple(t)
                checked_ids = list(set(id_ted_dict))
        
        #  - if t reaches tst, set endloop to True to stop the loop
        if FireObj.t_dif(t,tst)==0:
            endloop = True

        #  - update t with the previous time stamp
        t = FireObj.t_nb(t,nb='previous')

    return id_ted_dict, gdf_all

def create_final_allfires(id_ted_dict):
    '''Loop through the end-dates of all fires and extract the fire objects
    at their last active fire detection from the according allfires object
    
    Parameters
    ----------
    id_ted_dict : dictionary
        dictionary containing the end date of each fire (value), accessed by fire id (key)
    
    Returns
    -------
    allfires_final : Allfires obj
        Allfires obj containing all final fires at their last active time step
    '''
    
    import FireIO, FireObj
    
    # initialise new allfire object
    ted = max(id_ted_dict.values())
    allfires_final = FireObj.Allfires(ted)
    
    # find all dates for loop
    dates = list(set(id_ted_dict.values()))
    
    for t in dates:
        allfires = FireIO.load_fobj(t)
        
        # find all fires that end on that date
        fireids = [key for key in id_ted_dict if id_ted_dict[key] == t]
        
        # loop through fireids and append to new allfires object
        for fireid in fireids:
            allfires_final.fires.append(allfires.fires[fireid])
    
    return allfires_final

def make_gdf_fperim_simple(allfires):
    ''' Create geopandas DataFrame for fire basic attributes and fire perimeter (hull).
            Use saved gdf files for previous time step and update active fires only.

    Parameters
    ----------
    allfires : Allfires object
        the Allfires object to be used to modify gdf

    Returns
    -------
    gdf : geopandas DataFrame
        the gdf containing half daily fire basic attributes and fire perimeter
    '''
    import geopandas as gpd

    # initialize the GeoDataFrame with fire attribute names (in dd) as columns
    gdf = gpd.GeoDataFrame(columns=['fid','mergid'],crs='epsg:4326', geometry=[])
    heritage = dict(allfires.heritages)

    # update the hull of each active fire as the geometry column
    for fire_no in range(allfires.number_of_activefires):
        fid = allfires.fires[fire_no].id
        print(fid)
        gdf.loc[fid,'fid'] = fid
        if fid in heritage.keys():
            gdf.loc[fid,'mergid'] = heritage[fid]
        else:
            gdf.loc[fid,'mergid'] = fid
        fhull = allfires.fires[fire_no].hull
        if fhull is not None:
            if fhull.geom_type == 'MultiPolygon':
                gdf.loc[fid,'geometry'] = gpd.GeoDataFrame(geometry=[fhull]).geometry.values
            else:
                gdf.loc[fid,'geometry'] = fhull
    gdf = gdf.dissolve(by = 'mergid')
    # make sure the data types are correct
    gdf['fid'] = gdf.index

    return gdf

def Fire_merge_final(allfires):
    ''' For final fire objects that overlap in time and space, merge them
    
    Parameters
    ----------
    allfires : Allfires obj
        the existing Allfires object for the time step
    
    Returns
    -------
    allfires : Allfires obj
        Allfires obj after fire merging
    '''

    import FireClustering,FireVector,FireMain
    from FireConsts import CONNECTIVITY_THRESHOLD_KM
    
    # extract existing active fire data (use extending ranges)
    fids_active = range(len(allfires.fids_active))
    fires     = [allfires.fires[fid] for fid in fids_active]
    firehulls = [f.hull for f in fires]
    firerngs  = [FireVector.addbuffer(f.hull,CONNECTIVITY_THRESHOLD_KM[f.ftype]*1000) for f in fires]
    
    # inserting boxes into a spatial index
    idx = FireClustering.build_rtree(firerngs)

    # loop over all fire objects that have newly expanded or formed, record merging fire id pairs
    fids_merge = []  # initialize the merged fire id pairs {source id:target id}
    firedone = {i:False for i in fids_active}  # flag to mark an newly expanded fire obj that has been invalidated
    for id_a in range(len(fids_active)):
        fid_a = fids_active[id_a]
        if firedone[fid_a] == False: # skip objects that have been merged to others in earlier loop
            # potential neighbors
            id_cfs = FireClustering.idx_intersection(idx, firerngs[id_a].bounds)
            # loop over all potential neighbor fobj candidates
            for id_b in id_cfs:
                fid_b = fids_active[id_b]
                if (fid_a != fid_b):
                    cond1 = fires[fid_a].t_st <= fires[fid_b].t_ed
                    cond2 = fires[fid_a].t_ed >= fires[fid_b].t_st
                    if cond1 and cond2:
                        if firehulls[id_b].intersects(firerngs[id_a]):
                            true_a = fires[fid_a].id
                            true_b = fires[fid_b].id
                            if true_b > true_a:  # merge fid_ea to fid_ne
                                fids_merge.append((fid_b,fid_a))
                                if fid_b in firedone.keys():
                                    firedone[fid_b] = True  # remove fid_ea from the newly expanded fire list (since it has been invalidated)
                            else:
                                fids_merge.append((fid_a,fid_b))
       
    # loop over each pair in the fids_merge, and do modifications for both target and source objects
    #  - target: t_ed; pixels, newpixels, hull, extpixels
    #  - source: invalidated 
    if len(fids_merge) > 0:
        fids_merge =  FireMain.correct_nested_ids(fids_merge)
        
        for fid1,fid2 in fids_merge:
            # update source and target objects
            f_source = allfires.fires[fid1]
            f_target = allfires.fires[fid2]

            # - target fire t_ed set to current time
            f_target.t_ed = allfires.t

            # - target fire add source pixels to pixels and newpixels
            f_target.pixels = f_target.pixels + f_source.pixels
            f_target.newpixels = f_target.newpixels + f_source.newpixels

            # - update the hull using previous hull and previous exterior pixels
            phull = f_target.hull
            pextlocs = [p.loc for p in f_target.extpixels]
            newlocs = [p.loc for p in f_source.pixels]
            f_target.hull = FireVector.update_hull(phull,pextlocs+newlocs, sensor='viirs')

            # - use the updated hull to update exterior pixels
            f_target.extpixels = FireVector.cal_extpixels(f_target.extpixels+f_source.pixels,f_target.hull)

            # invalidate and deactivate source object
            f_source.mergeid = f_target.mergeid

            # record the heritages
            allfires.heritages.append((fires[fid1].id,fires[fid2].id))

        # remove duplicates and record fid change for merged and invalidated
        fids_invalid,fids_merged = zip(*fids_merge)
        fids_merged = sorted(set(fids_merged))
        fids_invalid = sorted(set(fids_invalid))
        allfires.record_fids_change(fids_merged = fids_merged, fids_invalid = fids_invalid)

    return allfires

def add_mcd64(id_ted_dict,ext):
    '''Retrieve final fire objects and add MCD64 pixels
    this works analoguous to the fire tracking algorithm:
        mcd64 pixels are clustered and used to extend existing viirs fire objects
        burned area pixels that are not close to a viirs fire are ignored
        fires that are close or overlapping after adding burned area data are merged
    
    Parameters
    ----------
    id_ted_dict : dict
        dictionary containing the end date of each fire (value), accessed by fire id (key)
    ext: list
        extent (geometry.bounds) of the output gdf
    
    Returns
    -------
    gdf : geopandas geodataframe
        geodataframe containing the updated fire perimeters based on viirs AF and mcd64
    
    '''
    
    import FireIO, FireMain
    import datetime
    
    # load fire objects based on VIIRS only
    allfires = create_final_allfires(id_ted_dict)
    
    year = allfires.t[0]
    
    # 1. read active fire and burned area pixels
    afp = FireIO.read_mcd64_pixels(year,ext)    
    
    # now we need to loop through fires to make sure only pixels from that timestep are added
    for fid in range(allfires.number_of_activefires):
        fire = allfires.fires[fid]
        ted_fire = datetime.datetime(*id_ted_dict[fire.id][:3]).timetuple().tm_yday
        tst_fire = datetime.datetime(*fire.t_st[:3]).timetuple().tm_yday
        # print(fid, fire.id)
        
        # filter active fire points
        afp_fire = afp.loc[(afp['doy'] > (tst_fire - 5)) & (afp['doy'] < (ted_fire + 5))]
        if len(afp_fire) > 0:
            afp_fire = [(row['lat'],row['lon'],0,0,0) for idx,row in afp_fire.iterrows()]
            
            # 4. expand that fire using afp
            allfires_new = FireMain.Fire_expand_rtree(allfires,afp_fire,[fid],sensor='mcd64',expand_only=True,log=False)
            #print(fire.id, allfires_new.fires[fid].id)
            
            # replace the old allfire with the new one
            allfires.fires[fid] = allfires_new.fires[fid]
            
    # next we merge fires with spatial and temporal overlap
    allfires = Fire_merge_final(allfires)     
    
    # 7. make fire perimeters
    gdf = make_gdf_fperim_simple(allfires)
    gdf = gdf.set_index('fid')
    gdf['mergid'] = gdf.index
    
    # 8. save
    FireIO.save_gdfobj(gdf,allfires.t,param='final_mcd64')
    
    return gdf
    
def save_gdf_trng(tst,ted):
    ''' Wrapper to create and save all ignitions as gdf
        1) final perimeters from VIIRS only
        2) final perimeters combining VIIRS and MCD64
        3) final perimeters with clipped lakes

    Parameters
    ----------
    tst : tuple, (int,int,int,str)
        the year, month, day and 'AM'|'PM' at start time
    ted : tuple, (int,int,int,str)
        the year, month, day and 'AM'|'PM' at end time
    '''
    import FireIO, PostProcess
    
    # find all final perimeters and write out to gdf
    id_ted_dict,gdf = find_all_end(tst, ted)
    FireIO.save_gdfobj(gdf,tst,param='final')
    
    # add mcd64 data to the final perimeters
    ext = gdf.geometry.total_bounds
    gdf_mcd = add_mcd64(id_ted_dict, ext)
    
    # clip lakes
    gdf_clip, gdf_lakes = PostProcess.clip_lakes_final(gdf_mcd, tst)
    FireIO.save_gdfobj(gdf_clip,tst,param='final_lakes')
    FireIO.save_gdfobj(gdf_lakes,tst,param='lakes')


if __name__ == "__main__":
    ''' The main code to record time series of geojson data for a fire
    '''

    import time
    t1 = time.time()
    # set the start and end time
    tst=(2020,6,1,'AM')
    ted=(2020,8,31,'PM')

    # for each day during the period,
    # create and save geojson files for temporal evolution of large fires
    save_gdf_trng(tst=tst,ted=ted)

    t2 = time.time()
    print(f'{(t2-t1)/60.} minutes used to run code')
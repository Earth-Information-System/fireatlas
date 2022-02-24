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
    # make sure the data types are correct
    gdf['fid'] = gdf['fid'].astype(int)
    gdf['mergid'] = gdf['mergid'].astype(int)

    return gdf

def add_mcd64(id_ted_dict,ext):
    
    import FireIO, FireMain
    import datetime
    # load viirs fireobjects
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
            
            # 4. rxpand that fire using afp
            allfires_new = FireMain.Fire_expand_rtree(allfires,afp_fire,[fid],sensor='mcd64',expand_only=True,log=False)
            #print(fire.id, allfires_new.fires[fid].id)
            
            # replace the old allfire with the new one
            allfires.fires[fid] = allfires_new.fires[fid]
        
    # # 5. do fire merging using updated fids_ne and fid_ea
    # fids_ne = allfires.fids_ne                         # new or expanded fires id
    # fids_ea = sorted(set(fids_ea+allfires.fids_new))   # existing active fires (new fires included)
    # allfires = FireMain.Fire_merge_rtree(allfires,fids_ne,fids_ea,sensor='mcd64')
    
    # # 6. correct heritages in case additional merges have happened
    # allfires.heritages = FireMain.correct_nested_ids(allfires.heritages)
    # # heritages have to me translated back to actual ids!!!
    # herit_corr = []
    # for tup in allfires.heritages:
    #     herit_corr.append((allfires.fires[tup[0]].id, allfires.fires[tup[1]].id))
    # allfires.heritages = herit_corr
    
    # 7. make fire perimeters
    gdf = make_gdf_fperim_simple(allfires)
    
    # 8. save
    FireIO.save_gdfobj(gdf,allfires.t,param='final_mcd64')
    
    return gdf
    
def save_gdf_trng(tst,ted):
    ''' Wrapper to create and save all ignitions as gdf

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
    
    # find water tiles
    #tiles = PostProcess.all_water_tiles(gdf)
    gdf_clip, gdf_lakes = PostProcess.clip_lakes_final(gdf, tst)

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
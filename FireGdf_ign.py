""" FireGdf_sfs
Module for creating geojson summary for temporal evolution of each single fire
        (large active fires  up to the present time step)
Large fires are fires with area > falim (defined in FireConsts)

List of geojson types
---------------------
* fperim(''): fire basic attributes and perimeter geometry
* fline('FL'): active fire line geometry
* NFP('NFP'): new fire pixels

List of functions
-----------------
* make_fire_history
* update_fire_history
* save_gdf_1t
* save_gdf_trng
* yrend_clean

Modules required
----------------
* FireObj
* FireIO
* FireConsts
"""

def clip_lakes(gdf,t,region):
    import shapely
    import FireIO,FireClustering,PostProcess
    
    # transform to lake crs
    gdf = gdf.to_crs('EPSG:3571')
    createds = True
    
    # loop through ignition points and clip
    for mergid in set(gdf.mergid):
        print('merge id:',mergid)
        lake_geoms = FireIO.load_lake_geoms(t,mergid,region=region) # load lake file
        if not isinstance(lake_geoms, type(None)): 
            # print(mergid)
            # lake_geoms is None if no lakes are within final perimeter
            if not isinstance(lake_geoms, shapely.geometry.multipolygon.MultiPolygon):
                # if only one lake is present we put it in a list for the loop
                lake_geoms = [lake_geoms]
            lake_idx = FireClustering.build_rtree(lake_geoms) # build index
            # loop through ips of that fire and clip lakes
            gdf_1f = gdf.loc[gdf.mergid == mergid] # single out the ignition points of that fire
            # print(len(gdf_1f))
            for fid in gdf_1f.index:
                gdf_1ip = gdf.loc[[fid]]
                gdf_1ip = PostProcess.clip_lakes_1fire_outer(gdf_1ip,lake_geoms,lake_idx)
                # append daily row to gdf_all
                if createds:
                    gdf_clip = gdf_1ip
                    createds = False
                else:
                    gdf_clip = gdf_clip.append(gdf_1ip)
        else:
            if createds:
                gdf_clip = gdf.loc[gdf.mergid == mergid]
                createds = False
            else:
                gdf_clip = gdf_clip.append(gdf.loc[gdf.mergid == mergid])
        print(len(gdf_clip))
    
    return gdf_clip

def find_all_ign(tst, ted,region=''):
    ''' find all ignition points in the half-daily snapshots and
    save them to one gdf

    Parameters
    ----------
    tst : Start date of fire season
    ted: End date of fire season

    Returns
    -------
    gdf_all : geopandas DataFrame
        the gdf containing half daily fire basic attributes and ignition perimeter
    '''
    import FireObj,FireIO

    endloop = False  # flag to control the ending of the loop
    creategdf = True
    t = list(tst)    # t is the time (year,month,day,ampm) for each step
    while endloop == False:
        #print(t)
        # read daily gdf
        if FireIO.check_gdfobj(t,op='',region=region):
            gdf = FireIO.load_gdfobj(t,region=region)
            gdf_ign = gdf[(gdf.isignition == 1)] # filter for ignitions
            
            # append daily row to gdf_all
            if creategdf:
                gdf_all = gdf_ign
                creategdf = False
            else:
                gdf_all = gdf_all.append(gdf_ign)

        #  - if t reaches ted, set endloop to True to stop the loop
        if FireObj.t_dif(t,ted)==0:
            endloop = True

        #  - update t with the next time stamp
        t = FireObj.t_nb(t,nb='next')
    
    # clip lakes
    gdf_all = clip_lakes(gdf_all,tst,region)
    
    return gdf_all

def save_gdf_trng(tst,ted,region=''):
    ''' Wrapper to create and save all ignitions as gdf

    Parameters
    ----------
    tst : tuple, (int,int,int,str)
        the year, month, day and 'AM'|'PM' at start time
    ted : tuple, (int,int,int,str)
        the year, month, day and 'AM'|'PM' at end time
    '''
    import FireIO
    gdf_ign = find_all_ign(tst,ted,region=region)
    FireIO.save_gdfobj(gdf_ign,tst,param='ign',region=region)
        


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

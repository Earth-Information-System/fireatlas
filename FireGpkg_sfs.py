""" FireGdf_sfs
Module for creating geojson summary for temporal evolution of each single fire
        (large active fires  up to the present time step)
Large fires are fires with area > falim

List of geojson types
---------------------
* fperim(''): fire basic attributes and perimeter geometry
* fline('FL'): active fire line geometry
* NFP('NFP'): new fire pixels

List of functions
-----------------
* make_fire_history
* merge_fires
* save_gdf_1fire
* save_gdf_trng

Modules required
----------------
* FireObj
* FireIO
* FireConsts
"""

def merge_fires(gdf_fid, fid, tst, t):
    ''' merge fire perimeters of all fires with the same merge id

    Parameters
    ----------
    gdf : geopandas DataFrame
        the gdf containing all entries that should be merged
    fid: the mergeid of the fire

    Returns
    -------
    gdf_diss : geopandas DataFrame
        the gdf containing a merged geometry and summary stats
    '''
    from datetime import date
    import FireObj

    dd = {'fireID':'int',                  # id
          'mergid':'int',               # this is the id in the large fire database
          #'clat':'float',               # centroid latitude   -> centroid[0]
          #'clon':'float',               # centroid longitude  -> centroid[1]
          'ftype':'int',                # fire type
          'isactive':'int',             # active status
          't_inactive':'int',           # how long has it been inactive
          'isignition':'int',           # is this a new ignition?
          'invalid':'int',              # invalid status
          'n_pixels':'int',             # number of total pixels
          'n_newpixels':'int',          # number of new pixels
          'farea':'float',              # fire size
          'fperim':'float',             # fire perimeter length
          'flinelen':'float',           # active fire front line length
          'duration':'float',           # fire duration
          'pixden':'float',             # fire pixel density
          'meanFRP':'float',            # mean FRP of the new fire pixels
          'tst_year':'int',             # t_st[0]
          'tst_month':'int',
          'tst_day':'int',
          'tst_ampm':'str',
          'ted_year':'int',             # t_ed[0]
          'ted_month':'int',
          'ted_day':'int',
          'ted_ampm':'str',
          'ted_doy':'int',
          }

    # dissolve the dataframe
    gdf_diss = gdf_fid.dissolve(by = 'mergid')

    # replace some values with sums
    if 'mergid' in dd.keys():
        gdf_diss['mergid'] = fid
    if 'n_pixels' in dd.keys():
        gdf_diss.loc[fid,'n_pixels'] = sum(gdf_fid.n_pixels)
    if 'n_newpixels' in dd.keys():
        gdf_diss.loc[fid,'n_newpixels'] = sum(gdf_fid.n_newpixels)
    if 'farea' in dd.keys():
        gdf_diss.loc[fid,'farea'] = sum(gdf_fid.farea)
    if 'fperim' in dd.keys():
        gdf_diss.loc[fid,'fperim'] = sum(gdf_fid.fperim)
    if 'flinelen' in dd.keys():
        gdf_diss.loc[fid,'flinelen'] = sum(gdf_fid.flinelen)
    if 'tst_year' in dd.keys():
        gdf_diss.loc[fid,'tst_year'] = int(tst[0])
    if 'tst_month' in dd.keys():
        gdf_diss.loc[fid,'tst_month'] = int(tst[1])
    if 'tst_day' in dd.keys():
        gdf_diss.loc[fid,'tst_day'] = int(tst[2])
    if 'tst_ampm' in dd.keys():
        gdf_diss.loc[fid,'tst_ampm'] = tst[3]
    if 'ted_year' in dd.keys():
        gdf_diss.loc[fid,'ted_year'] = int(t[0])
    if 'ted_month' in dd.keys():
        gdf_diss.loc[fid,'ted_month'] = int(t[0])
    if 'ted_day' in dd.keys():
        gdf_diss.loc[fid,'ted_day'] = int(t[0])
    if 'ted_ampm' in dd.keys():
        gdf_diss.loc[fid,'ted_ampm'] = t[3]
    if 'ted_doy' in dd.keys():
        gdf_diss.loc[fid,'ted_doy'] = date(*t[:-1]).timetuple().tm_yday
    if 'duration' in dd.keys():
        gdf_diss.loc[fid,'duration'] = FireObj.t_dif(tst,t) + 0.5

    # weighted average computed for averages
    if ('pixden' in dd.keys()) & ('farea' in dd.keys()):
        gdf_fid = gdf_fid.assign(pixweigh = gdf_fid['pixden']*gdf_fid['farea']) ## IS THIS CORRECT?
        gdf_diss.loc[fid,'pixden'] = sum(gdf_fid.pixweigh)/sum(gdf_fid.farea)
    if ('meanFRP' in dd.keys()) & ('n_newpixels' in dd.keys()):
        gdf_fid = gdf_fid.assign(FRPweigh = gdf_fid['meanFRP']*gdf_fid['n_newpixels'])
        gdf_diss.loc[fid,'meanFRP'] = sum(gdf_fid.FRPweigh)/sum(gdf_fid.n_newpixels)

    return(gdf_diss)

def make_fire_history(fid, start_end, fh, regnm):
    ''' derive time series of single fire attributes using the half-daily gdf summary

    Parameters
    ----------
    fire : Fire object
        the fire need to update
    start_end: tuple
        start and end date and AM/PM of the fire object

    Returns
    -------
    gdf_all : geopandas DataFrame
        the gdf containing half daily fire basic attributes and fire perimeter
    '''
    import math
    import FireObj,FireIO, PostProcess, FireClustering
    import shapely

    #fid = fire.id
    tst,ted = start_end  # start and ending time of the fire
    # lake_process = True
    #
    # # load lakes file and enter into index
    # lake_geoms = FireIO.load_lake_geoms(tst, fid)
    # if isinstance(lake_geoms, type(None)): # lake_geoms is None if no lakes are within final perimeter
    #     lake_process = False
    # else:
    #     # if only one lake is present we put it in a list for the loop
    #     if not isinstance(lake_geoms, shapely.geometry.multipolygon.MultiPolygon):
    #         lake_geoms = [lake_geoms]
    #         # here we could tell it to only build an index when number of features > X
    #     lake_idx = FireClustering.build_rtree(lake_geoms)

    endloop = False  # flag to control the ending of the loop
    t_inact = 0      # counts the days a fire is inactive before it spreads again
    t = list(tst)    # t is the time (year,month,day,ampm) for each step
    while endloop == False:
        #print(t)
        # read daily gdf
        gdf_perim = FireIO.load_gpkgobj(t, regnm,layer='perimeter')

        if gdf_perim is None:
            t = FireObj.t_nb(t,nb='next')
            t_inact += 0.5 # half-daily time steps
            continue

        gdf_1d = gdf_perim[gdf_perim.mergid.isin(fh)] # here we need to enter merge id instead of index

        # skip date if no new pixels
        gdf_1d = gdf_1d[gdf_1d.n_newpixels > 0]
        if len(gdf_1d) == 0:
            t = FireObj.t_nb(t,nb='next')
            t_inact += 0.5 # half-daily time steps
            continue

        # merge if several ids
        if len(gdf_1d) > 1:
            gdf_1d = merge_fires(gdf_1d, fid, tst, t)

        # # clip lakes and update attributes
        # if lake_process:
        #     gdf_1d = PostProcess.clip_lakes_1fire_outer(gdf_1d, lake_geoms, lake_idx)

        # change index to date
        gdf_1d.index = [FireObj.t2dt(t)]

        # update t_inactive
        gdf_1d.t_inactive = math.floor(t_inact) # inactivity is counted in full days

        # append daily row to gdf_all
        if FireObj.t_dif(t,tst)==0:
            gdf_all = gdf_1d
        else:
            gdf_all = gdf_all.append(gdf_1d)

        #  - if t reaches ted, set endloop to True to stop the loop
        if FireObj.t_dif(t,ted)==0:
            endloop = True

        #  - update t with the next time stamp
        t_inact = 0 # reset inactivity counter
        t = FireObj.t_nb(t,nb='next')


    # use correct time column as index
    # gdf_all.index = FireObj.ftrange(tst,ted)

    # remove fid column to save space
    # gdf_all = gdf_all.drop(columns='fid')

    return gdf_all

def make_fire_history_fline_NFP(fid, start_end,fh, regnm,layer='fireline'):
    ''' derive time series of single fire fireline or NFP using the half-daily gdf summary

    Parameters
    ----------
    fire : Fire object
        the fire need to update
    start_end: tuple
        start and end date and AM/PM of the fire object

    Returns
    -------
    gdf_all : geopandas DataFrame
        the gdf containing half daily fire basic attributes and fire perimeter
    '''
    import math
    import FireObj,FireIO, PostProcess, FireClustering
    import shapely

    tst,ted = start_end  # start and ending time of the fire

    endloop = False  # flag to control the ending of the loop
    t = list(tst)    # t is the time (year,month,day,ampm) for each step
    while endloop == False:
        #print(t)
        # read daily gdf
        gdf_perim = FireIO.load_gpkgobj(t, regnm,layer=layer)
        if gdf_perim is not None:
            # gdf_1d = gdf_perim[gdf_perim.mergid == fid] # here we need to enter merge id instead of index
            gdf_1d = gdf_perim[gdf_perim.mergid.isin(fh)]

            # skip date if no geometry is None
            gdf_1d = gdf_1d[gdf_1d.geometry != None]
            if len(gdf_1d) > 0:


                # # merge if several ids
                # if len(gdf_1d) > 1:
                #     gdf_1d = merge_fires(gdf_1d, fid, tst, t)

                # change index to date
                gdf_1d.index = [FireObj.t2dt(t)]


                # append daily row to gdf_all
                if FireObj.t_dif(t,tst)==0:
                    gdf_all = gdf_1d
                else:
                    gdf_all = gdf_all.append(gdf_1d)

        #  - if t reaches ted, set endloop to True to stop the loop
        if FireObj.t_dif(t,ted)==0:
            endloop = True

        #  - update t with the next time stamp
        t = FireObj.t_nb(t,nb='next')

    return gdf_all

def save_gdf_1fire(fid, start_end, fh, regnm):
    ''' derive and save time series of fire perimeter and fine line
        for each large active fire at a time step

    Parameters
    ----------
    fid : int
        the large fire id to be written out
    start_end: dict
        dictionary containing start and end dates for all large fires
    '''
    import FireIO

    gdf_sf = make_fire_history(fid, start_end, fh, regnm)
    gdf_sf_fline = make_fire_history_fline_NFP(fid, start_end,fh, regnm,layer='fireline')
    gdf_sf_NFP = make_fire_history_fline_NFP(fid, start_end,fh, regnm,layer='newfirepix')
    FireIO.save_gpkgsfs(gdf_sf,start_end[1],fid,regnm,
                        gdf_fline=gdf_sf_fline,gdf_NFP=gdf_sf_NFP)

def save_gdf_trng(ted,regnm,falim=4):
    ''' Wrapper to create and save gpkg files for all single large fires up to ted

    Parameters
    ----------
    ted : tuple, (int,int,int,str)
        the year, month, day and 'AM'|'PM' at end time
    falim : float
        the area threshold (in km2) for determining large fires
    '''
    import FireIO

    # read allfires
    allfires = FireIO.load_fobj(ted,regnm)

    # extract large valid fires (including active and inactive)
    large_ids = [f.fireID for f in allfires.validfires if f.farea > falim]

    # remove fires merged to other fires (is this really needed? all valid fires have not been merged to other fires... )
    skip,merge_ids = zip(*FireIO.load_fobj(ted,regnm).heritages)
    large_ids = [fire for fire in large_ids if fire not in skip]

    # start/end dates for all large fires (pre-calculate this list to save time-no need to do this within)
    # start_ends = {i:(allfires.fires[i].t_st, allfires.fires[i].t_ed) for i in large_ids}

    # loop over large fires
    for fid in large_ids:
        # start_end = start_ends[fid]
        f = allfires.fires[fid]

        fh = [fid]
        for h in allfires.heritages:
            if h[1] == fid:
                fh.append(h[0])

        start_end = (f.t_st,f.t_ed)
        save_gdf_1fire(fid, start_end, fh, regnm)


if __name__ == "__main__":
    ''' The main code to record time series of geojson data for a fire
    '''
    # import sys
    # sys.path.insert(1, 'D:/fire_atlas/1_code/fireatlas')

    import time
    t1 = time.time()

    # set the start and end time
    tst=(2021,7,13,'AM')
    ted=(2021,9,15,'PM')

    # for each day during the period,
    # create and save geojson files for temporal evolution of large fires
    save_gdf_trng(ted=ted,regnm='Dixie')

    t2 = time.time()
    print(f'{(t2-t1)/60.} minutes used to run code')

    # # only run at year end to clean the files
    # for year in range(2012,2020):
    #     print(year)
    #     yrend_clean(year)

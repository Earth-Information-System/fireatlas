""" FireGpkg
Module for creating regional geopackage summary at each time step

List of geojson types
---------------------
* fperim(''): fire basic attributes and perimeter geometry
* fline('FL'): active fire line geometry
* NFP('NFP'): new fire pixels

List of functions
-----------------
* make_gdf_fperim
* make_gdf_fline
* make_gdf_NFP
* save_gdf_1t
* save_gdf_trng

Modules required
----------------
* FireObj
* FireIO
* FireConsts
"""

def make_gdf_fperim(allfires, heritage, regnm):
    ''' Create gpd DataFrame for fire basic attributes and fire perimeter.

    Method:
    -------
    Extract attributes rom the allfires obj input, and create a gdf DataFrame.
    In order to save the running time, if the gpkg file for previous time step
    is available, read the data and only modify properties of previously active
    fires.

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
    # from FireConsts import dd
    import FireIO, FireObj

    # diagnostic data name and types saved in geojson files (in addition to geometries)
    # note this list can be expanded
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

    # read or initialize the gdf
    t_pt = FireObj.t_nb(allfires.t,nb='previous')
    gdf = None
    if allfires.t[0] == t_pt[0]: # don't read at year start
        gdf = FireIO.load_gpkgobj(t_pt,regnm,layer='perimeter')
    if gdf is None: # when no previous time step data, initialize the GeoDataFrame
        gdf = gpd.GeoDataFrame(columns=list(dd.keys()),crs='epsg:4326', geometry=[])
        gdf = gdf.set_index('fireID')   # set fid column as the index column
    else:
        gdf['isactive'] = 0  # pre-set all fires to inactive and update true status below

    # update data for each active fire using the allfires object
    for fid in (allfires.fids_active):
        #gdf.loc[fid,'clat'] = allfires.fires[fid].centroid[0]
        #gdf.loc[fid,'clon'] = allfires.fires[fid].centroid[1]
        if 'mergid' in dd.keys():
            gdf.loc[fid,'mergid'] = fid
        if 'ftype' in dd.keys():
            gdf.loc[fid,'ftype'] = allfires.fires[fid].ftype
        if 'isactive' in dd.keys():
            gdf.loc[fid,'isactive'] = 1
        if 'invalid' in dd.keys():
            gdf.loc[fid,'invalid'] = 0
        if 'isignition' in dd.keys():
            gdf.loc[fid,'isignition'] = allfires.fires[fid].isignition
        if 't_inactive' in dd.keys():
            gdf.loc[fid,'t_inactive'] = int(allfires.fires[fid].t_inactive)
        if 'n_pixels' in dd.keys():
            gdf.loc[fid,'n_pixels'] = int(allfires.fires[fid].n_pixels)
        if 'n_newpixels' in dd.keys():
            gdf.loc[fid,'n_newpixels'] = int(allfires.fires[fid].n_newpixels)
        if 'farea' in dd.keys():
            gdf.loc[fid,'farea'] = allfires.fires[fid].farea
        if 'fperim' in dd.keys():
            gdf.loc[fid,'fperim'] = allfires.fires[fid].fperim
        if 'flinelen' in dd.keys():
            gdf.loc[fid,'flinelen'] = allfires.fires[fid].flinelen
        if 'duration' in dd.keys():
            gdf.loc[fid,'duration'] = allfires.fires[fid].duration
        if 'pixden' in dd.keys():
            gdf.loc[fid,'pixden'] = allfires.fires[fid].pixden
        if 'meanFRP' in dd.keys():
            gdf.loc[fid,'meanFRP'] = allfires.fires[fid].meanFRP
        if 'tst_year' in dd.keys():
            gdf.loc[fid,'tst_year'] = int(allfires.fires[fid].t_st[0])
        if 'tst_month' in dd.keys():
            gdf.loc[fid,'tst_month'] = int(allfires.fires[fid].t_st[1])
        if 'tst_day' in dd.keys():
            gdf.loc[fid,'tst_day'] = int(allfires.fires[fid].t_st[2])
        if 'tst_ampm' in dd.keys():
            gdf.loc[fid,'tst_ampm'] = allfires.fires[fid].t_st[3]
        if 'ted_year' in dd.keys():
            gdf.loc[fid,'ted_year'] = int(allfires.fires[fid].t_ed[0])
        if 'ted_month' in dd.keys():
            gdf.loc[fid,'ted_month'] = int(allfires.fires[fid].t_ed[1])
        if 'ted_day' in dd.keys():
            gdf.loc[fid,'ted_day'] = int(allfires.fires[fid].t_ed[2])
        if 'ted_ampm' in dd.keys():
            gdf.loc[fid,'ted_ampm'] = allfires.fires[fid].t_ed[3]
        if 'ted_doy' in dd.keys():
            gdf.loc[fid,'ted_doy'] = allfires.fires[fid].cdoy

        # change mergeid in case fire has been merged
        # (this is probably only useful for retrospective calculation :
        #  For NRT run: all active fires should not be merged;
        #  maybe this should only be done in sfs recording)
        if 'mergid' in dd.keys():
            if fid in heritage.keys():
                gdf.loc[fid,'mergid'] = heritage[fid]

    # drop entries for newly invalidated fire objects
    if 'invalid' in dd.keys():
        for fid in (allfires.fids_invalid):
            gdf.loc[fid,'invalid'] = 1
        gdf.drop(gdf.index[gdf['invalid'] == 1], inplace = True)

    # # drop inactive fires
    # gdf.drop(gdf.index[gdf['isactive'] == 0], inplace = True)

    # make sure the attribute data formats follow that defined in dd (is this needed for gpkg?)
    for v,tp in dd.items():
        if v != 'fireID':
            gdf[v] = gdf[v].astype(tp)

    # update the hull of each active fire as the geometry column
    for fid in allfires.fids_active:
        fhull = allfires.fires[fid].hull
        if fhull is not None:
            if fhull.geom_type == 'MultiPolygon':
                gdf.loc[fid,'geometry'] = gpd.GeoDataFrame(geometry=[fhull]).geometry.values
            else:
                gdf.loc[fid,'geometry'] = fhull

    return gdf

def make_gdf_fline(allfires, heritage):
    ''' Create geopandas DataFrame for active fire line

    Parameters
    ----------
    allfires : Allfires object
        the Allfires object to be used to create gdf

    Returns
    -------
    gdf : geopandas DataFrame
        the gdf containing half daily fire line
    '''
    import geopandas as gpd

    # initialize the GeoDataFrame
    gdf = gpd.GeoDataFrame(columns=['fireID', 'mergid'],crs='epsg:4326', geometry=[])
    gdf['fireID'] = gdf['fireID'].astype('int')
    gdf['mergid'] = gdf['mergid'].astype('int')

    # for each active fire, record the fline to the geometry column of gdf
    for fid in allfires.fids_active:
        fline = allfires.fires[fid].fline
        if fline is not None:
            gdf.loc[fid,'fireID'] = fid

            # add and correct mergid
            gdf.loc[fid,'mergid'] = fid
            if fid in heritage.keys():
                gdf.loc[fid,'mergid'] = heritage[fid]

            # gdf.loc[fid,'geometry'] = gpd.GeoDataFrame(geometry=[fline]).geometry.values
            if fline.geom_type == 'MultiLineString':
                gdf.loc[fid,'geometry'] = gpd.GeoDataFrame(geometry=[fline]).geometry.values
            else:
                gdf.loc[fid,'geometry'] = fline

    if len(gdf) > 1:
        # dissolve by mergid
        gdf = gdf.dissolve(by = 'mergid')


    return gdf

def make_gdf_NFP(allfires, heritage):
    ''' Create geopandas DataFrame for new fire pixels (detected at this time step)

    Parameters
    ----------
    allfires : Allfires object
        the Allfires object to be used to create gdf

    Returns
    -------
    gdf : geopandas DataFrame
        the gdf containing half daily new fire pixels
    '''
    import geopandas as gpd
    #from FireConsts import dd
    from shapely.geometry import MultiPoint

    # initialize the GeoDataFrame
    gdf = gpd.GeoDataFrame(columns=['fireID', 'mergid'],crs='epsg:4326', geometry=[])
    gdf['fireID'] = gdf['fireID'].astype('int')
    gdf['mergid'] = gdf['mergid'].astype('int')

    # for each fire, record the newlocs to the geometry column of gdf
    for fid in allfires.fids_active:
        newlocs = allfires.fires[fid].newlocs
        if len(newlocs) > 0:
            gdf.loc[fid,'fireID'] = fid
            nfp = MultiPoint([(l[1],l[0]) for l in newlocs])
            gdf.loc[fid,'geometry'] = gpd.GeoDataFrame(geometry=[nfp]).geometry.values
            # add and correct mergid
            gdf.loc[fid,'mergid'] = fid
            if fid in heritage.keys():
                gdf.loc[fid,'mergid'] = heritage[fid]

    if len(gdf) > 1:
        # dissolve by mergid
        gdf = gdf.dissolve(by = 'mergid')

    return gdf

def make_gdf_NFPlist(allfires, heritage):
    ''' Create geopandas DataFrame for new fire pixels (detected at this time step)

    Parameters
    ----------
    allfires : Allfires object
        the Allfires object to be used to create gdf

    Returns
    -------
    gdf : geopandas DataFrame
        the gdf containing half daily new fire pixels
    '''
    import pandas as pd

    # initialize the GeoDataFrame
    # columns=['fireID', 'mergid','lon','lat','line','sample']
    columns=['fireID', 'mergid','lon','lat','frp','origin']
    df = pd.DataFrame(columns=columns)
    fid_df = 0

    # for each fire, record the newlocs to the geometry column of gdf
    for fid in allfires.fids_active:
        mergid = fid
        if fid in heritage.keys():
            mergid = heritage[fid]
        if len(allfires.fires[fid].newlocs)>0:
            lats, lons = zip(*allfires.fires[fid].newlocs)
            # line, sample, frp = zip(*allfires.fires[fid].newpixelatts)
            frp, origin = zip(*allfires.fires[fid].newpixelatts)
            for i,lat in enumerate(lats):
                df.loc[fid_df,'fireID'] = fid_df
                df.loc[fid_df,'mergid'] = mergid
                df.loc[fid_df,'lon'] = lons[i]
                df.loc[fid_df,'lat'] = lat
                # df.loc[fid_df,'line'] = line[i]
                # df.loc[fid_df,'sample'] = sample[i]
                df.loc[fid_df,'frp'] = frp[i]
                df.loc[fid_df,'origin'] = origin[i]
                fid_df += 1

    return df

def save_gdf_1t(allfires, heritage,regnm):
    ''' Creat gdf using one Allfires object and save it to a geopackage file at 1 time step.
            This can be used  for fperim, fline, and NFP files.

    Parameters
    ----------
    allfires : Allfires object
        the Allfires object to be used to create gdf
    '''
    import FireIO

    # create gdf using previous time step gdf values and new allfires object
    gdf_fperim = make_gdf_fperim(allfires, heritage,regnm)
    gdf_fline = make_gdf_fline(allfires, heritage)
    gdf_NFP = make_gdf_NFP(allfires, heritage)
    if len(gdf_fperim) > 0:
        FireIO.save_gpkgobj(gdf_fperim,gdf_fline,gdf_NFP,allfires.t,regnm)


    # gdf_NFPlist = make_gdf_NFPlist(allfires, heritage)
    # if len(gdf_NFPlist) > 0:
    #     FireIO.save_FP_txt(gdf, allfires.t)

def save_gdf_trng(tst,ted,regnm):
    ''' Wrapper to create and save gpkg files for a time period

    Parameters
    ----------
    tst : tuple, (int,int,int,str)
        the year, month, day and 'AM'|'PM' at start time
    ted : tuple, (int,int,int,str)
        the year, month, day and 'AM'|'PM' at end time
    '''
    import FireObj, FireIO

    # loop over all days during the period
    endloop = False  # flag to control the ending of the loop
    t = list(tst)    # t is the time (year,month,day,ampm) for each step
    heritage = dict(FireIO.load_fobj(ted,regnm).heritages)
    while endloop == False:
        print(t)

        # read Allfires object from the saved pkl file
        allfires = FireIO.load_fobj(t,regnm)

        # create and save gdfs according to input options
        save_gdf_1t(allfires,heritage,regnm)

        # time flow control
        #  - if t reaches ted, set endloop to True to stop the loop
        if FireObj.t_dif(t,ted)==0:
            endloop = True

        #  - update t with the next time stamp
        t = FireObj.t_nb(t,nb='next')

if __name__ == "__main__":
    ''' The main code to record daily geojson data for a time period
    '''

    import time
    t1 = time.time()
    # set the start and end time
    tst=(2021,7,13,'AM')
    ted=(2021,9,15,'PM')

    # for each day during the period, create and save geojson summary file
    #save_gdf_trng(tst=tst,ted=ted,fall=True)
    #save_gdf_trng(tst=tst,ted=ted,fperim=True)
    # save_gdf_trng(tst=tst,ted=ted,fline=True)
    save_gdf_trng(tst=tst,ted=ted,regnm='Dixie')

    t2 = time.time()
    print(f'{(t2-t1)/60.} minutes used to run code')

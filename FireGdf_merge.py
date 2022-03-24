""" FireGdf
Module for creating regional geojson summary at each time step

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



def make_gdf_fperim(allfires,heritage,region=''):
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
    from FireConsts import dd
    import FireIO, FireObj

    # initialize the gdf
    t_pt = FireObj.t_nb(allfires.t,nb='previous')
    if FireIO.check_gdfobj(t_pt,region=region) & (allfires.t[0]==t_pt[0]):
        # read gdf at previous day as initial gdf values if data file at same year exists
        gdf = FireIO.load_gdfobj(t_pt,region=region)
        gdf['isactive'] = 0   # pre-set all fires to inactive and update true status below
    else:
        # initialize the GeoDataFrame with fire attribute names (in dd) as columns
        gdf = gpd.GeoDataFrame(columns=list(dd.keys()),crs='epsg:4326', geometry=[])
        gdf = gdf.set_index('fireid')   # set fid column as the index column

    # update data for each active fire object
    for fid in (allfires.fids_active):
        #gdf.loc[fid,'clat'] = allfires.fires[fid].centroid[0]
        #gdf.loc[fid,'clon'] = allfires.fires[fid].centroid[1]
        gdf.loc[fid,'mergid'] = fid
        gdf.loc[fid,'ftype'] = allfires.fires[fid].ftype
        gdf.loc[fid,'isactive'] = 1
        gdf.loc[fid,'invalid'] = 0
        gdf.loc[fid,'isignition'] = allfires.fires[fid].isignition
        gdf.loc[fid,'t_inactive'] = int(allfires.fires[fid].t_inactive)
        gdf.loc[fid,'n_pixels'] = int(allfires.fires[fid].n_pixels)
        gdf.loc[fid,'n_newpixels'] = int(allfires.fires[fid].n_newpixels)
        gdf.loc[fid,'farea'] = allfires.fires[fid].farea
        gdf.loc[fid,'fperim'] = allfires.fires[fid].fperim
        gdf.loc[fid,'flinelen'] = allfires.fires[fid].flinelen
        gdf.loc[fid,'duration'] = allfires.fires[fid].duration
        gdf.loc[fid,'pixden'] = allfires.fires[fid].pixden
        gdf.loc[fid,'meanFRP'] = allfires.fires[fid].meanFRP
        gdf.loc[fid,'tst_year'] = int(allfires.fires[fid].t_st[0])
        gdf.loc[fid,'tst_month'] = int(allfires.fires[fid].t_st[1])
        gdf.loc[fid,'tst_day'] = int(allfires.fires[fid].t_st[2])
        gdf.loc[fid,'tst_ampm'] = allfires.fires[fid].t_st[3]
        gdf.loc[fid,'ted_year'] = int(allfires.fires[fid].t_ed[0])
        gdf.loc[fid,'ted_month'] = int(allfires.fires[fid].t_ed[1])
        gdf.loc[fid,'ted_day'] = int(allfires.fires[fid].t_ed[2])
        gdf.loc[fid,'ted_ampm'] = allfires.fires[fid].t_ed[3]
        gdf.loc[fid,'ted_doy'] = allfires.fires[fid].cdoy

        # change mergeid in case fire has been merged
        if fid in heritage.keys():
            gdf.loc[fid,'mergid'] = heritage[fid]
    
    # drop entries for newly invalidated fire objects
    for fid in (allfires.fids_invalid):
        gdf.loc[fid,'invalid'] = 1
    gdf.drop(gdf.index[gdf['invalid'] == 1], inplace = True)
    
    # drop inactive fires
    gdf.drop(gdf.index[gdf['isactive'] == 0], inplace = True)

    # make sure the attribute data formats follow that defined in dd
    for v,tp in dd.items():
        if v != 'fireid':
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
    #from FireConsts import dd

    # initialize the GeoDataFrame
    gdf = gpd.GeoDataFrame(columns=['fireid', 'mergid'],crs='epsg:4326', geometry=[])

    # for each active fire, record the fline to the geometry column of gdf
    for fid in allfires.fids_active:
        fline = allfires.fires[fid].fline
        if fline is not None:
            gdf.loc[fid,'fireid'] = fid
            
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
    gdf = gpd.GeoDataFrame(columns=['fireid', 'mergid'],crs='epsg:4326', geometry=[])

    # for each fire, record the newlocs to the geometry column of gdf
    for fid in allfires.fids_active:
        gdf.loc[fid,'fireid'] = fid
        newlocs = allfires.fires[fid].newlocs
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
    df = pd.DataFrame(columns=['fireid', 'mergid','lon','lat','line','sample'])
    fid_df = 0

    # for each fire, record the newlocs to the geometry column of gdf
    for fid in allfires.fids_active:
        mergid = fid
        if fid in heritage.keys():
            mergid = heritage[fid]
        if len(allfires.fires[fid].newlocs)>0:
            lats, lons = zip(*allfires.fires[fid].newlocs)
            line, sample, frp = zip(*allfires.fires[fid].newpixelatts)
            for i,lat in enumerate(lats):
                df.loc[fid_df,'fireid'] = fid_df
                df.loc[fid_df,'mergid'] = mergid
                df.loc[fid_df,'lon'] = lons[i]
                df.loc[fid_df,'lat'] = lat
                df.loc[fid_df,'line'] = line[i]
                df.loc[fid_df,'sample'] = sample[i]
                fid_df += 1
            
    return df

def save_gdf_1t(allfires,heritage,fperim=False,fline=False,NFP=False,NFP_txt=False,region=''):
    ''' Creat gdf using one Allfires object and save it to a geojson file at 1 time step.
            This can be used  for fperim, fline, and NFP files.

    Parameters
    ----------
    allfires : Allfires object
        the Allfires object to be used to create gdf
    '''
    import FireIO

    # create gdf using previous time step gdf values and new allfires object
    if fperim:
        gdf = make_gdf_fperim(allfires, heritage,region=region)
        if len(gdf) > 0:
            FireIO.save_gdfobj(gdf,allfires.t,param='',op='',region=region)
    elif fline:
        gdf = make_gdf_fline(allfires, heritage)
        if len(gdf) > 0:
            FireIO.save_gdfobj(gdf,allfires.t,param='',op='FL',region=region)
    elif NFP:
        gdf = make_gdf_NFP(allfires, heritage)
        if len(gdf) > 0:
            FireIO.save_gdfobj(gdf,allfires.t,param='',op='NFP',region=region)
    elif NFP_txt:
        gdf = make_gdf_NFPlist(allfires, heritage)
        FireIO.save_FP_txt(gdf, allfires.t,region=region)
    else:
        return

def save_gdf_trng(tst,ted,fperim=False,fline=False,NFP=False,NFP_txt=False,fall=False,region=''):
    ''' Wrapper to create and save gdf files for a time period

    Parameters
    ----------
    tst : tuple, (int,int,int,str)
        the year, month, day and 'AM'|'PM' at start time
    ted : tuple, (int,int,int,str)
        the year, month, day and 'AM'|'PM' at end time
    fperim : bool
        if set to true, save fire perimeters
    fline : bool
        if set to true, save fire lines
    NFP : bool
        if set to true, save new fire pixels
    fall : bool
        if set to true, save all
    '''
    import FireObj, FireIO

    # if fall is True, set all options True
    if fall:
        fperim,fline,NFP,NFP_txt = True,True,True,True

    # loop over all days during the period
    endloop = False  # flag to control the ending of the loop
    t = list(tst)    # t is the time (year,month,day,ampm) for each step
    heritage = dict(FireIO.load_fobj(ted,region=region).heritages)
    while endloop == False:
        #print(t)

        # read Allfires object from the saved pkl file
        allfires = FireIO.load_fobj(t,region=region)

        # create and save gdfs according to input options
        save_gdf_1t(allfires,heritage,fperim=fperim,fline=fline,NFP=NFP,NFP_txt=NFP_txt,region=region)

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
    tst=(2020,6,1,'AM')
    ted=(2020,8,31,'PM')

    # for each day during the period, create and save geojson summary file
    #save_gdf_trng(tst=tst,ted=ted,fall=True)
    #save_gdf_trng(tst=tst,ted=ted,fperim=True)
    # save_gdf_trng(tst=tst,ted=ted,fline=True)
    save_gdf_trng(tst=tst,ted=ted,NFP=True)

    t2 = time.time()
    print(f'{(t2-t1)/60.} minutes used to run code')

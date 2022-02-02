# -*- coding: utf-8 -*-
"""
Lakes workflow
Created on Wed Jan 19 14:14:52 2022

@author: rebec
"""

def retrieve_water_tile(geom, year):
    '''retrieve the Global Surface Water tile number
    from the bounding box of a geometry
    
    Parameters
    ----------
    geom : geometry
        the geometry of the fire to process
    year : year
        year from which GSW data should be taken

    Returns
    -------
    tiles : tuple of lists
        lists containing all overlapping x and y tiles
    '''
    import math
    import numpy as np
    
    # starting lats and lons for tiles
    tilesx = np.arange(-180, 180, 10) 
    tilesy = np.flip(np.arange(40, 80, 10))
    
    # extract bounding box of fire
    minx, miny, maxx, maxy = geom.bounds
    
    # extract nearest start lat and lon form bounding box
    minx = math.floor(minx/10)*10 
    maxx = math.floor(maxx/10)*10 
    miny = math.floor(miny/10)*10 
    maxy = math.floor(maxy/10)*10 
    
    # find indices
    tilex = np.where(tilesx == minx)[0].tolist()
    tiley = np.where(tilesy == miny)[0].tolist()
    if maxx - minx > 0:
        tilex.append(np.where(tilesx == maxx)[0].tolist()[0])
    if maxy - miny > 0:
        tiley.append(np.where(tilesy == maxy)[0].tolist()[0])
    
    tiles = (tilex, tiley)
    
    lake_files = tile_2_fnm(tiles, year)
    
    return lake_files


def tile_2_fnm(tiles, year):
    '''create filenames of all GSW tiles returned by 
    retrieval_water_tile
    
    Parameters
    ----------
    tiles : tuple of 2 lists
        2 lists containing the tile numbers in x and y direction
    year : year
        year from which GSW data should be taken

    Returns
    -------
    fnms : list of strings
        lists containing paths of all overlapping tile files
    '''
    import itertools, os
    
    tilex, tiley = tiles
    tilex = [str(tile*4).zfill(3) for tile in tilex]
    tiley = [str(tile*4).zfill(2) for tile in tiley]
    all_combinations = [list(zip(each_permutation, tiley)) for each_permutation in itertools.permutations(tilex, len(tiley))]
    
    fnms = []
    basepath = 'D:/fire_atlas/Data/GlobalSurfaceWater/vector/'+str(year)+'/area/GSW_300m_'
    for tile in all_combinations:
        
        fnm = basepath+str(tile[0][1])+'_'+str(tile[0][0])+'.gpkg'
        # check if path exists
        if os.path.exists(fnm):
            fnms.append(fnm)
        else:
           fnms = False
           # this means the tile does not contain any lakes
           # (or has not been processed for years < 2020)
    
    return fnms


def dissolve_lake_geoms(geom, year):
    '''read all lake geometries within the geometry boundaries and 
    if needed dissolve them to one MultiPolygon
    
    Parameters
    ----------
    geom : geometry
        the geometry of the fire to process
    year : year
        year from which GSW data should be taken

    Returns
    -------
    lakediss : geopandas dataframe
        can be empty (no overlap), or contains one entry (index = 1) with a 
        Polygon or MultiPolygon
    '''
    import geopandas as gpd
    import pandas as pd
    
    # retrieve GSW tile paths
    fnm_lakes = retrieve_water_tile(geom, year) # retrieve water tile paths
    
    # load and dissolve water features
    if fnm_lakes:
        if len(fnm_lakes) > 1:
            temp = []
            for lake in fnm_lakes:
                temp.append(gpd.read_file(lake, mask=geom))
            lakes_test = pd.concat(temp, ignore_index = True)
        else:
            lakes_test = gpd.read_file(fnm_lakes[0], mask=geom)
        
        lakes_test = lakes_test.assign(diss = 1)
        lakediss = lakes_test.dissolve(by = 'diss') # takes very long
    else:
        lakediss = None
    
    return lakediss
   

def clip_lakes_1fire_outer(gdf_1d, lakediss, lake_idx):
    '''clip lakes from a fire geometry (only outer lakes)
    
    Parameters
    ----------
    geom : geometry
        the geometry of the fire to process
    lakediss : geopandas dataframe
        can be empty (no overlap), or contains one entry (index = 1) with a 
        Polygon or MultiPolygon

    Returns
    -------
    tiles : tuple of lists
        lists containing all overlapping x and y tiles
    '''
    
    import pyproj
    import FireClustering
    #import copy
    
    
    # initialise
    #geom_nolakes = copy.deepcopy(geom) # geometry without lakes
    geod = pyproj.Geod(ellps='WGS84') # geoid used for computing distances in m
    intsct_length = 0 # length of perimeter bordering lake
    lake_area = 0 # total lake area within the fire
    
    # extract fire geometry
    geom = gdf_1d.iloc[0].geometry
    
    # extract potential lakes in geometry bounding box
    lake_fids = FireClustering.idx_intersection(lake_idx, geom.bounds)

    # loop over lakes
    for i in lake_fids:
        lake = lakediss[i]
        if lake.intersects(geom):
            # if lakes are within fire scar we just compute lake statistics
            if lake.within(geom):
                temp = geod.geometry_area_perimeter(lake)
                lake_area += abs(temp[0])
                intsct_length += temp[1] # should this be included?
            # if lakes are intersecting fire scar we clip
            else:
                geom = geom.difference(lake)
                true_intsct = lake.intersection(geom)
                intsct_length += geod.geometry_area_perimeter(true_intsct)[1]
    
    # update gdf
    gdf_1d = gdf_1d.assign(lake_area = lake_area)
    gdf_1d = gdf_1d.assign(lake_border = intsct_length)
    gdf_1d.set_geometry([geom], inplace = True)
    
    return gdf_1d

def clip_lakes_1fire(geom, lakediss):
    '''clip lakes from a fire geometry
    
    Parameters
    ----------
    geom : geometry
        the geometry of the fire to process
    lakediss : geopandas dataframe
        can be empty (no overlap), or contains one entry (index = 1) with a 
        Polygon or MultiPolygon

    Returns
    -------
    tiles : tuple of lists
        lists containing all overlapping x and y tiles
    '''
    
    import pyproj
    import copy
    
    geom_nolakes = copy.deepcopy(geom) # geometry without lakes
    geod = pyproj.Geod(ellps='WGS84') # geoid used for computing distances in m
    intsct_length = 0 # length of perimeter bordering lake
    lake_area = 0 # total lake area within the fire
    
    # clip lakes
    geom = geom.difference(lakediss.loc[1]['geometry'])

    # additional attributes
    intsct = geom_nolakes.difference(geom) # linestring of intersection
    
    # if not intsct.is_empty: # if there are lakes in the area
    if intsct.geom_type == 'Polygon': # if only one lake is present we put it in a list for the loop
        intsct = [intsct]
    # check which lakes are inside firescar and which are outside
    for i in range(len(intsct)):
        lake = intsct[i]
        if lake.within(geom_nolakes): # if lakes are within fire scar
            temp = geod.geometry_area_perimeter(lake)
            lake_area += abs(temp[0])
            intsct_length += temp[1]
        else: # if lakes are bordering fire scar
            true_intsct = lake.intersection(geom)
            intsct_length += geod.geometry_area_perimeter(true_intsct)[1]
    # else:
    #     print("surprise! this also doesn't contain lakes")

    lakes = lakediss.loc[1]['geometry'].difference(geom).intersection(geom_nolakes.buffer(0.0001))
    out = (geom, lakes, lake_area, intsct_length)

    return out


def clip_lakes_final(gdf, t):
    '''Wrapper to clip the lakes form final polygons'''
    import geopandas as gpd
    
    # extract year
    year = t[0]
    gdf_new = []
    lake_list = []
    # loop over all final perimeters
    for row in gdf.itertuples():
        fid = row.mergid
        print(fid)
        geom = row.geometry
        
        # read and dissolve lakes for the geom
        lakediss = dissolve_lake_geoms(geom, year) 

        if isinstance(lakediss, type(None)):
            #print('GSW tile missing or does not contain lakes!')
            gdf_new.append((fid, geom, 0, 0))
        else:
            if len(lakediss) > 0:
                clipped_geoms,lake_geom,lake_area,intsct_length = clip_lakes_1fire(geom, lakediss)
                gdf_new.append((fid, clipped_geoms, lake_area, intsct_length))
                lake_list.append((fid, lake_geom))
            else:
                #print('fire '+str(fid)+' does not contain lakes.')
                gdf_new.append((fid, geom, 0, 0))
    
    # turn fire geom and attributes into geopandas dataframe
    fids, geoms, lake_areas, lake_intscts = zip(*gdf_new)
    d = {'fid': fids, 'geometry': geoms, 'lake_area': lake_areas, 'lake_border': lake_intscts}
    gdf_new = gpd.GeoDataFrame(d, crs="EPSG:4326")
    gdf_new = gdf_new.assign(mergid = gdf_new['fid'])
    
    # same for lake geometries
    fids, geoms = zip(*lake_list)
    d = {'fid': fids, 'geometry': geoms, 'mergid': fids}
    gdf_lakes = gpd.GeoDataFrame(d, crs="EPSG:4326")
    
    #gdf_new = pd.DataFrame(gdf_new, columns=[fid,clipped_geoms])
    #gdf_new = gpd.GeoDataFrame(gdf_new, crs='epsg:4326', geometry=clipped_geoms)
    return gdf_new, gdf_lakes



def sort_by_tile():
    """ sort the Fire perimeters of an Allfires object by tile number"""
    pass


def process_1t():
    pass

if __name__ == "__main__":
    ''' The main code to record daily geojson data for a time period
    '''

    import time
    t1 = time.time()
    # set the start and end time
    tst=(2020,6,1,'AM')
    ted=(2020,8,31,'PM')

    # t = ted
    # allfires = FireIO.load_fobj(t)
    # tilelist = create_tile_dict(allfires)
    # fireids,tileids = zip(*tilelist)
    # tilesx
    
    # save_gdf_trng(tst=tst,ted=ted,fperim=True)

    t2 = time.time()
    print(f'{(t2-t1)/60.} minutes used to run code')
    
    # create dict of fires and tiles
    # fire = allfires.fires[0]
    # pick tile 1
    # process all fires
    # check if any fire stretches over 
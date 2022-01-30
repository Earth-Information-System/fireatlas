# -*- coding: utf-8 -*-
"""
Lakes workflow
Created on Wed Jan 19 14:14:52 2022

@author: rebec
"""

def retrieve_water_tile(geom):
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
    
    return (tilex, tiley)

def clip_lakes_1fire(geom, fnm_lakes):
    
    import geopandas as gpd
    
    for lake in fnm_lakes:
        lakes_test = gpd.read_file(lake, mask=geom)
        #lakes_test['geometry'] = lakes_test.buffer(0) # fix invalid geometries
        lakes_test = lakes_test.assign(diss = 1)
        lakediss = lakes_test.dissolve(by = 'diss') # takes very long

        # clip lakes
        if len(lakediss) == 1:
            geom = geom.difference(lakediss.loc[1]['geometry'])
    # envgdf = gpd.GeoDataFrame(geometry=gpd.GeoSeries(test_clip))
    # fnm = 'D:/fire_atlas/2_pipeline/2020/testclip.gpkg'
    # envgdf.to_file(fnm, driver='GPKG')
    return geom

def clip_lakes(gdf, t):
    '''Wrapper to clip the lakes form final polygons'''
    import PostProcess
    import copy
    import pandas as pd
    import geopandas as gpd
    
    # extract year
    year = t[0]
    gdf_new = []
    # loop over all final perimeters
    for row in gdf.itertuples():
        fid = row.mergid
        print(fid)
        geom = row.geometry
        tiles = PostProcess.retrieve_water_tile(geom)
        lake_files = PostProcess.tile_2_fnm(tiles, year)
        if lake_files:
            clipped_geoms = clip_lakes_1fire(geom, lake_files)
            gdf_new.append((fid,clipped_geoms))
        else:
            print(tiles+'for fid'+fid+' not processed.')
    
    fids, geoms = zip(*gdf_new)
    d = {'fid': fids, 'geometry': geoms}
    gdf_new = gpd.GeoDataFrame(d, crs="EPSG:4326")
    
    #gdf_new = pd.DataFrame(gdf_new, columns=[fid,clipped_geoms])
    #gdf_new = gpd.GeoDataFrame(gdf_new, crs='epsg:4326', geometry=clipped_geoms)
    return gdf_new



def tile_2_fnm(tiles, year):
    import itertools, os
    
    tilex, tiley = tiles
    tilex = [tile*4 for tile in tilex]
    tiley = [tile*4 for tile in tiley]
    all_combinations = [list(zip(each_permutation, tiley)) for each_permutation in itertools.permutations(tilex, len(tiley))]
    
    fnms = []
    basepath = 'D:/fire_atlas/Data/GlobalSurfaceWater/vector/'+str(year)+'/tile_'
    for tile in all_combinations:
        fnm = basepath+str(tile[0][1])+'_'+str(tile[0][0])+'.gpkg'
        # check if path exists
        if os.path.exists(fnm):
            fnms.append(fnm)
        else:
           fnms = False
           break
    
    return fnms

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
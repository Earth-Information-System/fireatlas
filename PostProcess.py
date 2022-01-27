# -*- coding: utf-8 -*-
"""
Lakes workflow
Created on Wed Jan 19 14:14:52 2022

@author: rebec
"""

def retrieve_water_tile(fire):
    import math
    import numpy as np
    
    # starting lats and lons for tiles
    tilesx = np.arange(-180, 180, 10) 
    tilesy = np.flip(np.arange(40, 80, 10))
    
    # extract bounding box of fire
    minx, miny, maxx, maxy = fire.hull.bounds
    
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
    
    return tilex, tiley

def create_tile_dict(allfires):
    tiledict = {}
    for fire_no in range(allfires.number_of_fires):
        fire = allfires.fires[fire_no]
        tiledict[fire_no] = list(retrieve_water_tile(fire))
    tiledict = list(tiledict.items())
    return tiledict

def sort_by_tile():
    """ sort the Fire perimeters of an Allfires object by tile number"""
    pass


def process_1t()

if __name__ == "__main__":
    ''' The main code to record daily geojson data for a time period
    '''

    import time
    t1 = time.time()
    # set the start and end time
    tst=(2020,6,1,'AM')
    ted=(2020,8,31,'PM')

    t = ted
    allfires = FireIO.load_fobj(t)
    tilelist = create_tile_dict(allfires)
    fireids,tileids = zip(*tilelist)
    tilesx
    
    # save_gdf_trng(tst=tst,ted=ted,fperim=True)

    t2 = time.time()
    print(f'{(t2-t1)/60.} minutes used to run code')
    
    # create dict of fires and tiles
    fire = allfires.fires[0]
    # pick tile 1
    # process all fires
    # check if any fire stretches over 
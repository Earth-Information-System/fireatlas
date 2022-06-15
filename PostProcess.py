# -*- coding: utf-8 -*-
"""
Lakes workflow
Created on Wed Jan 19 14:14:52 2022

@author: rebec
"""
def createProj(projWkt):
    '''creates an osr spatial reference from wkt
    *** not in use ***'''
    import osr
    
    src_srs = osr.SpatialReference()                                    # create the spatial reference
    src_srs.ImportFromWkt(projWkt)
    
    return src_srs

def polygonise_tif(fileout, proj, outband, outmask):
    '''polygonises a raster
    *** not in use ***'''
    
    import ogr, gdal
    # create output vector datasource in memory
    driver = ogr.GetDriverByName('GPKG')
    temp_ds = driver.CreateDataSource(fileout)

    # create temp dataset in memory
    temp_layer = temp_ds.CreateLayer('temp', srs = createProj(proj), geom_type = ogr.wkbPolygon)

    # store pixel value in first field (right now in int, could be changed to dbl) 
    gdal.Polygonize(outband, outmask, temp_layer, 0, [], callback=None )

def hull_from_pixels(coords_geo):
    '''
    Clusters a list of geographic coordinates and computes the convex hull
    of each cluster
    
    Parameters
    ----------
    coords_geo: list of tuples,
        geographic coordinates
    
    Returns
    -------
    gdf : geopandas geodataframe
        contains a convex hull for each cluster
    '''
    
    import FireClustering, FireVector
    import geopandas as gpd
    
    # cluster the pixels
    cid = FireClustering.do_clustering(coords_geo,1)
    
    # compute the hull of each cluster
    hulls = {}
    for fid in range(max(cid)+1):
        locs = [coords_geo[i] for i,x in enumerate(cid) if x == fid]
        hulls[fid] = FireVector.cal_hull(locs, 'none')
    
    # write out as geopandas and to file
    gdf = gpd.GeoDataFrame({'fireid': list(hulls.keys()), 'geometry': list(hulls.values())}, 
                           crs="EPSG:4326")
    gdf = gdf.to_crs('EPSG:3571')
    
    return gdf

def grid_viirs(year, region, geotrans, out_arr, clip=True):
    '''grid VIIRS pixels into polar LAEA:
        1) reads viirs burned area coordinates of a year
        2) transforms lat lon viirs coordinates to NP LAEA
        3) transforms NP LAEA to pixel coordinates of out_arr using geotransform
        4) adds burn pixels to the given numpy array
        optionally: clips the array to the extent of the viirs coordinates
    
    Parameters
    ----------
    year: int,
        year of processing
    geotrans: tuple (6 x float)
        geotransformation tuple of out_arr array from gdal
    out_arr: np.array
        circumpolar LAEA grid in 500m resolution
    clip: bool,
        should the output array and geotrans be clipped to the extent of 
        the available VIIRS pixels? default is True

    Returns
    -------
    geotrans : tuple (6 x float)
        geotransformation tuple of the returned array
    out_arr: numpy array,
        updated numpy array containing viirs active fire detections
    '''
    
    import FireIO
    import glob
    import numpy as np
    import pandas as pd
    
    path = FireIO.get_path(str(year), region)
    filelist_viirs = glob.glob(path+'/Snapshot/*_NFP.txt')
    df = pd.concat([pd.read_csv(i, dtype = {'lat':np.float64, 'lon':np.float64}) 
                    for i in filelist_viirs], ignore_index=True)

    # read lats, lons and pixel dimensions into separate lists
    lats = np.array(df['lat'])
    lons = np.array(df['lon'])
    
    # transform lat lon to polar projected
    x_proj, y_proj = FireIO.geo_to_polar(lons, lats)
    
    # optional: adjust geotrans to image
    if clip:
        geotrans = list(geotrans)
        geotrans[0] = min(x_proj)
        geotrans[3] = max(y_proj)
        
    # transform polar projected to pixel coordinates
    x_px, y_px = FireIO.world2Pixel(geotrans, x_proj, y_proj)
    coords = set(zip(x_px, y_px))
    
    # clip image extent 
    if clip:
        out_arr = np.zeros((max(y_px)+1, max(x_px)+1))
    
    # loop through coords and put in grid
    for coord in coords:
        out_arr[coord[1], coord[0]] = 1
     
    return geotrans, out_arr

def create_burnraster(year,region,reload_mcd64=True,fill_npixel=1,write=True):
    '''creates a binary burned/unburned raster in 500m resolution by merging
    mcd64 burned area and VIIRS active fires
    
    Parameters
    ----------
    year: int,
        processing year to take mcd64 and viirs data from
    reload_mcd64: bool,
        should mcd64 data be included in raster
    fill_npixel: int,
        for holefilling, what number of pixels should be considered
    
    Returns
    -------
    geotrans_out : geopandas geodataframe
        geotransformation tuple of the returned array
    arr_year: numpy array,
        500m-array containing viirs active fire and modis burned area
    '''
    
    from FireConsts import mcd64dir
    import FireIO
    import numpy as np
    import gdal
    from skimage import morphology

    # read example grid for mapping
    fnm = mcd64dir + 'mcd64_'+str(year)+'.tif'
    sample_ds = gdal.Open(fnm)
    geotrans = sample_ds.GetGeoTransform()
    proj = sample_ds.GetProjection()
    cols = sample_ds.RasterXSize
    rows = sample_ds.RasterYSize
    arr_year = np.zeros((cols, rows)).astype(int)

    # grid viirs pixels
    geotrans_new, arr_year = grid_viirs(year, region, geotrans, arr_year)
    
    # add mcd64 burn pixels to clipped grid
    if reload_mcd64:
        xoff, yoff = FireIO.world2Pixel(geotrans, geotrans_new[0], geotrans_new[3])
        mcd64_arr = FireIO.load_burnraster(year,region,'mcd64',int(xoff),int(yoff),
                               xsize=arr_year.shape[1],ysize=arr_year.shape[0])
        arr_year = ((arr_year + mcd64_arr) > 0)*1
    geotrans_out = geotrans_new
    
    # remove holes (of size fill_npixels) from burn scars
    if fill_npixel > 0:
        arr_rev = arr_year == 0
        arr_rev = morphology.remove_small_objects(arr_rev, fill_npixel+1)*1
        arr_year = 1-arr_rev
    
    if write:
        outpath = FireIO.get_path(str(year),region)
        # save to file 
        filename_year = outpath + '/Summary/viirs_' + str(year) + '.tif'
        FireIO.save2gtif(arr_year, filename_year, arr_year.shape[1], arr_year.shape[0], geotrans_out, proj)
    
    return geotrans_out, arr_year

def compute_unburned_area(year, region, gdf):
    '''Computes the unburned island from a fire scar'''
    
    import numpy as np
    import rasterio.mask
    import FireIO
    
    # create new attributes in gdf
    gdf = gdf.assign(burn_px = 0)
    gdf = gdf.assign(unburn_px = 0)
    
    # create or load the burn raster
    # if t == tst:
    geotrans, arr_year = create_burnraster(year,region,reload_mcd64=True,fill_npixel=2,write=True)
    burnfile = FireIO.load_burnraster(year,region,sensor='viirs')
    
    # loop through perimeters and clip
    for fid in gdf.index:
        geom = gdf.loc[fid].geometry
        if not geom.is_empty:
            with rasterio.open(burnfile) as src:
                out_image, out_transform = rasterio.mask.mask(src, [geom], crop=True)
            gdf.loc[fid,'burn_px'] = np.sum(out_image)
            gdf.loc[fid,'unburn_px'] = np.sum(out_image == 0)
    
    return gdf
        
    

def create_burnpoly(t, region,method='hull'):
    '''
    Creates polygons of unburned islands using VIIRS AF and MCD64
    1) create a burned/unburned raster for the year
    2) identify unburned islands using hole-filling
    3) clean data (remove holes in holes, remove small objects)
    4) extract pixel locations of the unburned islands and transform to geographic coordinates
    5) cluster and compute convex hull
    6) save as geopackage
    
    Parameters
    ----------
    t: tuple,
        time tuple
    method: string,
        method for polygonisation, either hull (convex hull, default)
        or pixelwise (all pixels are being polygonised)
    
    '''
    
    import FireIO
    import numpy as np
    from scipy import ndimage
    from skimage import morphology
    
    year = t[0]
    
    # create the burn raster
    geotrans, arr_year = create_burnraster(year,region,reload_mcd64=True,fill_npixel=2,write=True)
    
    # extract unburned islands
    arr_noholes = ndimage.binary_fill_holes(arr_year == 1).astype(int)  # fill all holes in burns
    arr_holes = arr_noholes-arr_year                                    # subtract original
    arr_holes = ndimage.binary_fill_holes(arr_holes == 1)               # fill burn pixels within holes
    arr_holes = morphology.remove_small_objects(arr_holes, 10)          # remove small unburned islands
    
    # now we extract the pixel locations using the geotrans
    coords = []
    index = np.where(arr_holes)
    for xpx, ypx in zip(index[0], index[1]):
        x, y = FireIO.pixel2World(geotrans, ypx, xpx)
        coords.append((x,y))
    x, y = zip(*coords)
    coords_geo = FireIO.polar_to_geo(x,y)
    coords_geo = list(zip(coords_geo[1], coords_geo[0])) # careful! this is lat,lon!!
    
    if method == 'pixelwise':
        pass
    else:
        # do clustering and computation of hull
        gdf = hull_from_pixels(coords_geo)
    
    # write to file
    FireIO.save_gdfobj(gdf,t,param='unburned')        

def retrieve_water_tile(ext, year):
    '''retrieve the Global Surface Water tile number
    from the bounding box of a geometry (in lat lon!)
    
    Parameters
    ----------
    ext : tuple
        the extent of the geometry of the fire to process
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
    minx, miny, maxx, maxy = ext
    
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
    from FireConsts import lakedir
    import itertools, os
    
    # there is no surface water informaion for 2021, so we replace that with 2020
    if year == 2021:
        year = 2020
    if year == 2012:
        year = 2013
    
    tilex, tiley = tiles
    tilex = [str(tile*4).zfill(3) for tile in tilex]
    tiley = [str(tile*4).zfill(2) for tile in tiley]
    all_combinations = [list(zip(each_permutation, tiley)) for each_permutation in itertools.permutations(tilex, len(tiley))]
    
    fnms = []
    basepath = lakedir+str(year)+'/GSW_300m_'
    for tile in all_combinations:
        
        fnm = basepath+str(tile[0][1])+'_'+str(tile[0][0])+'.gpkg'
        # check if path exists
        if os.path.exists(fnm):
            fnms.append(fnm)
        # missing file means the tile does not contain any lakes
        # (or has not been processed! make sure to run poygonise_water_test on the year before)
    
    return fnms

def dissolve_lake_geoms(geom, ext, year):
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
    from shapely.geometry import MultiPolygon, Polygon
    
    # retrieve GSW tile paths
    fnm_lakes = retrieve_water_tile(ext, year) # retrieve water tile paths
    
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
        lakediss = lakes_test.dissolve(by = 'diss')
        
        # at the edges of files we can have artificial holes in the lakes due to mosaicking
        if len(fnm_lakes) > 1:
            if len(lakediss) > 0:
                geom_orig = lakediss.loc[1].geometry
                # cnt = 0
                # for p in geom_orig:
                #     if len(p.interiors)>0:
                #         print(cnt)
                #         print(p)
                #     cnt+=1
                if type(geom_orig) == Polygon:
                    geom_fix = Polygon(geom_orig.exterior)
                else:
                    geom_fix = MultiPolygon(Polygon(p.exterior) for p in geom_orig)
                lakediss.set_geometry([geom_fix], inplace = True)
            
    else:
        lakediss = None
    
    return lakediss
   
def final_lakes(gdf, t, region):
    '''Wrapper to creade a gdf of final lakes
    
    Parameters
    ----------
    gdf : pandas geodataframe
        geodataframe of final fire perimeters
    t : tuple
        time tuple
    
    Returns
    -------
    gdf_lakes: geopandas geodataframe
        geodataframe containing only the lakes fore each fire perimeter id
        (this is used for clipping daily large fire perimeters)
    '''
    
    import geopandas as gpd
    from shapely.geometry import box
    import FireIO
    
    # extract year
    year = t[0]
    lake_list = []
    
    # convert to projected coordinate system (EPSG:3571)
    gdf_3571 = gdf.to_crs('EPSG:3571')
    
    # loop over all final perimeters
    for row in gdf_3571.itertuples():
        fireid = row.mergid
        # print(fid)
        geom = row.geometry
        
        # fetch bounding box in lat lon
        ext = gdf.loc[row.Index].geometry.bounds
        
        # read and dissolve lakes for the geom
        lakediss = dissolve_lake_geoms(geom, ext, year) 
        
        if not isinstance(lakediss, type(None)):
            if len(lakediss) > 0:
                lake_geom = box(*geom.bounds).buffer(2000).intersection(lakediss.loc[1].geometry)
                lake_list.append((fireid, lake_geom))
    
    # turn lake geom and fire id into geopandas dataframe
    fids, geoms = zip(*lake_list)
    d = {'fireid': fids, 'geometry': geoms, 'mergid': fids}
    gdf_lakes = gpd.GeoDataFrame(d, crs="EPSG:3571")
    
    # save to file
    FireIO.save_gdfobj(gdf_lakes,t,param='lakes',region=region)
    
    return gdf_lakes

def clip_lakes_1fire_outer(geom, lakediss):
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
    geom: clipped geometry (only outer lakes clipped)
    lake_area: area within fire covered by lakes
    lake_no: number of lakes within fire
    intsct_length: length of fire perimeter that ends at lake
    intsct_length: total length of fire perimeter bordering lake areas
    '''
    import shapely
    import FireClustering
    
    # initialise
    intsct_length = 0       # length of perimeter bordering lake for this time step
    intsct_length_tot = 0   # total length of intersection incl inner lakes
    lake_area = 0           # total lake area within the fire
    lake_no = 0             # number of lakes within perimeter
    
    # if its a polygon and not multipolygon we put it in a list
    if not isinstance(lakediss, shapely.geometry.multipolygon.MultiPolygon):
        lakediss = shapely.geometry.multipolygon.MultiPolygon([lakediss])
    
    # build an rtree index for fast seach of lakes
    lake_idx = FireClustering.build_rtree(lakediss)
    
    # extract potential lakes in geometry bounding box
    lake_fids = FireClustering.idx_intersection(lake_idx, geom.bounds)
    
    # loop over lakes
    for i in lake_fids:
        lake = lakediss[i]
        if lake.intersects(geom):
            # if lakes are within fire scar we just compute lake statistics
            if lake.within(geom):
                lake_no += 1
                lake_area += lake.area
                intsct_length_tot += lake.length
            # if lakes are intersecting fire scar we clip
            else:
                geom = geom.difference(lake)
                outer_perim = lake.intersection(geom).length
                intsct_length += outer_perim
                intsct_length_tot += outer_perim
    
    return geom, lake_area, lake_no, intsct_length, intsct_length_tot

def clip_lakes_1fire(geom, lakediss):
    '''clip lakes from a fire geometry
    
    Parameters
    ----------
    geom : geometry
        the geometry of the fire to process
    lakediss : geopandas dataframe
        can be empty (no overlap), or contains one entry (index = 1) with a 
        Polygon or MultiPolygon
    lakes_only: bool,
        return only lake r
    
    Returns
    -------
    tiles : tuple of lists
        lists containing all overlapping x and y tiles
    '''
    # import pyproj
    import copy
    
    # clip lakes
    geom = geom.difference(lakediss)
    
    geom_nolakes = copy.deepcopy(geom) # geometry without lakes
    #geod = pyproj.Geod(ellps='WGS84') # geoid used for computing distances in m (if EPSG:4326)
    intsct_length = 0 # length of perimeter bordering lake
    lake_area = 0 # total lake area within the fire
    
    # additional attributes
    intsct = geom_nolakes.difference(geom) # linestring of intersection
    
    # if not intsct.is_empty: # if there are lakes in the area
    if intsct.geom_type == 'Polygon': # if only one lake is present we put it in a list for the loop
        intsct = [intsct]
    # check which lakes are inside firescar and which are outside
    for i in range(len(intsct)):
        lake = intsct[i]
        if lake.within(geom_nolakes): # if lakes are within fire scar
            #temp = geod.geometry_area_perimeter(lake)
            #lake_area += abs(temp[0])
            # intsct_length += temp[1]
            lake_area = lake.area
            intsct_length += lake.length
        else: # if lakes are bordering fire scar
            true_intsct = lake.intersection(geom)
            intsct_length += true_intsct.length
            # intsct_length += geod.geometry_area_perimeter(true_intsct)[1]
    
    out = (geom, lake_area, intsct_length)
    
    return out

def sort_by_tile():
    """ sort the Fire perimeters of an Allfires object by tile number
    *** for circumpolar processing, currently not in use """
    pass

if __name__ == "__main__":
    ''' The main code to record daily geojson data for a time period
    '''
    import sys
    sys.path.insert(1, 'D:/fire_atlas/1_code/fireatlas')
    import os
    if 'GDAL_DATA' not in os.environ:
        os.environ['GDAL_DATA'] = r'C:/Users/rebec/anaconda3/envs/py3work/Library/share/gdal' 
        os.environ['PROJ_LIB'] = r'C:/Users/rebec/anaconda3/envs/fireatlas/Library/share/proj' 
    import time
    t1 = time.time()
    # set the start and end time
    tst=(2020,6,1,'AM')
    ted=(2020,9,30,'PM')
    
    create_burnpoly(ted)

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
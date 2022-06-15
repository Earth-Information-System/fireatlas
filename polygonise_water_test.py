# -*- coding: utf-8 -*-
"""
Transform JRC Global Surface Water tiles to geopackage

Created on Wed Jan 19 13:42:05 2022

@author: rebec
"""

import sys
import os
sys.path.insert(1, 'D:/fire_atlas/1_code/fireatlas')
if 'GDAL_DATA' not in os.environ:
    os.environ['GDAL_DATA'] = r'C:/Users/rebec/anaconda3/envs/py3work/Library/share/gdal' 
    os.environ['PROJ_LIB'] = r'C:/Users/rebec/anaconda3/envs/fireatlas/Library/share/proj' 

import glob
import numpy as np
import gdal#, ogr, osr
import time
import gc
#import matplotlib.pyplot as plt


# functions

def index_from_slice(all_labels, start_end, fid):
    '''creates indices of pixels with an id using the slices returned by 
    scipy's find_objects. This allows for fast identification of indices
    with very large arrays
    
    Parameters
    ----------
    all_labels : numpy array
        array containing all labelled objects
    start_end : slice
        a slice containing start and end x and y coordinates for the object of interest
    fid: int,
        id/label of the object of interest
    
    Returns
    -------
    [ind_x, ind_y] : list of lists
        x and y pixel coordinates of all pixels labelled with the fid of interest
        
    '''
    mask = all_labels[start_end] == fid
    ind_x = np.arange(start_end[1].start,start_end[1].stop)
    ind_y = np.arange(start_end[0].start,start_end[0].stop)
    grid = np.meshgrid(ind_x, ind_y)
    ind_y = grid[0][mask]
    ind_x = grid[1][mask]
    # for i in range(len(ind_x)):
    #     print(arr[ind_x[i], ind_y[i]])
    return [ind_x, ind_y]

def find_edges(arr):
    '''find the boundaries of all features to speed up the calculation of hulls'''
    
    from skimage import segmentation
    
    arr_bin = arr > 0 # create a binary array for the edge detection
    
    inner_edges = segmentation.find_boundaries(arr_bin, mode='inner',background=0)
    # single layer of pixels is not enough, creates splintered hulls!
    inner_edges = segmentation.find_boundaries(inner_edges, mode='thick',background=0)
    inner_edges2 = segmentation.find_boundaries(inner_edges, mode='thick',background=0)
    inner_edges = np.logical_or(inner_edges, inner_edges2)
    inner_edges = np.logical_and(arr_bin, inner_edges)
    
    # fill non-edges with zeros in original array
    arr = np.ma.masked_where(~inner_edges, arr)
    arr = np.ma.filled(arr, fill_value = 0)
    
    return arr

def hull_from_pixels(arr, geoTrans):
    ''' Computes convex hull from outper pixels of each lake object
    
    This version uses ndimage.label for seqmentation (returns more features than measure.label)
    and ndimage.find_objects to extract indices (quite fast)
    Uses only the outer pixels of each lake to speed up the calculation of the convex hull
    Processing steps include:
        1) remove small lakes (< 5 90m pixels)
        2) fill up islands in the lakes
        3) label the lakes (give each lake a unique id)
        4) identify outer pixels of each lake, set inner pixels to zero
        5) loop through all lake ids, transform pixel coordinates to geographic
            coordinates and compute the convex hull of the outer lake pixels
        6) reproject the lakes to a projected coordinate system (North Pole LAEA)
    
    Parameters
    ----------
    arr : numpy array
        bool, lake or no lake
    geoTrans : tuple
        affine georeferencing transform, used to extract pixel coordinates

    Returns
    -------
    gdf : geopandas geodataframe
        gdf containing convex hulls of all lakes
    
    '''
    import FireVector, FireIO
    import geopandas as gpd
    from scipy import ndimage#, sparse
    from skimage import morphology#, measure
    # from shapely.geometry import MultiPolygon
    from shapely.ops import unary_union
    
    t01 = time.time()
    # remove small lakes pixels (< 36 30m pixels)
    arr = morphology.remove_small_objects(arr, 5)
    # fill holes (lake islands) --> helps with segementation and boundary extraction
    arr_lakes = ndimage.binary_fill_holes(arr)
    
    # we keep large islands since they can also burn (see NWT slave lake 2014)
    islands = arr_lakes.astype(int) - arr.astype(int)
    islands = morphology.remove_small_objects(islands == 1, 251)
    # we also keep lakes on islands if they fit the size requirements
    island_lakes = ndimage.binary_fill_holes(islands)
    island_lakes = island_lakes.astype(int) - islands.astype(int)
    island_lakes = morphology.remove_small_objects(island_lakes == 1, 5)
    
    # cluster the pixels: all neighbouring pixel receive one id
    all_labels, maxlab = ndimage.label(arr_lakes)
    all_labels = find_edges(all_labels) # reduce to edges only for speed
    
    # extract the lake labels of lakes containing islands
    island_labels = np.unique(all_labels[islands])
    
    # label the islands and island-lakes
    islands_labels, maxlab_islands = ndimage.label(islands)
    islands_labels = find_edges(islands_labels) # reduce to edges only for speed
    island_lakes_labels, maxlab_island_lakes = ndimage.label(island_lakes)
    
    # use scipy find_objects!
    start_end_idx = ndimage.find_objects(all_labels, max_label=maxlab) # index for lakes
    start_end_idx_il = ndimage.find_objects(island_lakes_labels, max_label=maxlab_island_lakes) # index for lakes on islands
    t0 = time.time()
    print('Finished raster operations.', round((t0-t01)/60, 2))
    print(maxlab, 'lakes to be processed. Estimated time:', int(all_labels.max()/1000*110/60)+1, 'min')
    
    # loop through clusters and extract hulls
    hulls = {}
    # extract lat/lon coordinates from projected image using geotrans
    for fid in range(1,maxlab+1):
        coords = []
        start_end = start_end_idx[fid-1]
        index = index_from_slice(all_labels, start_end, fid)
        if len(index[0]) < 1:
            print('no indices')
            break
        for xpx, ypx in zip(index[0]+0.5, index[1]+0.5): # add 0.5 for center of pixel
            x, y = FireIO.pixel2World(geoTrans, ypx, xpx)
            coords.append((x,y))
        x, y = zip(*coords)
        coords_geo = FireIO.polar_to_geo(x,y)
        coords_geo = list(zip(coords_geo[1], coords_geo[0])) # careful! this is lat,lon!!
        
        # compute the hull of each cluster
        hulls[fid] = FireVector.cal_hull(coords_geo, 'none')
        
        # include islands in lakes if there are any
        if fid in island_labels:
            coords = []
            index = index_from_slice(islands.astype(int), start_end, 1) # technically there could be islands here that are not in the lake!
            for xpx, ypx in zip(index[0], index[1]): 
                x, y = FireIO.pixel2World(geoTrans, ypx+0.5, xpx+0.5)# add 0.5 for center of pixel
                coords.append((x,y,islands_labels[xpx,ypx]))
            x, y, isl_lab = zip(*coords)
            coords_geo = FireIO.polar_to_geo(x,y)
            coords_geo = list(zip(coords_geo[1], coords_geo[0])) # careful! this is lat,lon!!
            
            # compute the hull of each island
            island_hulls = []
            for island_id in np.unique(isl_lab):
                coords_island = [coord for i,coord in enumerate(coords_geo) if isl_lab[i] == island_id]
                island_hulls.append(FireVector.cal_hull(coords_island, 'none'))
                
            # grab the lake geometry and clip the islands
            if len(island_hulls) > 1:
                island_hulls = unary_union(island_hulls)
            else:
                island_hulls = island_hulls[0]
            hulls[fid] = hulls[fid].difference(island_hulls)
        
    # finally: include lakes on islands
    print(maxlab_island_lakes, 'lakes on islands.')
    for il_fid in range(1,maxlab_island_lakes+1):
        coords = []
        start_end = start_end_idx_il[il_fid-1]
        index = index_from_slice(island_lakes_labels, start_end, il_fid) # technically there could be islands here that are not in the lake!
        for xpx, ypx in zip(index[0], index[1]): 
            x, y = FireIO.pixel2World(geoTrans, ypx+0.5, xpx+0.5)# add 0.5 for center of pixel
            coords.append((x,y,islands_labels[xpx,ypx]))
        x, y, isl_lab = zip(*coords)
        coords_geo = FireIO.polar_to_geo(x,y)
        coords_geo = list(zip(coords_geo[1], coords_geo[0])) # careful! this is lat,lon!!
        
        # compute the hull of each cluster
        hulls[fid+il_fid] = FireVector.cal_hull(coords_geo, 'none')
        
        # if fid in np.arange(100,np.max(all_labels),100):
        #     print(fid, 'processed in', time.time()-t0)
        #     t0 = time.time()
    print('True processing time:', (time.time()-t0)/60)
    
    # write out as geopandas and to file
    gdf = gpd.GeoDataFrame({'fid': list(hulls.keys()), 'geometry': list(hulls.values())}, 
                           crs="EPSG:4326")
    gdf = gdf.to_crs('EPSG:3571')
    
    return gdf

def hull_from_pixels1(arr, geoTrans):
    '''this version uses sklearn measure for labeling (less features than scipy ndimage.label)
    and scipy sparse to compute a sparse matrix in which to search for indices
    *** not in use ***
    '''
    import FireVector, FireIO
    import geopandas as gpd
    from scipy import ndimage, sparse
    from skimage import morphology, measure#, segmentation
    
    # remove small lakes pixels (< 36 30m pixels)
    arr = morphology.remove_small_objects(arr, 5)
    # fill holes (lake islands) --> helps with segementation and boundary extraction
    arr = ndimage.binary_fill_holes(arr)
    
    # cluster the pixels: all neighbouring pixel receive one id
    all_labels = measure.label(arr,background=0) # or ndimage.label
    # find the boundaries of all lakes to speed up the calculation of hulls!
    # edges = segmentation.find_boundaries(arr)
    # inner_edges = np.logical_and(arr, edges)

    # all_labels = np.ma.masked_where(~inner_edges, all_labels)
    # all_labels = np.ma.filled(all_labels, fill_value = 0)
    
    # sparse.find is much faster than np.where!!
    mask_sparse = sparse.csc_matrix(all_labels)
    
    # loop through clusters and extract hulls
    hulls = {}
    # extract lat/lon coordinates from projected image using geotrans
    t0 = time.time()
    print(all_labels.max(), 'lakes to be processed. Estimated time:', int(all_labels.max()/1000/30)+1, 'min')
    for fid in range(1,all_labels.max()+1):
        coords = []
        # t0 = time.time()
        index = sparse.find(mask_sparse==fid)
        # index = np.where(all_labels_masked==fid) # much slower!!
        # print(time.time()-t0)
        for xpx, ypx in zip(index[0]+0.5, index[1]+0.5): # add 0.5 for center of pixel
            x, y = FireIO.pixel2World(geoTrans, ypx, xpx)
            coords.append((x,y))
        x, y = zip(*coords)
        coords_geo = FireIO.polar_to_geo(x,y)
        coords_geo = list(zip(coords_geo[1], coords_geo[0])) # careful! this is lat,lon!!
        
        # compute the hull of each cluster
        #locs = [coords_geo[i] for i,x in enumerate(cid) if x == fid]
        hulls[fid] = FireVector.cal_hull(coords_geo, 'none')
        # if fid in np.arange(100,np.max(all_labels),100):
        #     print(fid, 'processed in', time.time()-t0)
        #     t0 = time.time()
    print('True processing time:', (time.time()-t0)/60)
    
    # write out as geopandas and to file
    gdf = gpd.GeoDataFrame({'fid': list(hulls.keys()), 'geometry': list(hulls.values())}, 
                           crs="EPSG:4326")
    gdf = gdf.to_crs('EPSG:3571')
    
    return gdf

def save_gdf(gdf,year,tilex,tiley):
    ''' Save geopandas to a gpgk file
    year, tilex and tiley are used for the filename and directory
    
    Parameters
    ----------
    gdf : geopandas geodataframe
        geodataframe to write to file
    year : int
        year used for lake data
    tilex: int,
        water tile number in x direction
    tiley: int,
        water tile number in y direction
        
    '''
    from FireConsts import lakedir
    
    path_out = lakedir+str(year)+'/'
    fnm = path_out+'GSW_300m_'+tilex+'_'+tiley+'.gpkg'
    # create a new id column used for indexing 
    #(id column will be dropped when reading gpkg)
    gdf["id"] = gdf.index 
    gdf = gdf.set_index('id')

    # save file
    gdf.to_file(fnm, driver='GPKG')


### OLD: functions for polygonising the raster pixels
# def createProj(projWkt):
#     '''creates an osr spatial reference from wkt'''
#     src_srs = osr.SpatialReference()                                    # create the spatial reference
#     src_srs.ImportFromWkt(projWkt)
#     return src_srs

# def polygonise_tif(fileout, proj, outband, outmask):
#     # create output vector datasource in memory
#     driver = ogr.GetDriverByName('GPKG')
#     temp_ds = driver.CreateDataSource(fileout)

#     # create temp dataset in memory
#     temp_layer = temp_ds.CreateLayer('temp', srs = createProj(proj), geom_type = ogr.wkbPolygon)

#     # store pixel value in first field (right now in int, could be changed to dbl) 
#     gdal.Polygonize(outband, outmask, temp_layer, 0, [], callback=None )


#%% check if all geometries of end product are valid
import geopandas as gpd
from FireConsts import lakedir

year = 2012
path_out = lakedir+str(year)+'/'
files_year = glob.glob(path_out+'*.gpkg')

for file in files_year:
    gdf = gpd.read_file(file)
    if all(gdf.is_valid):
        continue
    else:
        print(file)
        break

#%%

# the minimum output resolution of GSW is 30m,
# however since it is in 4326 the pixels are smaller in the north
# longitudinal resolution can easily be 10m
# we regrid to square 90m pixels first (North Pole LAEA)

# specs
year = 2018
#coarseness = 5 # has to be multiple of 2 or 5
path_in = 'D:/fire_atlas/Data/GlobalSurfaceWater/'+str(year)+'/'
tempfile = 'D:/fire_atlas/Data/GlobalSurfaceWater/vector/temp18.tif'

files_year = glob.glob(path_in+'*0.tif')
# files_year = files_year[:42]

# loop over files in folder
for file in files_year:
    #t0 = time.time()
    #file = files_year[69]
    
    # set output filename
    tilex,tiley =file.split('-')[1:3]
    tilex = tilex[4:6]
    tiley = tiley[3:6]
    
    # check if file should be processed (>50N, not in ocean)
    cond1 = int(tilex) > 8 # for now we skip everything <50N
    cond2 = 52 <= int(tiley) <= 64 # skip Atlantic
    if cond1 or cond2:
        continue
    
    print('processing tile: '+tilex+tiley)
    # open input tif file
    id_ds = gdal.Open(file)
    
    # warp to LAEA, 90m resolution
    if os.path.exists(tempfile):
        os.remove(tempfile)
    warp = gdal.Warp(tempfile,id_ds,dstSRS='EPSG:3571',xRes=90,yRes=90)#, options="-overwrite")
    src1 = warp.ReadAsArray()
    
    # compute binary image (value 3 = permanent water)
    arr = src1 > 2
     
    geoTrans = warp.GetGeoTransform()
    # proj = warp.GetProjection()
    
    # now we extract the pixel locations using the geotrans
    gdf = hull_from_pixels(arr, geoTrans)
    if len(gdf)>0:
        save_gdf(gdf,year,tilex,tiley)
    warp = id_ds = None
    gc.collect()
    
    
    ### OLD: regrid to coarser resolution using coarseness (no warp)
    # src1 = id_ds.ReadAsArray()
    # temp = src1.reshape((src1.shape[0] // coarseness, coarseness,
    #                      src1.shape[1] // coarseness, coarseness))
    # src1 = np.median(temp, axis=(1,3)) # median means there is a tendency to underestimate water?
    
    ### OLD: polygonise the actual pixels!
    
    # correct resolution of output file geotransform
    # [cols, rows] = src1.shape
    # geoTrans = list(id_ds.GetGeoTransform()) # get geotransform to convert pixel to world
    # geoTrans[1] = geoTrans[1]*coarseness
    # geoTrans[5] = geoTrans[5]*coarseness
    # geoTrans = tuple(geoTrans)
    # proj = id_ds.GetProjection() 
    
    # # output the labelled tif as well for comparison
    # drv = gdal.GetDriverByName('GTiff')
    # outds = drv.Create('D:/fire_atlas/Data/GlobalSurfaceWater/vector/labelled.tif', rows, cols, 1, gdal.GDT_Int16)
    # outds.SetGeoTransform(geoTrans)   # sets geotransform from input
    # outds.SetProjection(proj)       # sets projection from input
    # outband = outds.GetRasterBand(1)
    # outband.WriteArray(all_labels)
    # outds = outband = None
    
    # # create output tif in memory (don't put in function - much slower!)
    # drv = gdal.GetDriverByName('MEM')
    # outds = drv.Create('', rows, cols, 1, gdal.GDT_Byte)
    # outds.SetGeoTransform(geoTrans)   # sets geotransform from input
    # outds.SetProjection(proj)       # sets projection from input
    
    # # write np.array to rasterband
    # outband = outds.GetRasterBand(1)
    # outband.WriteArray(src1)
    # outband.SetNoDataValue(0)
    # outmask = outband.GetMaskBand()
    
    # polygonise_tif(fileout, proj, outband, outmask)
    
    # # close files
    # warp = outband = temp_layer = temp_ds = None
    
    #print(time.time() - t0)

#%% OLD: after polygonise filter by area
#from pyproj import Geod
# import os
# import geopandas as gpd
# import pyproj
# import shapely.ops as ops
# from functools import partial

# def geom_area(geom):
#     project = partial(
#         pyproj.transform,
#         pyproj.Proj('EPSG:4326'),
#         pyproj.Proj(
#             proj='aea',
#             lat_1=geom.bounds[1],
#             lat_2=geom.bounds[3]),
#         always_xy=True) # destination coordinate system

#     geom_aea = ops.transform(project, geom)
#     area = geom_aea.area
#     return area

# # specify a named ellipsoid
# # geod = Geod(ellps="WGS84")

# year = 2020
# path_in = 'D:/fire_atlas/Data/GlobalSurfaceWater/vector/'+str(year)+'/fixed/'
# path_out = 'D:/fire_atlas/Data/GlobalSurfaceWater/vector/'+str(year)+'/area/'

# # loop over files in folder
# for filename in os.listdir(path_in)[40:45]:
    
#     if filename == 'GSW_300m_00_128.gpkg':
#         continue
#     file = path_in + "/" + filename
#     outname = path_out + '/' + filename
#     t0 = time.time()

#     # read file
#     gdf = gpd.read_file(file)
#     if len(gdf)==0:
#         continue
#     gdf = gdf.reset_index()
    
#     #gdf2 = gdf.iloc[0:10000]
#     # compute all areas
#     # takes about 30 secs per 1000 features (15 min for 28000)
#     t0 = time.time()
#     areas = gdf['geometry'].apply(geom_area)
#     print('Computed area for '+str(len(gdf))+' geometries in '+str((time.time()-t0)/60)+' min')
    
#     gdf['area'] = areas
#     #gdf = gdf.assign(area = lambda x: geom_area(x['geometry']))
#     #gdf = gdf.assign(area = lambda x: geod.geometry_area_perimeter(x.geometry)[0]) # in geodesic area (only simple polygons)
    
    
#     # filter out small polygons
#     gdf = gdf.loc[gdf['area'] > 140625]
    
#     # write out
#     gdf.to_file(outname, driver='GPKG')




#%% OLD rasterio version of polygonise
# import rasterio
# from rasterio.features import shapes
# from shapely.geometry import shape
# import geopandas as gp

# file = path_in+'yearlyClassification2020-0000040000-0001320000.tif'
# mask = None
# with rasterio.Env():
#     with rasterio.open(file) as src:
#         image = src.read(1) # first band
#         results = ({'properties': {'raster_val': v}, 'geometry': s}
#                    for i, (s, v) in enumerate(shapes(image, mask=mask, transform=src.transform)))

# geoms = list(results)
# gpd_polygonized_raster  = gp.GeoDataFrame.from_features(geoms)





### Other stuff
# pix_size = 40000 # x and y dimension of a single tile
# tiles_mult = 4 # multiplicor for tile number (tiles are numbered in steps of 4)

# NA_out = np.zeros((12*pix_size, 3*pix_size), dtype = np.byte)
# RU_out = np.zeros((20*pix_size, 3*pix_size), dtype = np.byte)

# files = glob.glob(path_in + '*')
# for file in files:
    
#     tiley = int(str.split(file, '-')[1])
#     tilex = int(str.split(file, '-')[2][:-4]) #.tif must be removed
    
#     cond1 = tiley > 8 # for now we skip everything <50N
#     cond2 = 13*tiles_mult <= tilex <=16*tiles_mult # skip Atlantic
#     if cond1 or cond2:
#         continue
    
#     # open file
#     src1 = gdal.Open(file).ReadAsArray()
#     # create a binary image (value 3 = permanent water)
#     src1 = (src1 == 3)*1
    
#     condNA = 1*tiles_mult <= tilex <=12*tiles_mult
#     if condNA:
        

# -*- coding: utf-8 -*-
"""
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
import gdal, ogr, osr
import time
import matplotlib.pyplot as plt


# functions

def index_from_slice(all_labels, start_end, fid):
    mask = all_labels[start_end] == fid
    ind_x = np.arange(start_end[1].start,start_end[1].stop)
    ind_y = np.arange(start_end[0].start,start_end[0].stop)
    grid = np.meshgrid(ind_x, ind_y)
    ind_y = grid[0][mask]
    ind_x = grid[1][mask]
    # for i in range(len(ind_x)):
    #     print(arr[ind_x[i], ind_y[i]])
    return [ind_x, ind_y]

def hull_from_pixels(arr, geoTrans):
    '''This version uses ndimage.label for seqmentation (more features than measure.label)
    and ndimage.find_objects to extract indices (quite fast)
    tested: only using outer boundaries of pixels of each lake, this does not
    work well with calc_hull (creates multipolygons and drops several pixels in calculation)
    '''
    import FireVector, FireIO
    import geopandas as gpd
    from scipy import ndimage#, sparse
    from skimage import morphology, measure, segmentation
    
    # remove small lakes pixels (< 36 30m pixels)
    arr = morphology.remove_small_objects(arr, 5)
    # fill holes (lake islands) --> helps with segementation and boundary extraction
    arr = ndimage.binary_fill_holes(arr)
    
    # cluster the pixels: all neighbouring pixel receive one id
    # all_labels = measure.label(arr,background=0) # or ndimage.label
    all_labels, maxlab = ndimage.label(arr)
    
    # find the boundaries of all lakes to speed up the calculation of hulls!
    inner_edges = segmentation.find_boundaries(arr, mode='inner',background=0)
    inner_edges = segmentation.find_boundaries(inner_edges, mode='thick',background=0)
    inner_edges2 = segmentation.find_boundaries(inner_edges, mode='thick',background=0)
    inner_edges = np.logical_or(inner_edges, inner_edges2)
    inner_edges = np.logical_and(arr, inner_edges)

    all_labels = np.ma.masked_where(~inner_edges, all_labels)
    all_labels = np.ma.filled(all_labels, fill_value = 0)
    
    # use scipy find_objects!
    start_end_idx = ndimage.find_objects(all_labels, max_label=maxlab)
    
    # loop through clusters and extract hulls
    hulls = {}
    # extract lat/lon coordinates from projected image using geotrans
    t0 = time.time()
    print(all_labels.max(), 'lakes to be processed. Estimated time:', int(all_labels.max()/1000*110/60)+1, 'min')
    for fid in range(1,all_labels.max()+1):
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

def hull_from_pixels1(arr, geoTrans):
    '''this version uses sklearn measure for labelling (less features than scipy ndimage.label)
    and scipy sparse to compute a sparse matrix in which to search for indices'''
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
    print(all_labels.max(), 'lakes to be processed. Estimated time:', int(all_labels.max()/1000*110/60)+1, 'min')
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

def save_gdf(gdf, fnm):
    ''' Save geopandas to a gpgk file
    '''
    
    # create a new id column used for indexing 
    #(id column will be dropped when reading gpkg)
    gdf["id"] = gdf.index 
    gdf = gdf.set_index('id')

    # save file
    gdf.to_file(fnm, driver='GPKG')


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

#%%

# the minimum output resolution of GSW is 30m,
# however since it is in 4326 the pixels are smaller in the north
# longitudinal resolution can easily be 10m
# we regrid to tenfold resolution and filter by area afterwards
# (regridding to square pixels would take way too long!)

# specs
year = 2020
#coarseness = 5 # has to be multiple of 2 or 5
path_in = 'D:/fire_atlas/Data/GlobalSurfaceWater/'+str(year)+'/'
# path_out = 'D:/fire_atlas/Data/GlobalSurfaceWater/vector/'+str(year)+'/raw'+'_'+str(coarseness)+'/'
path_out = 'D:/fire_atlas/Data/GlobalSurfaceWater/vector/'+str(year)+'/raw/'
tempfile = 'D:/fire_atlas/Data/GlobalSurfaceWater/vector/temp.tif'

files_year = glob.glob(path_in+'*0.tif')
# files_year = files_year[10:]

# loop over files in folder
for file in files_year:
    #t0 = time.time()
    #file = files_year[69]
    
    # set output filename
    tilex,tiley =file.split('-')[1:3]
    tilex = tilex[4:6]
    tiley = tiley[3:6]
    fileout = path_out+'GSW_300m_'+tilex+'_'+tiley+'.gpkg'
    
    # check if file should be processed (>50N, not in ocean)
    cond1 = int(tilex) > 8 # for now we skip everything <50N
    cond2 = 52 <= int(tiley) <=64 # skip Atlantic
    if cond1 or cond2:
        continue
    
    print('processing tile: '+tilex+tiley)
    # open input tif file
    id_ds = gdal.Open(file)
    # src1 = id_ds.ReadAsArray()
    
    # # regrid to coarser resolution (300m)
    # temp = src1.reshape((src1.shape[0] // coarseness, coarseness,
    #                      src1.shape[1] // coarseness, coarseness))
    # src1 = np.median(temp, axis=(1,3)) # median means there is a tendency to underestimate water?
    
    # warp to LAEA, 1/2viirs resolution
    if os.path.exists(tempfile):
        os.remove(tempfile)
    warp = gdal.Warp(tempfile,id_ds,dstSRS='EPSG:3571',xRes = 90, yRes = 90)#, options="-overwrite")
    src1 = warp.ReadAsArray()
    
    # compute binary image (value 3 = permanent water)
    arr = src1 > 2 # since we also include values of 2.5 we counter the underestimations?
    
    # remove single pixels
    #src1 = morphology.remove_small_objects(src1, 1+np.floor(10/coarseness))
    # src1 = morphology.remove_small_objects(src1, 36)
    
    # fill holes (will make polygonise faster and less error-prone in case of raster.polygonise)
    # src1 = ndimage.binary_fill_holes(src1)
    
    # correct resolution of output file geotransform
    # [cols, rows] = src1.shape
    # geoTrans = list(id_ds.GetGeoTransform()) # get geotransform to convert pixel to world
    # geoTrans[1] = geoTrans[1]*coarseness
    # geoTrans[5] = geoTrans[5]*coarseness
    # geoTrans = tuple(geoTrans)
    # proj = id_ds.GetProjection()  
    geoTrans = warp.GetGeoTransform()
    # proj = warp.GetProjection()
    
    # now we extract the pixel locations using the geotrans
    gdf = hull_from_pixels(arr, geoTrans)
    save_gdf(gdf, fileout)
    warp = id_ds = None
    
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

#%% after correction filter by area
#from pyproj import Geod
import os
import geopandas as gpd
import pyproj
import shapely.ops as ops
from functools import partial

def geom_area(geom):
    project = partial(
        pyproj.transform,
        pyproj.Proj('EPSG:4326'),
        pyproj.Proj(
            proj='aea',
            lat_1=geom.bounds[1],
            lat_2=geom.bounds[3]),
        always_xy=True) # destination coordinate system

    geom_aea = ops.transform(project, geom)
    area = geom_aea.area
    return area

# specify a named ellipsoid
# geod = Geod(ellps="WGS84")

year = 2020
path_in = 'D:/fire_atlas/Data/GlobalSurfaceWater/vector/'+str(year)+'/fixed/'
path_out = 'D:/fire_atlas/Data/GlobalSurfaceWater/vector/'+str(year)+'/area/'

# loop over files in folder
for filename in os.listdir(path_in)[40:45]:
    
    if filename == 'GSW_300m_00_128.gpkg':
        continue
    file = path_in + "/" + filename
    outname = path_out + '/' + filename
    t0 = time.time()

    # read file
    gdf = gpd.read_file(file)
    if len(gdf)==0:
        continue
    gdf = gdf.reset_index()
    
    #gdf2 = gdf.iloc[0:10000]
    # compute all areas
    # takes about 30 secs per 1000 features (15 min for 28000)
    t0 = time.time()
    areas = gdf['geometry'].apply(geom_area)
    print('Computed area for '+str(len(gdf))+' geometries in '+str((time.time()-t0)/60)+' min')
    
    gdf['area'] = areas
    #gdf = gdf.assign(area = lambda x: geom_area(x['geometry']))
    #gdf = gdf.assign(area = lambda x: geod.geometry_area_perimeter(x.geometry)[0]) # in geodesic area (only simple polygons)
    
    
    # filter out small polygons
    gdf = gdf.loc[gdf['area'] > 140625]
    
    # write out
    gdf.to_file(outname, driver='GPKG')




#%% rasterio version
import rasterio
from rasterio.features import shapes
from shapely.geometry import shape
import geopandas as gp

file = path_in+'yearlyClassification2020-0000040000-0001320000.tif'
mask = None
with rasterio.Env():
    with rasterio.open(file) as src:
        image = src.read(1) # first band
        results = ({'properties': {'raster_val': v}, 'geometry': s}
                   for i, (s, v) in enumerate(shapes(image, mask=mask, transform=src.transform)))

geoms = list(results)
gpd_polygonized_raster  = gp.GeoDataFrame.from_features(geoms)





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
        

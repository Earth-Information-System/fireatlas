# -*- coding: utf-8 -*-
"""
Created on Tue Jul 12 20:42:12 2022

@author: rebec
"""

import os
if 'GDAL_DATA' not in os.environ:
    os.environ['GDAL_DATA'] = r'C:/Users/rebec/anaconda3/envs/py3work/Library/share/gdal' 
    os.environ['PROJ_LIB'] = r'C:/Users/rebec/anaconda3/envs/fireatlas/Library/share/proj' 
import numpy as np
import glob
import math
import copy as cp
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import gdal,osr
import affine
from astropy.convolution import convolve
from skimage import segmentation

def build_rtree(coords):
    '''Builds Rtree from a shapely shape
    and optionally uses list of fids as identifier'''
    import rtree
    
    idx = rtree.index.Index() # create new index
    for ind, coord in enumerate(coords):
        # left, bottom, right, top
        idx.insert(ind, (coord[0],coord[1],coord[0],coord[1]), ind)
    
    return idx

wd = 'D:/fire_atlas/2_pipeline/'
res = 200000    # pixel resolution in m
fltr_ma = 2     # moving average filter size in pixels (RADIUS!)
min_fno = 4     # minimum number of fires per pixel (applied after ma filter)

#%% parameters extracted from final fire perimeters

# read all large final fires (preprocessed in R)
df = gpd.read_file(wd+'summary/final_large_3571.gpkg')
df["x"] = df.centroid.x
df["y"] = df.centroid.y
df = pd.DataFrame(df.drop(columns='geometry'))

### assign a 10km grid cell to all fires
# extract minimum and maximum x and y coordinates for new grid
minx = math.floor(min(df.x))
maxx = math.ceil(max(df.x))
miny = math.floor(min(df.y))
maxy = math.ceil(max(df.y))

# we make the raster symmetric to simplify grid creation
xsize = math.ceil(max([abs(minx), maxx])/res)*2
ysize = math.ceil(max([abs(miny), maxy])/res)*2
xabs = xsize*res/2 # this is the corner x value (absolute value)
yabs = ysize*res/2 # this is the corner y value (absolute value)

# group in intervals for aggregation
df['y_int'] = pd.cut(df.y, bins=np.arange(-yabs,yabs+res,res))
df['x_int'] = pd.cut(df.x, bins=np.arange(-xabs,xabs+res,res))
df['x'] = df['x_int'].apply(lambda x: x.mid)
df['y'] = df['y_int'].apply(lambda x: x.mid)


# compute doy from start date
df['year'] = df.tst_year
df['month'] = df.tst_month
df['day'] = df.tst_day
df['tst_doy'] = pd.to_datetime(df[['year', 'month', 'day']]).dt.dayofyear


## aggregate by grid cell
# fire number, size, duration, etc.
df_fno = df.groupby(['x','y']).size()
df_agg = df.groupby(['x','y']).agg({'farea_3571':['sum','mean'],'duration':'mean',
                                    'ted_doy':'mean','tst_doy':'mean'})
df_agg.columns = df_agg.columns.map('_'.join).str.strip('_')

# for land cover take maximum area
df_lcc = df.groupby(['x','y','lcc_coarse'],as_index=False)['farea_3571'].agg(sum).dropna()
df_lcc = df_lcc.sort_values('farea_3571').groupby(['x','y']).tail(1).drop(columns=['farea_3571']).set_index(['x','y'])
lcc_dict = {'Snow/Ice':254,'Water':253,'Barren':252,'Wetlands':40,'Urban':30,
            'Tundra':10,'Sparse boreal forest':11,'Forest boreal':12,'Croplands boreal':13,
            'Forest temperate':22,'Croplands temperate':23,'Grasslands temperate':24,'Temperate mosaic':25}
df_lcc['lcc_code'] = df_lcc.apply(lambda row: lcc_dict[row['lcc_coarse']],axis=1)

# combine all
df_agg = df_agg.join(df_fno.rename('fno')).join(df_lcc)
df_agg = df_agg.reset_index()
df_agg['fno_rel'] = df_agg.fno/((res/1000)**2)/10        # in fires/km2/yr
df_agg['ba_rel'] = (df_agg.farea_3571_sum/((res/1000)**2))/10   # total annual burned area (%)

## assign values to raster

# filter out empty grid cells (these will stay nan)
df_agg_fltr = df_agg[df_agg['fno']>0]
# df_agg_fltr.to_csv(wd+'summary/final_large_agg.csv')

# create array
arr_all = np.full((ysize,xsize,8),fill_value = np.nan)
xint = np.arange(-xabs+res/2, xabs, res)    # x coordinate of pixel centers for array
yint = np.arange(yabs-res/2, -yabs, -res)   # y coordinate of pixel centers for array

# not elegant, but we can loop through the columns and add to values to raster manually
for idx in df_agg_fltr.index:
    x,y = (df_agg_fltr.loc[idx].x, df_agg_fltr.loc[idx].y)
    idx_y = np.where(xint == x)[0]
    idx_x = np.where(yint == y)[0]
    arr_all[idx_x,idx_y,0] = df_agg_fltr.loc[idx].fno
    arr_all[idx_x,idx_y,1] = df_agg_fltr.loc[idx].fno_rel
    arr_all[idx_x,idx_y,2] = df_agg_fltr.loc[idx].ba_rel
    arr_all[idx_x,idx_y,3] = df_agg_fltr.loc[idx].farea_3571_mean
    arr_all[idx_x,idx_y,4] = df_agg_fltr.loc[idx].duration_mean
    arr_all[idx_x,idx_y,5] = df_agg_fltr.loc[idx].tst_doy_mean
    arr_all[idx_x,idx_y,6] = df_agg_fltr.loc[idx].ted_doy_mean
    arr_all[idx_x,idx_y,7] = df_agg_fltr.loc[idx].lcc_code


# visualise to check if grid size is useful
plt.imshow(arr_all[:,:,0])
plt.colorbar()
# fig, ax = plt.subplots(nrows=7, ncols=1)
# for idx,row in enumerate(ax):
#     row.imshow(arr_all[:,:,idx])
# plt.show()

# optional: run moving average
if fltr_ma > 0:
    def circular_filter(radius):
        y, x = np.ogrid[-radius:radius+1, -radius:radius+1]
        mask = x**2 + y**2 <= radius**2           
        return mask*1
    
    kernel = circular_filter(fltr_ma)
    
    # arr_ma = cp.deepcopy(arr_all)
    for band in range(arr_all.shape[2]):
        if band == 7:
            pass
        else:
            # astropy convolve replaces nans by interpolating
            arr_all[:,:,band] = convolve(arr_all[:,:,band],kernel,preserve_nan=True)
    plt.figure()
    plt.imshow(arr_all[:,:,0])
    plt.colorbar()

# optional: remove pixels with low fire number
if min_fno > 1:
    # remove cells with less fires than set in min_fno
    mask2d = arr_all[:,:,0] >= min_fno
    mask = np.broadcast_to(np.expand_dims(mask2d,2), arr_all.shape)
    arr_all[~mask] = np.nan
    plt.figure()
    plt.imshow(arr_all[:,:,0])

# save to geotif
fname = wd+'summary/summary_'+str(int(res/1000))+'km_large_ma'+str(fltr_ma)+'_minfno'+str(min_fno)+'.tif'
nrows,ncols = np.shape(arr_all[:,:,0])
nbands = arr_all.shape[2]
geotransform=(-xabs,res,0,yabs,0, -res)   
output_raster = gdal.GetDriverByName('GTiff').Create(fname,ncols,nrows,nbands,gdal.GDT_Float32)  # Open the file
output_raster.SetGeoTransform(geotransform)
srs = osr.SpatialReference()
srs.ImportFromEPSG(3571)
output_raster.SetProjection(srs.ExportToWkt())
for band in range(nbands):
    output_raster.GetRasterBand(band+1).WriteArray(arr_all[:,:,band])
output_raster = None

## also output the df with all fires and fire types
# fname_df = wd + 'summary/final_ftype.csv'
# df.to_csv(fname_df)

#%% rasterise tree line and distance from tree line
import rasterio
infile = 'D:/shared_data/ecoregions/cav/Raster CAVM/raster_cavm_v1.tif'
cavm_bin_file = 'D:/shared_data/ecoregions/cav/Raster CAVM/raster_binary.tif'
outfile_temp = wd + 'summary/treeline' + str(int(res/1000/4)) + '.tif'
outfile_temp1 = wd + 'summary/treeline' + str(int(res/1000)) + '.tif'
outfile_temp2 = wd + 'summary/treeline_dist' + str(int(res/1000)) + '.tif'
outname = wd + 'summary/treeline_dist.tif'

# create a binary image of CAVM
with rasterio.open(infile, 'r') as src:
    arr = src.read()
    arr = np.squeeze(arr) > 90
    arr = (arr*1).astype(np.byte)
    with rasterio.open(cavm_bin_file,'w',height=src.height,width=src.width,driver='Gtiff',
                       count=1,crs=src.crs,transform=src.transform,dtype=arr.dtype) as dst_out:
        dst_out.write(arr, 1)

# regrid cavm to output grid
warp = gdal.Warp(outfile_temp,cavm_bin_file,dstSRS='EPSG:3571',
                 xRes=res/4, yRes=res/4, outputBounds=(-xabs,-yabs,xabs,yabs))#, resampleAlg = 'lanczos')
# pretty much all resampling algorithms underestimate the treeline...
# maybe there is a better way?

# aggregate to (coarser) output resolution
arr = warp.ReadAsArray()
warp = None
arr_coarse = arr.reshape(int(arr.shape[0]/4), 4, int(arr.shape[1]/4), 4).sum(axis=(1, 3))
arr_coarse = (arr_coarse < 16).astype(np.uint8)
with rasterio.open(outfile_temp, 'r') as src:
    trans = src.transform
    trans = affine.Affine(trans.a*4, trans.b, trans.c, trans.d, trans.e*4, trans.f)
    with rasterio.open(outfile_temp1,'w',height=src.height/4,width=src.width/4,driver='Gtiff',
                       count=1,crs=src.crs,transform=trans,dtype=arr.dtype) as dst_out:
        dst_out.write(arr_coarse, 1)

# compute distance from treeline for each pixel
tl_raster = cp.deepcopy(arr_coarse).astype(np.int)
tl_raster[(tl_raster == 0)] = -99
tl_raster[(tl_raster == 1)] = 0

value = 0
newsum = 1
while newsum > 0:
    value += 1
    edges = segmentation.find_boundaries((tl_raster > -99)*1, mode = 'outer')
    tl_raster[edges] = value
    newsum = np.sum(edges)
    
# compute centroid distance
# insert all tundra pixels into rtree index
ds = rasterio.open(outfile_temp1, 'r')
x,y = np.where(arr_coarse == 1)
coords_tundra = ds.xy(x, y)
idx = build_rtree(coords_tundra)

# loop through boreal pixel, find nearest 
x,y = np.where(arr_coarse == 0)
xx,yy = ds.xy(x, y)
coords_taiga = list(zip(xx,yy))

tl_raster[tl_raster == -99] = 0

result = []
for ind,coord in enumerate(coords_taiga):
    # find nearest neighbour
    neighbours = list(idx.nearest(coord, 1, objects='raw'))
    if len(neighbours) > 1:
        print('more than one neighbour')
    tl_coord = coords_tundra[neighbours[0]]
    # compute centroid distance
    dist = math.sqrt((coord[0]-tl_coord[0])**2 + (coord[1]-tl_coord[1])**2)
    # save to list
    result.append((coord[0], coord[1], dist))
    # also write to raster
    tl_raster[x[ind],y[ind]] = dist
   
with rasterio.open(outfile_temp1, 'r') as src:
    with rasterio.open(outfile_temp2,'w',height=src.height,width=src.width,driver='Gtiff',
                       count=1,crs=src.crs,transform=src.transform,dtype=tl_raster.dtype) as dst_out:
        dst_out.write(tl_raster, 1)
df = pd.DataFrame(result, columns=['x','y','dist_tl'])
df.to_csv(wd + 'summary/distance_treeline_centroids.csv')


#%% number of fires contributing to 50% of burned area in grid cell


#%% spread parameters (extracted from daily perimeters)

files = glob.glob(wd + '20??/Summary/fsummary_sf.csv', recursive=True)
df = pd.concat(map(pd.read_csv, files), ignore_index=True)


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
# import copy as cp
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import gdal,osr
from haversine import haversine

wd = 'D:/fire_atlas/2_pipeline/'
res = 5

# create a 5 degree grid
lats = np.arange(90- res / 2, 40, -res)
lons = np.arange(-180 + res / 2, 180, res)
# arr_blank = np.zeros((len(lons),len(lats)))
# lons_mesh,lats_mesh = np.meshgrid(lons, lats)
# plt.imshow(lats_mesh)

# arr_fno = cp.deepcopy(arr_blank)    # number of fires
# arr_ignno = cp.deepcopy(arr_blank)  # number of ignitions
# arr_fsize = cp.deepcopy(arr_blank)  # average fire size
# arr_ignno = cp.deepcopy(arr_blank)  # number of ignitions
# arr_maxmon = cp.deepcopy(arr_blank) # average start date of fires
# arr_minmon = cp.deepcopy(arr_blank) # average end date of fires

# read all final fires
for year in range(2012,2022):
    gdf_year = gpd.read_file(wd+str(year)+'/Summary/final_viirs' + str(year) + '.gpkg')
    gdf_year["x"] = gdf_year.centroid.x
    gdf_year["y"] = gdf_year.centroid.y
    gdf_year = pd.DataFrame(gdf_year.drop(columns='geometry'))
    
    if year == 2012:
        df = gdf_year
    else:
        df = pd.concat([df, gdf_year])
df = df.reset_index()
df['year'] = df.tst_year
df['month'] = df.tst_month
df['day'] = df.tst_day
df['tst_doy'] = pd.to_datetime(df[['year', 'month', 'day']]).dt.dayofyear

# assign a grid cell to all fires
df['lat_int'] = pd.cut(df.y, bins=np.arange(40,90,res))
df['lon_int'] = pd.cut(df.x, bins=np.arange(-180,180,res))
df['lon'] = df['lon_int'].apply(lambda x: x.mid)
df['lat'] = df['lat_int'].apply(lambda x: x.mid)

# aggregate by grid cell
df_fno = df.groupby(['lat','lon']).size()
df_agg = df.groupby(['lat','lon']).agg({'farea':'mean','duration':'mean',
                                        'lake_area':'mean','lake_no':'mean',
                                        'ted_doy':'mean','tst_doy':'mean'})
df_lcc = df.groupby(['lat','lon'])['lcc_final'].agg(pd.Series.mode)
df_agg = df_agg.join(df_fno.rename('fno')).join(df_lcc.rename('lcc'))
# df_agg['avg_lakesize'] = df_agg.lake_area/df_agg.lake_no
df_agg = df_agg.reset_index()
df_agg['px_size'] = df_agg.apply(lambda row: haversine((row.lat,0),(row.lat,res))*haversine((0,0),(res,0)), axis=1)

# assign values to raster
df_agg_fltr = df_agg[df_agg['fno']>0]

# not elegant, but we can loop through the columns and add to values to raster manually
arr_all = np.full((len(lats),len(lons),8),fill_value = np.nan)
for idx in df_agg_fltr.index:
    lat,lon = (df_agg_fltr.loc[idx].lat, df_agg_fltr.loc[idx].lon)
    idx_y = np.where(lons == lon)[0]
    idx_x = np.where(lats == lat)[0]
    arr_all[idx_x,idx_y,0] = df_agg_fltr.loc[idx].fno/df_agg_fltr.loc[idx].px_size*1000000/12 # in fires/km2/yr
    arr_all[idx_x,idx_y,1] = df_agg_fltr.loc[idx].farea
    arr_all[idx_x,idx_y,2] = df_agg_fltr.loc[idx].duration
    arr_all[idx_x,idx_y,3] = df_agg_fltr.loc[idx].tst_doy
    arr_all[idx_x,idx_y,4] = df_agg_fltr.loc[idx].ted_doy
    arr_all[idx_x,idx_y,5] = df_agg_fltr.loc[idx].lake_area
    arr_all[idx_x,idx_y,6] = df_agg_fltr.loc[idx].lake_no
    # arr_all[idx,idy,7] = df_agg_fltr.loc[latlon].lcc

# visualise
# plt.imshow(arr_all[:,:,0])
fig, ax = plt.subplots(nrows=7, ncols=1)
for idx,row in enumerate(ax):
    row.imshow(arr_all[:,:,idx])
plt.show()

# save to geotif
nrows,ncols = np.shape(arr_all[:,:,0])
nbands = arr_all.shape[2]
geotransform=(-180,res,0,90,0, -res)   
output_raster = gdal.GetDriverByName('GTiff').Create(wd+'summary_'+str(res)+'deg.tif',ncols,nrows,nbands,gdal.GDT_Float32)  # Open the file
output_raster.SetGeoTransform(geotransform)
srs = osr.SpatialReference()
srs.ImportFromEPSG(4326)
output_raster.SetProjection(srs.ExportToWkt())
for band in range(nbands):
    output_raster.GetRasterBand(band+1).WriteArray(arr_all[:,:,band])
output_raster = None

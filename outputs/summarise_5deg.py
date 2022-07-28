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
# import copy as cp
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import gdal,osr
from haversine import haversine, Unit

wd = 'D:/fire_atlas/2_pipeline/'
res = 5
exclude_small = True

#%% parameters extracted from final fire perimeters
# create a n degree grid (n set by res parameter)
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
    # final fires
    gdf_year = gpd.read_file(wd+str(year)+'/Summary/final_viirs' + str(year) + '.gpkg')
    gdf_year["x"] = gdf_year.centroid.x
    gdf_year["y"] = gdf_year.centroid.y
    gdf_year = pd.DataFrame(gdf_year.drop(columns='geometry'))
    
    # ignitions
    gdf_ign_year = gpd.read_file(wd+str(year)+'/Summary/ign' + str(year) + '.gpkg')
    gdf_ign_year["x"] = gdf_ign_year.centroid.x
    gdf_ign_year["y"] = gdf_ign_year.centroid.y
    gdf_ign_year = pd.DataFrame(gdf_ign_year.drop(columns='geometry'))
    
    if year == 2012:
        df = gdf_year
        df_ign = gdf_ign_year
    else:
        df = pd.concat([df, gdf_year])
        df_ign = pd.concat([df_ign, gdf_ign_year])

# reset the index to make all columns accessible
df = df.reset_index()

# optional: filter out small fires
if exclude_small:
    df = df[df.farea > 1.890771] # based on 95th percentile of size/area distribution

# compute doy from start date
df['year'] = df.tst_year
df['month'] = df.tst_month
df['day'] = df.tst_day
df['tst_doy'] = pd.to_datetime(df[['year', 'month', 'day']]).dt.dayofyear

# assign a 5 degree grid cell to all fires
df['lat_int'] = pd.cut(df.y, bins=np.arange(40,90+res,res))
df['lon_int'] = pd.cut(df.x, bins=np.arange(-180,180+res,res))
df['lon'] = df['lon_int'].apply(lambda x: x.mid)
df['lat'] = df['lat_int'].apply(lambda x: x.mid)

df_ign['lat_int'] = pd.cut(df_ign.y, bins=np.arange(40,90+res,res))
df_ign['lon_int'] = pd.cut(df_ign.x, bins=np.arange(-180,180+res,res))
df_ign['lon'] = df_ign['lon_int'].apply(lambda x: x.mid)
df_ign['lat'] = df_ign['lat_int'].apply(lambda x: x.mid)


# aggregate by grid cell
df_ignno = df_ign.groupby(['lat','lon']).size()
df_fno = df.groupby(['lat','lon']).size()

df_agg = df.groupby(['lat','lon']).agg({'farea':['sum','mean'],'duration':'mean',
                                        'ted_doy':'mean','tst_doy':'mean'})
df_agg.columns = df_agg.columns.map('_'.join).str.strip('_')

# df_lcc = df.groupby(['lat','lon'])['lcc_coarse'].agg(pd.Series.mode)
df_lcc = df.groupby(['lat','lon','lcc_coarse'],as_index=False)['farea'].agg(sum).dropna()
df_lcc = df_lcc.sort_values('farea').groupby(['lat','lon']).tail(1).drop(columns=['farea']).set_index(['lat','lon'])
lcc_dict = {'Snow/Ice':254,'Water':253,'Barren':252,'Wetlands':40,'Urban':30,
            'Tundra':10,'Sparse boreal forest':11,'Forest boreal':12,'Croplands boreal':13,
            'Forest temperate':22,'Croplands temperate':23,'Grasslands temperate':24,'Temperate mosaic':25}
df_lcc['lcc_code'] = df_lcc.apply(lambda row: lcc_dict[row['lcc_coarse']],axis=1)

df_agg = df_agg.join(df_ignno.rename('ignno')).join(df_fno.rename('fno')).join(df_lcc)
df_agg = df_agg.reset_index()
df_agg['px_size'] = df_agg.apply(lambda row: haversine((row.lat,0),(row.lat,res),unit=Unit.KILOMETERS)* haversine((0,0),(res,0),unit=Unit.KILOMETERS), axis=1)
df_agg['fno_rel'] = df_agg.fno/df_agg.px_size/12        # in fires/km2/yr
df_agg['ignno_rel'] = df_agg.ignno/df_agg.px_size/12    # in fires/km2/yr
df_agg['ba_rel'] = df_agg.farea_sum/df_agg.px_size/12   # total annual burned area (%)

# assign values to raster
df_agg_fltr = df_agg[df_agg['fno']>0]

# not elegant, but we can loop through the columns and add to values to raster manually
arr_all = np.full((len(lats),len(lons),8),fill_value = np.nan)
for idx in df_agg_fltr.index:
    lat,lon = (df_agg_fltr.loc[idx].lat, df_agg_fltr.loc[idx].lon)
    idx_y = np.where(lons == lon)[0]
    idx_x = np.where(lats == lat)[0]
    arr_all[idx_x,idx_y,0] = df_agg_fltr.loc[idx].fno_rel
    arr_all[idx_x,idx_y,1] = df_agg_fltr.loc[idx].ignno_rel
    arr_all[idx_x,idx_y,2] = df_agg_fltr.loc[idx].ba_rel
    arr_all[idx_x,idx_y,3] = df_agg_fltr.loc[idx].farea_mean
    arr_all[idx_x,idx_y,4] = df_agg_fltr.loc[idx].duration_mean
    arr_all[idx_x,idx_y,5] = df_agg_fltr.loc[idx].tst_doy_mean
    arr_all[idx_x,idx_y,6] = df_agg_fltr.loc[idx].ted_doy_mean
    # arr_all[idx_x,idx_y,6] = df_agg_fltr.loc[idx].lake_area/df_agg_fltr.loc[idx].px_size # in area/km2
    # arr_all[idx_x,idx_y,7] = df_agg_fltr.loc[idx].lake_no
    arr_all[idx_x,idx_y,7] = df_agg_fltr.loc[idx].lcc_code

# visualise
plt.imshow(arr_all[:,:,7])
# fig, ax = plt.subplots(nrows=7, ncols=1)
# for idx,row in enumerate(ax):
#     row.imshow(arr_all[:,:,idx])
# plt.show()

# save to geotif
if exclude_small:
    fname = wd+'summary/summary_'+str(res)+'deg_large.tif'
else:
    fname = wd+'summary/summary_'+str(res)+'deg.tif'
nrows,ncols = np.shape(arr_all[:,:,0])
nbands = arr_all.shape[2]
geotransform=(-180,res,0,90,0, -res)   
output_raster = gdal.GetDriverByName('GTiff').Create(fname,ncols,nrows,nbands,gdal.GDT_Float32)  # Open the file
output_raster.SetGeoTransform(geotransform)
srs = osr.SpatialReference()
srs.ImportFromEPSG(4326)
output_raster.SetProjection(srs.ExportToWkt())
for band in range(nbands):
    output_raster.GetRasterBand(band+1).WriteArray(arr_all[:,:,band])
output_raster = None


## also output the df with all fires and fire types
fname_df = wd + 'summary/final_ftype.csv'
df.to_csv(fname_df)

#%% number of fires contributing to 50% of burned area in grid cell


#%% spread parameters (extracted from daily perimeters)

files = glob.glob(wd + '20??/Summary/fsummary_sf.csv', recursive=True)
df = pd.concat(map(pd.read_csv, files), ignore_index=True)


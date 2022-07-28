# -*- coding: utf-8 -*-
"""
Created on Thu Jul 14 19:22:37 2022

@author: rebec
"""
import os
if 'GDAL_DATA' not in os.environ:
    os.environ['GDAL_DATA'] = r'C:/Users/rebec/anaconda3/envs/py3work/Library/share/gdal' 
    os.environ['PROJ_LIB'] = r'C:/Users/rebec/anaconda3/envs/fireatlas/Library/share/proj' 
import glob
import rasterio
import rasterio.mask
from rasterio.merge import merge

for year in range(2012,2021):
    pathLCT = 'D:/fire_atlas/Data/modis_lc_dave/' + str(year) + '_LCT/' + str(year) + '/'
    files = glob.glob(pathLCT+'*')
    
    # merge sinusoidal modis land cover tiles to mosaic
    filelist = []
    for file in files:
        src = rasterio.open(file)
        filelist.append(src)
    
    # create and plot a mosiac
    mosaic, out_trans = merge(filelist)
    
    # write out mosiac band
    out_fp = pathLCT + 'merged.tif'
    out_meta = src.meta.copy()
    out_meta.update({"driver": "GTiff",
                     "height": mosaic.shape[1], "width": mosaic.shape[2],
                     "transform": out_trans})
    with rasterio.open(out_fp, "w", **out_meta) as dest:
        dest.write(mosaic)
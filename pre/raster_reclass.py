# -*- coding: utf-8 -*-
"""
Created on Wed Aug 18 14:00:23 2021

@author: rebec
"""


import os
import numpy as np
import copy
import matplotlib.pyplot as plt

# GDAL needs an environment variable GDAL_DATA to find and use projection info
# this has to be added before gdal is imported
if 'GDAL_DATA' not in os.environ:
    # os.environ['GDAL_DATA'] = r'C:/Users/rebec/Anaconda3/Library/share/gdal'              # old anaconda version
    os.environ['GDAL_DATA'] = r'C:/Users/rebec/Anaconda3/envs/py3work/Library/share/gdal'   # new anaconda version (Aug 2020)
from osgeo import gdal#, osr

''' Supporting functions '''
def save2gtif(arr, outfile, cols, rows, geotrans, proj):
    driver = gdal.GetDriverByName('GTiff')
    ds = driver.Create(outfile, cols, rows, 1, gdal.GDT_Byte)
    ds.SetGeoTransform(geotrans)
    ds.SetProjection(proj)
    band = ds.GetRasterBand(1)
    band.WriteArray(arr)

#%% Reclassifying Hugelius Histel and Histosol fractions to binary raster

path = r'D:\fire_atlas\Data\peat\Hugelius_etal_2020_PNAS_grids\Grids_TIFF_WGS84'

# load histel data (permafrost peatlands)
ds_histel = gdal.Open(path+'/Histel_fraction_WGS84.tif')
arr = ds_histel.ReadAsArray()

# load histosol data (non-frozen peatlands)
ds_histosol = gdal.Open(path+'/Histosol_fraction_WGS84.tif')
arr1 = ds_histosol.ReadAsArray()

# combine both
arr = arr + arr1
arr[arr < 0] = 0

# plt.imshow(arr)

# set a threshold and turn into binary image
arr_bin = arr > 10

# write out to file
gt = ds_histel.GetGeoTransform()
proj = ds_histel.GetProjection()
outfile = path+'/peat_greater_10.tif'
save2gtif(arr_bin, outfile, arr_bin.shape[1], arr_bin.shape[0], gt, proj)


#%% Reclassifying CAVM raster version to coarse classes

# original classes
class_orig = [1,2,3,4,5,21,22,23,24,31,32,33,34,41,42,43,91,92,93,99,127]
class_new = [1,2,1,1,1,2,2,2,2,3,3,3,3,4,4,4,0,0,0,0,0,0]
lc_dict = dict(zip(class_orig, class_new))

# load data
path = 'D:/shared_data/ecoregions/cav/Raster CAVM/'
path_out = 'D:/fire_atlas/Data/cavm/'
ds = gdal.Open(path+'raster_cavm_v1.tif')
arr = ds.ReadAsArray()
arr_out = copy.deepcopy(arr)

for lcc in class_orig:
    arr_out[arr_out == lcc] = lc_dict[lcc]

# write out to file
gt = ds.GetGeoTransform()
proj = ds.GetProjection()
outfile = path_out+'/cavm_coarse.tif'
save2gtif(arr_out, outfile, arr_out.shape[1], arr_out.shape[0], gt, proj)


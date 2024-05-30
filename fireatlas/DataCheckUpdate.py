""" DataUpdate
This module include functions used to check and update needed data files
"""
import urllib.request
import os
import fsspec
import xarray as xr
import tempfile
import pandas as pd
import warnings

from datetime import date, timedelta

from fireatlas import settings
from fireatlas.FireLog import logger
from fireatlas.preprocess import preprocess_input_file

# ------------------------------------------------------------------------------
# update external dataset
# ------------------------------------------------------------------------------
def wget(url, **kwargs):
    target_dir = "."
    if "locdir" in kwargs:
        target_dir = kwargs.pop("locdir")
    target_file = os.path.join(target_dir, os.path.basename(url))
    logger.info(f"Downloading {url} to {target_file}")
    opener = urllib.request.build_opener(urllib.request.HTTPCookieProcessor())
    request = urllib.request.Request(url)
    if "header" in kwargs:
        header = kwargs.pop("header")
        assert header == "NASA", f"Non-standard header is not implemented: {header}"
        request.add_header("Authorization", "Bearer ZW9ybGFuZDpaV3hwYW1Gb0xtOXliR0Z1WkVCdVlYTmhMbWR2ZGc9PToxNjQyNzE4ODAyOjQyYzMzM2ViODViOWI3OTVlYzAyYTdmYWE2ZjYwYjFjZTc5MGJmNDg")
    if len(kwargs) > 0:
        logger.info(f"WARNING: Ignoring unused wget arguments: {list(kwargs.keys())}")
    with opener.open(request) as response, fsspec.open(target_file, "wb") as f:
        f.write(response.read())
    return target_file

def update_VNP14IMGTDL():
    ''' Batch read and extract update_VJ114IMGTDL data'''
    # The directory to save VNP14IMGTDL data
    data_dir = os.path.join(settings.dirextdata, "VIIRS", "VNP14IMGTDL/")

    # Derive the date periods needed to download
    today = date.today()
    fnms = settings.fs.glob(os.path.join(data_dir, f'SUOMI_VIIRS_C2_Global_VNP14IMGTDL_NRT_{str(today.year)}*.txt'))
    if len(fnms)==0:
        ndays = 0
    else:
        doys = [int(d[-7:-4]) for d in fnms]
        ndays = max(doys)
    dstart = date(today.year,1,1) + timedelta(days=ndays)
    # to avoid any incomplete data just check the last 10 days
    # this shoudl be very quick and it's okay if we are duplicating efforts
    dstart = dstart - timedelta(days=10)

    # Do the download process
    urldir = "https://nrt4.modaps.eosdis.nasa.gov/api/v2/content/archives/FIRMS/suomi-npp-viirs-c2/Global/"
    for d in pd.date_range(dstart,today):
        urlfnm = urldir + "SUOMI_VIIRS_C2_Global_VNP14IMGTDL_NRT_"+d.strftime('%Y%j')+".txt"
        try:
            downloaded_filepath = wget(url=urlfnm,locdir=data_dir,robots_off=True,no_wget=False,timestamping=True,header='NASA')
        except Exception as e:
            print("\nCould not download VNP14IMGTDL data for",d)
            print('Error message:',e)
            continue
        preprocess_input_file(downloaded_filepath)


def update_VJ114IMGTDL():
    ''' Batch read and extract update_VJ114IMGTDL data'''
    # The directory to save VNP14IMGTDL data
    data_dir = os.path.join(settings.dirextdata, 'VIIRS', 'VJ114IMGTDL/')
    # Derive the date periods needed to download
    today = date.today()
    fnms = settings.fs.glob(os.path.join(data_dir, f'J1_VIIRS_C2_Global_VJ114IMGTDL_NRT_{str(today.year)}*.txt'))
    if len(fnms)==0:
        ndays = 0
    else:
        doys = [int(d[-7:-4]) for d in fnms]
        ndays = max(doys)
    dstart = date(today.year,1,1) + timedelta(days=ndays)
    # to avoid any incomplete data just check the last 10 days
    # this shoudl be very quick and it's okay if we are duplicating efforts
    dstart = dstart - timedelta(days=10)

    # Do the download process
    urldir = "https://nrt4.modaps.eosdis.nasa.gov/api/v2/content/archives/FIRMS/noaa-20-viirs-c2/Global/"
    for d in pd.date_range(dstart,today):
        urlfnm = urldir + "J1_VIIRS_C2_Global_VJ114IMGTDL_NRT_"+d.strftime('%Y%j')+".txt"
        try:
            downloaded_filepath = wget(url=urlfnm,locdir=data_dir,robots_off=True,no_wget=False,timestamping=True,header='NASA')
        except Exception as e: 
            print("\nCould not download VJ114IMGTDL data for",d)
            print('Error message:',e)
            continue
        preprocess_input_file(downloaded_filepath)

def update_GridMET_fm1000():
    ''' Get updated GridMET data (including fm1000)
    '''
    # The directory to save GridMET data
    data_dir = os.path.join(settings.dirextdata, 'GridMET/')

    today = date.today()

    # Do the download process
    urldir = "http://www.northwestknowledge.net/metdata/data/"
    # strvars = ['vpd','pr','tmmn','tmmx','vs','fm100','fm1000','bi','pdsi']
    strvars = ['fm1000']
    for strvar in strvars:
        target_file = strvar + '_' + str(today.year) + '.nc'
        urlfnm = urldir + target_file
        with tempfile.TemporaryDirectory() as tempdir:
            wget(urlfnm, locdir=tempdir)
            file_name = os.path.join(tempdir, target_file)
            # Convert to Zarr
            zarrfile = target_file.replace(".nc", ".zarr")
            print(f"Converting {target_file} to {zarrfile}.")
            dat = xr.open_dataset(file_name)
            dat.to_zarr(os.path.join(data_dir, zarrfile), mode="w")

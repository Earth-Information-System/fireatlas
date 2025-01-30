""" DataUpdate
This module include functions used to check and update needed data files
"""
import os
import fsspec
import xarray as xr
import tempfile
import requests

from datetime import date, date

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

    headers = {}
    if "header" in kwargs:
        header = kwargs.pop("header")
        assert header == "NASA", f"Non-standard header is not implemented: {header}"
        headers["Authorization"] = "Bearer eyJ0eXAiOiJKV1QiLCJvcmlnaW4iOiJFYXJ0aGRhdGEgTG9naW4iLCJzaWciOiJlZGxqd3RwdWJrZXlfb3BzIiwiYWxnIjoiUlMyNTYifQ.eyJ0eXBlIjoiVXNlciIsInVpZCI6InRtY2NhYmUiLCJleHAiOjE3NDE0NDUxMjMsImlhdCI6MTczNjI2MTEyMywiaXNzIjoiaHR0cHM6Ly91cnMuZWFydGhkYXRhLm5hc2EuZ292IiwiaWRlbnRpdHlfcHJvdmlkZXIiOiJlZGxfb3BzIiwiYWNyIjoiZWRsIiwiYXNzdXJhbmNlX2xldmVsIjozfQ.Wm8LrY4wvAty2q-EH00vdMOiB1KKOzXIaWn_WAqpeOTVudVLbNEZN8RF2dJVICtkZB3cvmK4e6O0lsejqhWMyg_HkLcioyf8hKxTJHVRsLLKD10sXbXZFQI0f_BsrH4akMidj_PlYY7xCzfHAguEydnxzkC26dJLLWpls6Gcoioba1AjFbGNdVJTsT5kAwsh0-i-Oz_EqdDO6L7Juteg4xrL6bGfyZJTB8I33NgpQPqiCWOK1QNgSItT0JMiYauCBWnyQlUgaPquPMPA-b9CKWvBiw9X-7IbZ-GDIa4Oy6b_2qjEN2oyBR5_slBDGWjI5CFotc0471lWLAc9Z-LwLQ"

    if len(kwargs) > 0:
        logger.debug(f"WARNING: Ignoring unused wget arguments: {list(kwargs.keys())}")

    response = requests.get(url, headers=headers)
    response.raise_for_status()  # This will raise an HTTPError for bad requests (4XX or 5XX)

    with fsspec.open(target_file, "wb") as f:
        f.write(response.content)
    return target_file


def update_VNP14IMGTDL(d: date):
    ''' Batch read and extract update_S-NPP data'''
    # The directory to save VNP14IMGTDL data
    data_dir = os.path.join(settings.dirextdata, "VIIRS", "VNP14IMGTDL/")

    # Do the download process
    urldir = "https://nrt3.modaps.eosdis.nasa.gov/api/v2/content/archives/FIRMS/suomi-npp-viirs-c2/Global/"
    urlfnm = urldir + "SUOMI_VIIRS_C2_Global_VNP14IMGTDL_NRT_"+d.strftime('%Y%j')+".txt"
    try:
        downloaded_filepath = wget(url=urlfnm,locdir=data_dir,robots_off=True,no_wget=False,timestamping=True,header='NASA')
        preprocess_input_file(downloaded_filepath)
    except Exception as e:
        logger.warning(f"Could not download VNP14IMGTDL data for {d}")
        logger.warning(f"Error message: {str(e)}")


def update_VJ114IMGTDL(d: date):
    ''' Batch read and extract update_NOAA20 data'''
    # The directory to save VJ114IMGTDL data
    data_dir = os.path.join(settings.dirextdata, 'VIIRS', 'VJ114IMGTDL/')

    # Do the download process
    urldir = "https://nrt3.modaps.eosdis.nasa.gov/api/v2/content/archives/FIRMS/noaa-20-viirs-c2/Global/"
    urlfnm = urldir + "J1_VIIRS_C2_Global_VJ114IMGTDL_NRT_"+d.strftime('%Y%j')+".txt"
    try:
        downloaded_filepath = wget(url=urlfnm,locdir=data_dir,robots_off=True,no_wget=False,timestamping=True,header='NASA')
        preprocess_input_file(downloaded_filepath)
    except Exception as e:
        logger.warning(f"Could not download VJ114IMGTDL data for {d}")
        logger.warning(f"Error message: {str(e)}")


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

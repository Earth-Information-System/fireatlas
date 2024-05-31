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
        headers["Authorization"] = "Bearer eyJ0eXAiOiJKV1QiLCJvcmlnaW4iOiJFYXJ0aGRhdGEgTG9naW4iLCJzaWciOiJlZGxqd3RwdWJrZXlfb3BzIiwiYWxnIjoiUlMyNTYifQ.eyJ0eXBlIjoiVXNlciIsInVpZCI6ImdzZmNfbGFuZHNsaWRlcyIsImV4cCI6MTcyMDgxNzI5OSwiaWF0IjoxNzE1NjMzMjk5LCJpc3MiOiJFYXJ0aGRhdGEgTG9naW4ifQ.mgwyy26sW8i6cEopY2QOdt27Mb0e_kbLOO9aEFn_yhqpOqHwuvV7-RSv1w5P0p3xtVM_jV-Lji0dL9plvtNCrx1rlew3Re7orjCF7cwfcBgjHmvA3Qdf9SR-JEsXTeJUmCgXFElAEYXj8jGtIR6r4jGHylfcOveZqUebVMtBlk_nJnavRfJrsg5Jq-t2GDFDaYXvm4tbtPtoXSZBsHIiuaXnxtTX1MXC_KZZv5sLT9vw-iYNGdm_64DVnY10nPLsfv9AoEKHL81EWbTTcs8LDJIqMVaMF7X9sHxPidF0qF0S2gEKkCU5G5Thl0c7ijQzUorZ878L_HbTX3ezjm_TwQ"

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
    urldir = "https://nrt4.modaps.eosdis.nasa.gov/api/v2/content/archives/FIRMS/suomi-npp-viirs-c2/Global/"
    urlfnm = urldir + "SUOMI_VIIRS_C2_Global_VNP14IMGTDL_NRT_"+d.strftime('%Y%j')+".txt"
    try:
        downloaded_filepath = wget(url=urlfnm,locdir=data_dir,robots_off=True,no_wget=False,timestamping=True,header='NASA')
        preprocess_input_file(downloaded_filepath)
    except Exception as e:
        logger.warn("\nCould not download VNP14IMGTDL data for",d)
        logger.warn('Error message:', str(e))


def update_VJ114IMGTDL(d: date):
    ''' Batch read and extract update_NOAA20 data'''
    # The directory to save VJ114IMGTDL data
    data_dir = os.path.join(settings.dirextdata, 'VIIRS', 'VJ114IMGTDL/')

    # Do the download process
    urldir = "https://nrt4.modaps.eosdis.nasa.gov/api/v2/content/archives/FIRMS/noaa-20-viirs-c2/Global/"
    urlfnm = urldir + "J1_VIIRS_C2_Global_VJ114IMGTDL_NRT_"+d.strftime('%Y%j')+".txt"
    try:
        downloaded_filepath = wget(url=urlfnm,locdir=data_dir,robots_off=True,no_wget=False,timestamping=True,header='NASA')
        preprocess_input_file(downloaded_filepath)
    except Exception as e:
        logger.warn("\nCould not download VJ114IMGTDL data for",d)
        logger.warn('Error message:', str(e))


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

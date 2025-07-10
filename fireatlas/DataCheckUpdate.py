""" DataUpdate
This module include functions used to check and update needed data files
"""
import os
import fsspec
import time
import xarray as xr
import tempfile
import requests
import pandas as pd

from datetime import date
from typing import Literal

from fireatlas import settings
from fireatlas.FireLog import logger
from fireatlas.preprocess import preprocess_input_file

MAP_KEY = "0e50658bd8e8ea368db7379b0be28630"
N_MAX_RETRIES = 30

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
        headers["Authorization"] = "Bearer eyJ0eXAiOiJKV1QiLCJvcmlnaW4iOiJFYXJ0aGRhdGEgTG9naW4iLCJzaWciOiJlZGxqd3RwdWJrZXlfb3BzIiwiYWxnIjoiUlMyNTYifQ.eyJ0eXBlIjoiVXNlciIsInVpZCI6InpiZWNrZXIiLCJleHAiOjE3NTQyMjY5MDMsImlhdCI6MTc0OTA0MjkwMywiaXNzIjoiaHR0cHM6Ly91cnMuZWFydGhkYXRhLm5hc2EuZ292IiwiaWRlbnRpdHlfcHJvdmlkZXIiOiJlZGxfb3BzIiwiYWNyIjoiZWRsIiwiYXNzdXJhbmNlX2xldmVsIjozfQ.n45K2oSM8E_6IAQsB4relFwQjX-QI3D9GlUpBrzaz3P2pQHObUbxyrsz0Y_LnV4a5_r1QfztbEA8nk0DhyF8UOc6RkbRcybxw7LGdRg1lWgGYeDhDf4EQ_HQxZghO4qM3CUmvt4v1Srvep2RjeJXJ6RKaenXgOgpojg3NAqvHRzI52W7X03y_rsNTxONvH9DS0GPdIkHAwvX5ZedDinypGoS6DbnDcVJHhW7UV0gkwGU6XqXWlNS1YmzMBsCds83gCp2JbflhQr81woQeyGIO6YpxriGx-oVhkeM22qEfhEF9BnOfZSm_dRit2x5uWgNrnvL8wXLGtO89HF3abWahw"

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

def update_FIRMS(d:date, sat: Literal["SNPP", "NOAA20", "NOAA21"], product: Literal["SP", "NRT"]):
    """
    Get 1 day of global active fire detections from the FIRMS API.
    If a file already exists for that day, it will be overwritten
    by the new data. 

    sat: 
        satellite name e.g. "SNPP", "NOAA20", "NOAA21" 
    product: 
        "NRT": near real time 
        "SP": standard product 

    If approaching API download rate limits, will back off automatically. 
    """

    if (sat == "NOAA21") and (product == "SP"):
        raise ValueError("NOAA21 standard product is not available. Use NOAA21 NRT.")

    data_dir = os.path.join(settings.dirextdata, "VIIRS", f"FIRMS_VIIRS_{sat}_{product}/")
    status_url = 'https://firms.modaps.eosdis.nasa.gov/mapserver/mapkey_status/?MAP_KEY=' + MAP_KEY
    firms_api = "https://firms.modaps.eosdis.nasa.gov/api/area/csv/" 
    query = f"/VIIRS_{sat}_{product}/world/1/" + d.strftime("%Y-%m-%d")
    url = firms_api + MAP_KEY + query

    retries = 0 
    while retries < N_MAX_RETRIES:

        retries += 1
        if retries >= N_MAX_RETRIES:
            logger.warning(f"Could not download {product} {sat} data for {d}")
            logger.warning(f"Error message: Max retries exceeded.")
            return
        
        resp = pd.read_json(status_url, typ='series')
        count = resp['current_transactions']
        limit = resp['transaction_limit'] 

        if (limit - count > limit * .1):
            try:
                logger.info(f"Downloading {sat} {product} for {d}")
                df = pd.read_csv(url) 
                break
            except Exception as e:
                logger.warning(f"Error while downloading {sat} {product} for {d}: {e}. Retrying download.")
             
        else: 
            logger.warning(
                f"Current FIRMS API transactions ({count}) approaching account limit ({limit}). Sleeping 60 seconds. Retry #{retries}"
            )
            time.sleep(60)
        
    if len(df) < 1: 
        logger.warning(
            f"{product} {sat} data is empty for {d}. This date may be outside range of data availability."
        )
        return 

    daterange = pd.to_datetime(df['acq_date'])
    tst, ted = daterange.min(), daterange.max() 

    if tst.date() != ted.date():
        raise ValueError(f"Unexpected date range for single day file: {tst} to {ted}")

    filename_out = f"FIRMS_VIIRS_{sat}_{product}_{tst.strftime('%Y%m%d')}.csv" 
    downloaded_filepath = os.path.join(data_dir, filename_out)  
    os.makedirs(os.path.dirname(downloaded_filepath), exist_ok=True)
    df.to_csv(downloaded_filepath)
    
    preprocess_input_file(downloaded_filepath)
    return

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


def get_FIRMS_data_availability(sat: Literal["SNPP", "NOAA20", "NOAA21"]):
    """Get current date range of data available via FIRMS API for each VIIRS sensor. 

    Parameters
    ----------
    sat : Literal["SNPP", "NOAA20", "NOAA21"]

    Returns
    -------
    sp_start : pd Timestamp or None 
        First date for which standard product (SP) data is available
        or None if SP data is not available for this satellite
    sp_end : pd Timestamp or None 
        Last date for which standard product (SP) data is available
        or None if SP data is not available for this satellite
    nrt_start : pd Timestamp 
        First date for which near real time (NRT) data is available 
    nrt_end : pd Timestamp 
        Last date for which near real time (NRT) data is available
    
    """
    da_url = 'https://firms.modaps.eosdis.nasa.gov/api/data_availability/csv/' + MAP_KEY + '/all'
    df = pd.read_csv(da_url, index_col='data_id')

    if sat == "SNPP":
        sp_start = pd.to_datetime(df.loc["VIIRS_SNPP_SP"].min_date)
        sp_end = pd.to_datetime(df.loc["VIIRS_SNPP_SP"].max_date)
        nrt_start = pd.to_datetime(df.loc["VIIRS_SNPP_NRT"].min_date)
        nrt_end = pd.to_datetime(df.loc["VIIRS_SNPP_NRT"].max_date)
    elif sat == "NOAA20": 
        sp_start = pd.to_datetime(df.loc["VIIRS_NOAA20_SP"].min_date)
        sp_end = pd.to_datetime(df.loc["VIIRS_NOAA20_SP"].max_date)
        nrt_start = pd.to_datetime(df.loc["VIIRS_NOAA20_NRT"].min_date)
        nrt_end = pd.to_datetime(df.loc["VIIRS_NOAA20_NRT"].max_date)
    elif sat == "NOAA21":
        sp_start = None 
        sp_end = None 
        nrt_start = pd.to_datetime(df.loc["VIIRS_NOAA21_NRT"].min_date)
        nrt_end = pd.to_datetime(df.loc["VIIRS_NOAA21_NRT"].max_date)

    return sp_start, sp_end, nrt_start, nrt_end

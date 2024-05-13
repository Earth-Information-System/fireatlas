""" DataUpdate
This module include functions used to check and update needed data files
"""
import os
import fsspec
import xarray as xr
import tempfile
import pandas as pd
import warnings
import requests

from datetime import date, timedelta

from fireatlas import settings
from fireatlas.FireLog import logger
from fireatlas.preprocess import preprocess_input_file

# ------------------------------------------------------------------------------
# Check external dataset
# ------------------------------------------------------------------------------
def check_VNP14IMGML_avail(year, month, ver="C1.05"):
    """ Check if the monthly VIIRS data exist (C1.05 or C1.04)

    Parameters
    ----------
    year : int
        the year
    month : int
        the month
    """
    # derive monthly file name with path
    t = date(year, month, 1)
    dirFC = os.path.join(settings.dirextdata, "VIIRS", "VNP14IMGML") + "/"
    fnmFC = os.path.join(dirFC, "VNP14IMGML." + t.strftime("%Y%m") + "." + ver + ".txt")

    # check if the file exist
    if settings.fs.exists(fnmFC):
        print("[Yes] Monthly VNP14IMGML exists at ", dirFC)
        fnms = settings.fs.glob(dirFC + "VNP14IMGML." + t.strftime("%Y") + "*.txt")
        mons = max([int(os.path.basename(fnm)[15:17]) for fnm in fnms])
        print(f"The latest VNP14IMGML available month for {year} is {mons}")
    else:
        print("[No] Monthly VNP14IMGML is not available at ", dirFC)
        # if the file does not exist, also find the last available month
        fnms = settings.fs.glob(dirFC + "VNP14IMGML." + t.strftime("%Y") + "*.txt")
        if len(fnms) == 0:
            print("No any monthly VNP14IMGML is not available at ", dirFC)
        else:
            mons = max([int(os.path.basename(fnm)[15:17]) for fnm in fnms])
            print(f"The latest VNP14IMGML available month for {year} is {mons}")


def check_VNP14IMGTDL_avail(year, month, day):
    """ Check if the daily NRT VIIRS data exist

    Parameters
    ----------
    year : int
        the year
    month : int
        the month
    day : int
        the day
    """
    # derive monthly file name with path
    t = date(year, month, day)
    dirFC = os.path.join(settings.dirextdata, "VIIRS", "VNP14IMGTDL") + "/"
    fnmFC = os.path.join(
        dirFC, "SUOMI_VIIRS_C2_Global_VNP14IMGTDL_NRT_" + t.strftime("%Y%j") + ".txt"
    )

    # check if the file exist
    if settings.fs.exists(fnmFC):
        print("[Yes] Daily VNP14IMGTDL exists!")
    else:
        print("[No] Daily VNP14IMGTDL is not available!")

        # if the file does not exist, also find the last available day
        fnms = settings.fs.glob(
            dirFC
            + "SUOMI_VIIRS_C2_Global_VNP14IMGTDL_NRT_"
            + t.strftime("%Y")
            + "*.txt"
        )
        if len(fnms) > 0:
            days = max([int(os.path.basename(fnm)[42:45]) for fnm in fnms])
            tmax = date(t.year, 1, 1) + timedelta(days=days - 1)
            print(
                "The latest VNP14IMGTDL available date is "
                + tmax.strftime("%Y%m%d")
                + "(doy = "
                + tmax.strftime("%j")
                + ")"
            )

def check_VJ114IMGTDL_avail(year, month, day):
    """ Check if the daily NRT VIIRS data exist

    Parameters
    ----------
    year : int
        the year
    month : int
        the month
    day : int
        the day
    """
    # derive monthly file name with path
    t = date(year, month, day)
    dirFC = os.path.join(settings.dirextdata,'VIIRS', "VJ114IMGTDL") + "/"
    fnmFC = os.path.join(
        dirFC, "J1_VIIRS_C2_Global_VJ114IMGTDL_NRT_" + t.strftime("%Y%j") + ".txt"
    )

    # check if the file exist
    if settings.fs.exists(fnmFC):
        print("[Yes] Daily VJ114IMGTDL exists!")
    else:
        print("[No] Daily VJ114IMGTDL is not available!")

        # if the file does not exist, also find the last available day
        fnms = settings.fs.glob(
            dirFC
            + "J1_VIIRS_C2_Global_VJ114IMGTDL_NRT_"
            + t.strftime("%Y")
            + "*.txt"
        )
        if len(fnms) > 0:
            days = max([int(os.path.basename(fnm)[-7:-4]) for fnm in fnms])
            tmax = date(t.year, 1, 1) + timedelta(days=days - 1)
            print(
                "The latest VNP14IMGTDL available date is "
                + tmax.strftime("%Y%m%d")
                + "(doy = "
                + tmax.strftime("%j")
                + ")"
            )
                        
            
def check_GridMET_avail(year, month, day):
    """ Check if the GridMET fm1000 data exist

    Parameters
    ----------
    year : int
        the year
    month : int
        the month
    day : int
        the day
    """
    warnings.simplefilter("ignore")

    t = date(year, month, day)
    dirGridMET = os.path.join(settings.dirextdata, "GridMET") + "/"
    fnmFM1000 = dirGridMET + "fm1000_" + t.strftime("%Y") + ".nc"

    if settings.fs.exists(fnmFM1000):
        print("[Yes] Annual GridMET FM1000 exists!")

        # if the file exist, check whether the date is beyond the time range in the data
        ds = xr.open_dataarray(fnmFM1000)
        tmax = pd.to_datetime(ds["day"][-1].values).date()
        print("The last time in GridMET data is " + tmax.strftime("%Y-%m-%d"))
    else:
        print("[No] Annual GridMET FM1000 is not available!")


def check_data_avail(year, month, day):
    """ A tool used to check the availability of NRT fire and FM1000 data:
    1. monthly suomi VIIRS active fire (VNP14IMGML)
    2. daily suomi VIIRS NRT active fire (VNP14IMGTDL)
    3. annual GridMET FM1000 data

    Parameters
    ----------
    year : int
        the year
    month : int
        the month
    day : int
        the day
    """

    # VIIRS VNP14IMGML data
    check_VNP14IMGML_avail(year, month)

    # VIIRS VNP14IMGTDL data
    check_VNP14IMGTDL_avail(year, month, day)
    
    # VIIRS VJ114IMGTDL data
    check_VJ114IMGTDL_avail(year, month, day)

    # GridMET fm1000 data
    check_GridMET_avail(year, month, day)


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
    dstart = dstart - timedelta(days=5)

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
    dstart = dstart - timedelta(days=5)

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


if __name__ == "__main__":
    """ Do check or update the data
    """
    # check_data_avail(year,month,day)
    check_VNP14IMGML_avail(2021, 1, ver="C1.05")

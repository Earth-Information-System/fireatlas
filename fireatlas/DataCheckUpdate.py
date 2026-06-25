"""DataUpdate
This module include functions used to check and update needed data files
"""

import os
import time
import pandas as pd

from datetime import date
from typing import Literal

from fireatlas import settings
from fireatlas.FireLog import logger

MAP_KEY = "3cb8ce1d0094e20f07f8697df832da3a"
N_MAX_RETRIES = 30


def update_FIRMS(
    d: date, sat: Literal["SNPP", "NOAA20", "NOAA21"], product: Literal["SP", "NRT"]
):
    """
    Get 1 day of global active fire detections from the FIRMS API.
    If a file already exists for that day, it will be overwritten
    by the new data.

    sat:
        satellite name e.g. "SNPP", "NOAA20", "NOAA21"
    product:
        "NRT": near real time. NRT will automatically include both NRT and Ultra Real Time (URT) data.
        "SP": standard product

    If approaching API download rate limits, will back off automatically.

    Returns:
    --------
    downloaded_filepath: str
        Location of downloaded data file
    """

    if (sat == "NOAA21") and (product == "SP"):
        raise ValueError("NOAA21 standard product is not available. Use NOAA21 NRT.")

    data_dir = os.path.join(
        settings.dirextdata, "VIIRS", f"FIRMS_VIIRS_{sat}_{product}/"
    )
    logger.info(f"Running update_FIRMS and saving to {data_dir}")
    status_url = (
        "https://firms.modaps.eosdis.nasa.gov/mapserver/mapkey_status/?MAP_KEY="
        + MAP_KEY
    )
    firms_api = "https://firms.modaps.eosdis.nasa.gov/api/area/csv/"
    query = f"/VIIRS_{sat}_{product}/world/1/" + d.strftime("%Y-%m-%d")
    url = firms_api + MAP_KEY + query

    retries = 0
    while retries < N_MAX_RETRIES:
        retries += 1
        if retries >= N_MAX_RETRIES:
            logger.warning(f"Could not download {product} {sat} data for {d}")
            logger.warning("Error message: Max retries exceeded.")
            return

        resp = pd.read_json(status_url, typ="series")
        count = resp["current_transactions"]
        limit = resp["transaction_limit"]

        if limit - count > limit * 0.1:
            try:
                logger.info(f"Downloading {sat} {product} for {d}")
                df = pd.read_csv(url)
                break
            except Exception as e:
                logger.warning(
                    f"Error while downloading {sat} {product} for {d}: {e}. Retrying download."
                )

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

    daterange = pd.to_datetime(df["acq_date"])
    tst, ted = daterange.min(), daterange.max()

    if tst.date() != ted.date():
        raise ValueError(f"Unexpected date range for single day file: {tst} to {ted}")

    filename_out = f"FIRMS_VIIRS_{sat}_{product}_{tst.strftime('%Y%m%d')}.csv"
    downloaded_filepath = os.path.join(data_dir, filename_out)
    os.makedirs(os.path.dirname(downloaded_filepath), exist_ok=True)
    df.to_csv(downloaded_filepath)
    logger.info(f"Saved df to {downloaded_filepath}")

    return downloaded_filepath


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
    da_url = (
        "https://firms.modaps.eosdis.nasa.gov/api/data_availability/csv/"
        + MAP_KEY
        + "/all"
    )
    df = pd.read_csv(da_url, index_col="data_id")

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

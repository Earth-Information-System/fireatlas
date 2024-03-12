import os
import uuid
import pandas as pd
from typing import Literal, Optional
from shapely import to_geojson, from_geojson
import sys

from tqdm import tqdm
import rasterio
import rasterio.warp

import FireConsts
import FireClustering
from FireLog import logger
import FireIO
import FireTime
from FireTypes import Region, TimeStep
from utils import timed


def preprocessed_region_filename(region: Region):
    return os.path.join(FireConsts.dirprpdata, f"{region[0]}.json")


@timed
def preprocess_region(region: Region, force=False):
    from FireMain import maybe_remove_static_sources

    output_filepath = preprocessed_region_filename(region)
    if not force and os.path.exists(output_filepath):
        logger.info("Preprocessing has already occurred for this region.")
        logger.debug("Use `force=True` to rerun this preprocessing step.")
        return output_filepath

    # make path if necessary
    os.makedirs(os.path.dirname(output_filepath), exist_ok=True)

    region = maybe_remove_static_sources(region, FireConsts.dirextdata)

    with open(output_filepath, "w") as f:
        f.write(to_geojson(region[1], indent=2))
    
    return output_filepath

@timed
def read_region(region: Region):
    filepath = preprocessed_region_filename(region)

    with open(filepath, "r") as f:
        shape = from_geojson(f.read())
    return (region[0], shape)


def preprocessed_landcover_filename(filename="nlcd_export_510m_simplified"):
    return os.path.join(FireConsts.dirprpdata, f"{filename}_latlon.tif")


@timed
def preprocess_landcover(filename="nlcd_export_510m_simplified", force=False):
    # if landcover output already exists, exit early so we don't reprocess
    output_filepath = preprocessed_landcover_filename(filename)
    if not force and os.path.exists(output_filepath):
        logger.info("Preprocessing has already occurred for this landcover file.")
        logger.debug("Use `force=True` to rerun this preprocessing step.")
        return output_filepath
    
    fnmLCT = os.path.join(FireConsts.dirextdata, "NLCD", f"{filename}.tif")
    
    # make nested path if necessary
    os.makedirs(os.path.dirname(output_filepath), exist_ok=True)

    dst_crs = f"EPSG:4326"

    with rasterio.open(fnmLCT) as src:
        transform, width, height = rasterio.warp.calculate_default_transform(
            src.crs, dst_crs, src.width, src.height, *src.bounds
        )
        kwargs = src.meta.copy()
        kwargs.update(
            {"crs": dst_crs, "transform": transform, "width": width, "height": height}
        )

        with rasterio.open(output_filepath, "w", **kwargs) as dst:
            for i in range(1, src.count + 1):
                rasterio.warp.reproject(
                    source=rasterio.band(src, i),
                    destination=rasterio.band(dst, i),
                    src_transform=src.transform,
                    src_crs=src.crs,
                    dst_transform=transform,
                    dst_crs=dst_crs,
                    resampling=rasterio.warp.Resampling.nearest,
                )
    return output_filepath


def preprocessed_filename(
    t: TimeStep,
    sat: Literal["NOAA20", "SNPP", "VIIRS", "TESTING123"],
    *,
    region: Optional[Region] = None,
    suffix="",
):
    return os.path.join(
        FireConsts.dirprpdata,
        *([] if region is None else [region[0]]),
        sat,
        f"{t[0]}{t[1]:02}{t[2]:02}_{t[3]}{suffix}.txt",
    )


def NRT_filepath(t: TimeStep, sat: Literal["SNPP", "NOAA20"]):
    """ Filepath for NRT VIIRS data

    Parameters
    ----------
    t : tuple, (int,int,int,str)
        the year, month, day and 'AM'|'PM' during the initialization
    sat: Literal["SNPP", "NOAA20"]
        which satellite to use

    Returns
    -------
    filepath : str
        Path to input data or None if file does not exist
    """
    if sat == "SNPP":
        filepath = FireIO.VNP14IMGTDL_filepath(t)
    elif sat == "NOAA20":
        filepath = FireIO.VJ114IMGTDL_filepath(t)
    else:
        raise ValueError("Please set SNPP or NOAA20 for sat")
    return filepath


def monthly_filepath(t: TimeStep, sat: Literal["NOAA20", "SNPP"]):
    """ Filepath for monthly VIIRS data

    Parameters
    ----------
    t : tuple, (int,int,int,str)
        the year, month, day and 'AM'|'PM' during the initialization
    sat: Literal["SNPP", "NOAA20"]
        which satellite to use

    Returns
    -------
    filepath : str
        Path to input data or None if file does not exist
    """
    if sat == "SNPP":
        filepath = FireIO.VNP14IMGML_filepath(t)
    elif sat == "NOAA20":
        filepath = FireIO.VJ114IMGML_filepath(t)
    else:
        raise ValueError("please set SNPP or NOAA20 for sat")
    return filepath


def check_preprocessed_file(list_of_ts, sat):
    """Before running preprocess_monthly_file, check if the file already exists
    for that satellite using a list of time steps

    Parameters
    ----------
    list_of_ts : list of time tuples
    sat : str, 'SNPP' or 'NOAA20'

    Returns
    -------
    list of unique combos of years and months that need to be processed
    """
    # check that there's viirs data for these dates and if not, keep track
    needs_processing = []
    for ts in list_of_ts:
        filepath = preprocessed_filename(ts, sat)
        if not os.path.exists(filepath):
            needs_processing.append(ts)

    # get unique combos of year and month. make sure sorted
    unique_ym = list(set([(ts[0], ts[1]) for ts in needs_processing]))
    unique_ym = sorted(unique_ym, key=lambda x: (x[0], x[1]))
    return unique_ym


@timed
def preprocess_input_file(filepath: str):
    """
    Preprocess monthly or daily NRT file of fire location data.

    NOTE: Satellite is deduced from the filepath.

    Parameters
    ----------
    filepath : str
        Path to input data. Can be local or s3.

    Returns
    -------
    output_paths : list[str]
        List of filepaths that this function has written to.
    """
    if filepath is None:
        raise ValueError("Please provide a valid filepath")
    
    logger.info(f"preprocessing {filepath.split('/')[-1]}")

    if "VNP14IMGTDL" in filepath:
        sat = "SNPP"
        df = FireIO.read_VNP14IMGTDL(filepath)
    elif "VJ114IMGTDL" in filepath:
        sat = "NOAA20"
        df = FireIO.read_VJ114IMGTDL(filepath)
    elif "VNP14IMGML" in filepath:
        sat = "SNPP"
        df = FireIO.read_VNP14IMGML(filepath)
        df = df.loc[df["Type"] == 0]  # type filtering
    elif "VJ114IMGML" in filepath:
        sat = "NOAA20"
        df = FireIO.read_VJ114IMGML(filepath)
        df = df.loc[df["mask"] >= 7]
    else:
        raise ValueError("please set SNPP or NOAA20 for sat")

    # set ampm
    df = FireIO.AFP_setampm(df)

    # add the satellite information
    df["Sat"] = sat
    df["input_filename"] = filepath.split("/")[-1]

    # return selected columns
    df = df[["Lat", "Lon", "FRP", "Sat", "DT", "DS", "input_filename", "datetime", "ampm"]]

    output_paths = []

    # groupby days and if there are more than 1 days, include a progress bar
    gb = df.groupby(df["datetime"].dt.date)
    if gb.ngroups > 1:
        gb = tqdm(gb, "Processing days", file=sys.stdout)

    for day, data in gb:
        for ampm in ["AM", "PM"]:
            time_filtered_df = data.loc[df["ampm"] == ampm]

            output_filepath = preprocessed_filename(
                (day.year, day.month, day.day, ampm), sat
            )

            # make nested path if necessary
            os.makedirs(os.path.dirname(output_filepath), exist_ok=True)

            # save active pixels at this time step (day and ampm filter)
            time_filtered_df.to_csv(output_filepath, index=False)
            
            output_paths.append(output_filepath)
    
    return output_paths


def preprocess_monthly_file(t: TimeStep, sat: Literal["NOAA20", "SNPP"]):
    filepath = monthly_filepath(t, sat)
    return preprocess_input_file(filepath)


def preprocess_NRT_file(t: TimeStep, sat: Literal["NOAA20", "SNPP"]):
    filepath = NRT_filepath(t, sat)
    return preprocess_input_file(filepath)


@timed
def read_preprocessed_input(
    t: TimeStep,
    sat: Literal["NOAA20", "SNPP", "VIIRS", "TESTING123"],
):
    filename = preprocessed_filename(t, sat)
    df = pd.read_csv(filename)
    return df


@timed
def read_preprocessed(
    t: TimeStep,
    sat: Literal["NOAA20", "SNPP", "VIIRS", "TESTING123"],
    region: Region,
):
    filename = preprocessed_filename(t, sat, region=region)
    df = pd.read_csv(filename).set_index("uuid").assign(t=FireTime.t2dt(t))
    return df


@timed
def preprocess_region_t(t: TimeStep, sensor: Literal['SNPP', 'NOAA20', 'VIIRS', "TESTING123"], region: Region, force=False):
    
    # if regional output already exists, exit early so we don't reprocess
    output_filepath = preprocessed_filename(t, sensor, region=region)
    if not force and os.path.exists(output_filepath):
        logger.info(
            "Preprocessing has already occurred for this combination of "
            "timestep, sensor, and region." 
        )
        logger.debug("Use `force=True` to rerun this preprocessing step.")
        return output_filepath
    
    region = read_region(region)
    logger.info(
        f"filtering and clustering {t[0]}-{t[1]}-{t[2]} {t[3]}, {sensor}, {region[0]}"
    )
    if sensor == "VIIRS":
        df = pd.concat(
            [
                read_preprocessed_input(t, sat="SNPP"),
                read_preprocessed_input(t, sat="NOAA20"),
            ],
            ignore_index=True,
        )
    else:
        df = read_preprocessed_input(t, sat=sensor)

    # read in the preprocessed region
    region = read_region(region)

    # do regional filtering
    shp_Reg = FireIO.get_reg_shp(region[1])
    df = FireIO.AFP_regfilter(df, shp_Reg)

    columns = ["Lat", "Lon", "FRP", "Sat", "DT", "DS", "input_filename", "datetime", "ampm", "x", "y"]

    if not df.empty:
        # return selected columns
        df = df[columns]

        # do preliminary clustering using new active fire locations (assign cid to each pixel)
        df = FireClustering.do_clustering(df, FireConsts.CONNECTIVITY_CLUSTER_KM)
        
        # assign a uuid to each pixel and put it as the first column
        df.insert(0, "uuid" , [uuid.uuid4() for _ in range(len(df.index))])
    else: 
        # make a dummy DataFrame with the right columns so that we know later that 
        # we don't need to do this step again.
        df = pd.DataFrame(columns = ["uuid", *columns, "initial_cid"])

    # make nested path if necessary
    os.makedirs(os.path.dirname(output_filepath), exist_ok=True)

    df.to_csv(output_filepath, index=False)

    return output_filepath
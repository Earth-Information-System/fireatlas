import os
import uuid
import pandas as pd
from typing import Literal, Optional
from shapely import to_geojson, from_geojson

import rasterio
import rasterio.warp

import FireConsts
import FireClustering
from FireLog import logger
import FireIO
import FireTime
from FireTypes import Region, TimeStep
from utils import timed


OUTPUT_DIR = FireConsts.dirprpdata


def preprocessed_region_filename(region: Region):
    return os.path.join(OUTPUT_DIR, f"{region[0]}.json")


@timed
def preprocess_region(region: Region):
    from FireMain import maybe_remove_static_sources

    output_filepath = preprocessed_region_filename(region)
    
    # make path if necessary
    os.makedirs(os.path.dirname(output_filepath), exist_ok=True)

    region = maybe_remove_static_sources(region, FireConsts.dirextdata)

    with open(output_filepath, "w") as f:
        f.write(to_geojson(region[1], indent=2))


@timed
def read_region(region: Region):
    filepath = preprocessed_region_filename(region)

    with open(filepath, "r") as f:
        shape = from_geojson(f.read())
    return (region[0], shape)


def preprocessed_landcover_filename(filename="nlcd_export_510m_simplified"):
    return os.path.join(OUTPUT_DIR, f"{filename}_latlon.tif")


@timed
def preprocess_landcover(filename="nlcd_export_510m_simplified"):
    fnmLCT = os.path.join(FireConsts.dirextdata, "NLCD", f"{filename}.tif")
    
    output_filepath = preprocessed_landcover_filename(filename)
    
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


def preprocessed_filename(
    t: TimeStep,
    sat: Literal["NOAA20", "SNPP", "VIIRS", "TESTING123"],
    *,
    region: Optional[Region] = None,
    suffix="",
    preprocessed_data_dir=FireConsts.dirprpdata,
):
    return os.path.join(
        preprocessed_data_dir,
        sat,
        *([] if region is None else [region[0]]),
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
    Preprocess monthly or daily NRT file of fire location data

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

    # return selected columns
    df = df[["Lat", "Lon", "FRP", "Sat", "DT", "DS", "datetime", "ampm"]]

    output_paths = []

    for day, data in df.groupby(df["datetime"].dt.date):
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

@timed
def read_preprocessed(
    t: TimeStep,
    sat: Literal["NOAA20", "SNPP", "VIIRS", "TESTING123"],
    *,
    region: Optional[Region] = None,
):
    filename = preprocessed_filename(t, sat, region=region)
    if os.path.exists(filename):
        df = pd.read_csv(filename)
        if region:
            df = df.set_index("uuid").assign(t=FireTime.t2dt(t))
        return df

@timed
def preprocess_region_t(t: TimeStep, sensor: Literal["VIIRS", "TESTING123"], region: Region):
    logger.info(
        f"filtering and clustering {t[0]}-{t[1]}-{t[2]} {t[3]}, {sensor}, {region[0]}"
    )
    if sensor == "VIIRS":
        df = pd.concat(
            [
                read_preprocessed(t, sat="SNPP"),
                read_preprocessed(t, sat="NOAA20"),
            ],
            ignore_index=True,
        )
    else:
        df = read_preprocessed(t, sat=sensor)

    # read in the preprocessed region
    region = read_region(region)

    # if regional output already exists, exit early so we don't reprocess
    output_filepath = preprocessed_filename(t, sensor, region=region)
    if os.path.exists(output_filepath):
        return output_filepath

    # do regional filtering
    shp_Reg = FireIO.get_reg_shp(region[1])
    df = FireIO.AFP_regfilter(df, shp_Reg)

    columns = ["Lat", "Lon", "FRP", "Sat", "DT", "DS", "datetime", "ampm", "x", "y"]

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
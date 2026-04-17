import os
import uuid
import fsspec
import pandas as pd
import datetime as dt
from typing import Literal, Optional
from shapely import to_geojson, from_geojson
import sys

from tqdm import tqdm
import rasterio
import rasterio.warp

from fireatlas.FireLog import logger
from fireatlas.FireTypes import Region, TimeStep, Location
from fireatlas.utils import timed
from fireatlas.FireClustering import do_clustering
from fireatlas.FireTime import t_generator, t2dt, t_nb, t_nd, t_nm
from fireatlas import FireIO, FireMain, settings, FireTime


def preprocessed_region_filename(region: Region, location: Location = None):
    return os.path.join(
        settings.get_path(location), settings.PREPROCESSED_DIR, region[0], f"{region[0]}.json"
    )


@timed
def preprocess_region(region: Region, force=False):

    output_filepath = preprocessed_region_filename(region, location="local")
    if not force and os.path.exists(output_filepath):
        logger.info("Preprocessing has already occurred for this region.")
        logger.debug("Use `force=True` to rerun this preprocessing step.")
        return output_filepath

    # make path if necessary
    os.makedirs(os.path.dirname(output_filepath), exist_ok=True)

    region = FireMain.maybe_remove_static_sources(region)

    with open(output_filepath, "w") as f:
        f.write(to_geojson(region[1], indent=2))

    return output_filepath


@timed
def read_region(region: Region, location: Location = None):
    filepath = preprocessed_region_filename(region, location=location)

    # use fsspec here b/c it could be s3 or local
    with fsspec.open(filepath, "r") as f:
        shape = from_geojson(f.read())
    return (region[0], shape)


def preprocessed_landcover_filename(
    filename="nlcd_export_510m_simplified",
    location: Location = None,
):
    return os.path.join(settings.get_path(location), settings.PREPROCESSED_DIR, f"{filename}_latlon.tif")


@timed
def preprocess_landcover(filename="nlcd_export_510m_simplified", force=False):
    # if landcover output already exists, exit early so we don't reprocess
    output_filepath = preprocessed_landcover_filename(filename, location="local")
    if not force and os.path.exists(output_filepath):
        logger.info("Preprocessing has already occurred for this landcover file.")
        logger.debug("Use `force=True` to rerun this preprocessing step.")
        return output_filepath

    fnmLCT = os.path.join(settings.dirextdata, "NLCD", f"{filename}.tif")

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
    sat: Optional[Literal["NOAA20", "NOAA21", "SNPP"]] = None,
    region: Optional[Region] = None,
    suffix="",
    location: Location = None
):
    if sat is None:
        sat = settings.FIRE_SOURCE

    return os.path.join(
        settings.get_path(location),
        settings.PREPROCESSED_DIR,
        *([] if region is None else [region[0]]),
        sat,
        f"{t[0]}{t[1]:02}{t[2]:02}_{t[3]}{suffix}.txt",
    )


def NRT_filepath(t: TimeStep, sat: Literal["SNPP", "NOAA20", "NOAA21"]):
    """Filepath for NRT VIIRS data

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

def FIRMS_NRT_filepath(t: TimeStep, sat: Literal["SNPP", "NOAA20", "NOAA21"]):
    """Filepath for daily NRT VIIRS data from FIRMS

    Parameters
    ----------
    t : tuple, (int,int,int,str)
        the year, month, day and 'AM'|'PM' during the initialization
    sat: Literal["SNPP", "NOAA20", "NOAA21"]
        which satellite to use

    Returns
    -------
    filepath : str
        Path to input data or None if file does not exist
    """
    if sat == "SNPP":
        filepath = FireIO.FIRMS_VIIRS_SNPP_NRT_filepath(t)
    elif sat == "NOAA20":
        filepath = FireIO.FIRMS_VIIRS_NOAA20_NRT_filepath(t)
    elif sat == "NOAA21":
        filepath = FireIO.FIRMS_VIIRS_NOAA21_NRT_filepath(t)
    else:
        raise ValueError("Please set SNPP, NOAA20, or NOAA21 for sat")

    if not settings.fs.exists(filepath):
        return None
    else:
        return filepath

def FIRMS_SP_filepath(t: TimeStep, sat: Literal["SNPP", "NOAA20", "NOAA21"]):
    """Filepath for daily SP VIIRS data from FIRMS

    Parameters
    ----------
    t : tuple, (int,int,int,str)
        the year, month, day and 'AM'|'PM' during the initialization
    sat: Literal["SNPP", "NOAA20", "NOAA21"]
        which satellite to use

    Returns
    -------
    filepath : str
        Path to input data or None if file does not exist
    """
    if sat == "SNPP":
        filepath = FireIO.FIRMS_VIIRS_SNPP_SP_filepath(t)
    elif sat == "NOAA20":
        filepath = FireIO.FIRMS_VIIRS_NOAA20_SP_filepath(t)
    else:
        raise ValueError("Please set SNPP or NOAA20 for sat")

    if not settings.fs.exists(filepath):
        return None
    else:
        return filepath


def monthly_filepath(t: TimeStep, sat: Literal["NOAA20", "SNPP"]):
    """Filepath for monthly VIIRS data

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
        raise ValueError(f"sat={sat} not recognized: please set SNPP or NOAA20 for sat")
    return filepath


def check_preprocessed_file(
    tst: TimeStep,
    ted: TimeStep,
    sat: Literal["SNPP", "NOAA20", "NOAA21"],
    freq: Literal["monthly", "NRT"] = "monthly",
    location: Location = None,
):
    """Before running preprocess_monthly_file, check if the preprocessed files already exist
    for that satellite using a list of time steps

    Parameters
    ----------
    tst : tuple, (int,int,int,str)
        the year, month, day and 'AM'|'PM' to start checking for files
    ted : tuple, (int,int,int,str)
        the year, month, day and 'AM'|'PM' to end checking for files
    sat: Literal["SNPP", "NOAA20", "NOAA21"]
        which satellite to use
    freq: Literal["monthly", "NRT"]
        which files to use - monthly or daily (NRT)
    location: optional Literal["s3", "local"] 
        where to check for files

    Returns
    -------
    list of unique combos of years and months (and days if NRT) that need to be processed
    """
    location = location or settings.READ_LOCATION
    fs = fsspec.filesystem(location, use_listings_cache=False)
    # check that there is preprocessed data for these dates and if not, keep track

    needs_processing = []
    for t in t_generator(tst, ted):
        filepath = preprocessed_filename(t, sat=sat, location=location)
        if not fs.exists(filepath):
            needs_processing.append(t)

    if freq == "monthly":
        return list(set([(t[0], t[1]) for t in needs_processing]))
    else:
        return list(set([(t[0], t[1], t[2], t[3]) for t in needs_processing]))


def get_date_from_input_filename(filepath: str):
    """Extract year, month, and optionally day from input filename.
    
    Supports MODAPS (monthly and daily) and FIRMS formats.
    
    Returns: (year, month, day) where day is None for monthly files
    """
    filename = os.path.basename(filepath)
    day = None
    if "14IMGML" in filename:
        # monthly file e.g. VNP14IMGML.201201.C2.05.csv
        datestring = filename.split(".")[1]
        year = datestring[:4]
        month = datestring[-2:]
    elif "IMGTDL" in filename:
        # daily file e.g. SUOMI_VIIRS_C2_Global_VNP14IMGTDL_NRT_2012360.txt
        datestring = filename.split("_")[6]
        year = datestring[:4]
        julian_day = datestring[4:7]
        d = dt.date(int(year), 1, 1) + dt.timedelta(days=int(julian_day) - 1)
        year = d.year
        month = d.month
        day = d.day
    elif "FIRMS_VIIRS" in filename:
        # FIRMS downloaded daily file, e.g. FIRMS_VIIRS_SNPP_SP_20120101.csv
        datestring = filename.split("_")[4]
        year = datestring[:4]
        month = datestring[4:6]
        day = datestring[6:8]
    else:
        raise ValueError(f"Could not infer date from input filename {filename}")
    
    return (int(year), int(month), day)


@timed
def preprocess_input_file(filepath: str, filepath_prev: str | None, filepath_next: str | None):
    """
    Preprocess monthly or daily NRT file of fire location data.

    NOTE: Input files are named by UTC date or month. Output files are named
    with the aprox local solar date/time for each observation. Pixels in the
    output file YYYYMMDD_AM.txt are for the AM overpass for that date as defined
    by aprox local solar time. This means that they may come from the previous
    or next UTC date.

    NOTE: Satellite is deduced from the filepath.

    Parameters
    ----------
    filepath : str
        Path to input data. Can be local or s3.
    filepath_prev : str | None
        Path to input data for previous timestep. If None, this function will simply not
        check the input file for the previous UTC timestep. This can lead to
        missing values that are within the current timestep in local time but not
        UTC time.
    filepath_next : str | None
        Path to input data for next timestep.

    Returns
    -------
    output_paths : list[str]
        List of filepaths that this function has written to.
    """
    if filepath is None:
        raise ValueError("Please provide a valid filepath")

    logger.info(f"preprocessing {filepath.split('/')[-1]}")
    dfs = []
    sat = None
    for f in [filepath_prev, filepath, filepath_next]:
        if not f:
            # it can be valid to have no prev or next file
            # move on to next file
            continue

        # read file
        if "VNP14IMGTDL" in f:
            sat = "SNPP"
            df = FireIO.read_VNP14IMGTDL(f)
        elif "VJ114IMGTDL" in f:
            sat = "NOAA20"
            df = FireIO.read_VJ114IMGTDL(f)
        elif "VNP14IMGML" in f:
            sat = "SNPP"
            df = FireIO.read_VNP14IMGML(f)
            df = df.loc[df["Type"] == 0]
        elif "VJ114IMGML" in f:
            sat = "NOAA20"
            df = FireIO.read_VJ114IMGML(f)
        elif "FIRMS_VIIRS_SNPP_NRT" in f:
            sat = "SNPP"
            df = FireIO.read_FIRMS_VIIRS_NRT(f)
        elif "FIRMS_VIIRS_SNPP_SP" in f:
            sat = "SNPP"
            df = FireIO.read_FIRMS_VIIRS_SP(f)
            df = df.loc[df["Type"] == 0]
            # Type filter: inferred hot spot type == presumed vegetation fire
        elif "FIRMS_VIIRS_NOAA20_NRT" in f:
            sat = "NOAA20"
            df = FireIO.read_FIRMS_VIIRS_NRT(f)
        elif "FIRMS_VIIRS_NOAA20_SP" in f:
            sat = "NOAA20"
            df = FireIO.read_FIRMS_VIIRS_SP(f)
            df = df.loc[df["Type"] == 0]
            # Type filter: inferred hot spot type == presumed vegetation fire
        elif "FIRMS_VIIRS_NOAA21_NRT" in f:
            sat = "NOAA21"
            df = FireIO.read_FIRMS_VIIRS_NRT(f)
        else:
            raise ValueError(f"Filepath {f} not recognized during preprocessing.")

        # add file retrieval information
        df["input_filename"] = f.split("/")[-1]

        dfs.append(df)

    df = pd.concat(dfs)

    # Convert from UTC to aprox local time
    df["local_datetime"] = (pd.to_timedelta(df.Lon / 15, unit="hours") + df["datetime"])

    # get the date of the main input file
    query_year, query_month, query_day = get_date_from_input_filename(filepath)
    # Select only observations that are on the date of the main input file in the local timezone
    if ("VJ114IMGML" in filepath) or ("VNP14IMGML" in filepath):
        df = df[(df.local_datetime.dt.year == query_year) & (df.local_datetime.dt.month == query_month)]
    else:
        df = df[(df.local_datetime.dt.day == query_day) &
                (df.local_datetime.dt.month == query_month) &
                (df.local_datetime.dt.year == query_year)
            ]

    df = FireIO.AFP_setampm(df)
    df["Sat"] = sat
    # groupby days and if there are more than 1 days, include a progress bar
    gb = df.groupby(df["local_datetime"].dt.date)

    # return selected columns

    if settings.FIRE_NRT == True: # preserve version code if working with NRT data
        df = df[
            ["Lat", "Lon", "FRP", "Sat", "DT", "DS", "input_filename", "datetime", "ampm", "version"]
        ]
    else:
        df = df[
            ["Lat", "Lon", "FRP", "Sat", "DT", "DS", "input_filename", "datetime", "ampm"]
        ]

    output_paths = []

    if gb.ngroups > 1:
        gb = tqdm(gb, "Processing days", file=sys.stdout)

    for day, data in gb:
        for ampm in ["AM", "PM"]:
            time_filtered_df = data.loc[df["ampm"] == ampm]

            output_filepath = preprocessed_filename(
                (day.year, day.month, day.day, ampm), sat=sat, location="local"
            )

            # make nested path if necessary
            os.makedirs(os.path.dirname(output_filepath), exist_ok=True)

            # save active pixels at this time step (day and ampm filter)
            time_filtered_df.to_csv(output_filepath, index=False)

            output_paths.append(output_filepath)

    return output_paths


def preprocess_monthly_file(t: TimeStep, sat: Literal["NOAA20", "SNPP"]):
    filepath = monthly_filepath(t, sat=sat)
    return preprocess_input_file(filepath, None, None)


def preprocess_NRT_file(t: TimeStep, sat: Literal["NOAA20", "SNPP"]):
    filepath = NRT_filepath(t, sat=sat)
    return preprocess_input_file(filepath, None, None)


def preprocess_daily_file(filepath, t: TimeStep, sat: Literal["SNPP", "NOAA20", "NOAA21"]):
    """Find previous and next daily input files, then preprocess this timestep.
    Prefers FIRMS standard product (SP) over FIRMS NRT if we have both.

    Parameters
    ----------
    filepath : str
        path to the daily input file to be preprocessed
    t : TimeStep
        time of the input daily file
    sat : Literal["SNPP", "NOAA20", "NOAA21"]
        which satellite the input file is from

    Returns
    -------
    output_paths : list[str]
        List of filepaths that preprocess_input_file function has written to.
    """
    day_prev = t_nd(t, "previous")
    day_next = t_nd(t, "next")

    if settings.FIRE_NRT == True:
        filepath_prev = FIRMS_NRT_filepath(day_prev, sat)
        filepath_next = FIRMS_NRT_filepath(day_next, sat=sat)

    else:
        filepath_prev = FIRMS_SP_filepath(day_prev, sat=sat)
        filepath_next = FIRMS_SP_filepath(day_next, sat=sat)

    return preprocess_input_file(filepath, filepath_prev, filepath_next)


@timed
def read_preprocessed_input(
    t: TimeStep,
    sat: Literal["NOAA20", "NOAA21", "SNPP"],
    location: Location = None,
):
    filename = preprocessed_filename(t, sat=sat, location=location)
    df = pd.read_csv(filename)
    return df


@timed
def read_preprocessed(
    t: TimeStep,
    region: Region,
    location: Location = None,
):
    filename = preprocessed_filename(t, region=region, location=location)
    df = pd.read_csv(filename).set_index("uuid").assign(t=t2dt(t))
    df["datetime"] = pd.to_datetime(df["datetime"], format='ISO8601')
    return df


@timed
def preprocess_region_t(
    t: TimeStep,
    region: Region,
    force: bool = False,
    read_location: Location = None,
    read_region_location: Location = None,
):

    # if regional output already exists, exit early so we don't reprocess
    output_filepath = preprocessed_filename(t, region=region, location="local")
    if not force and os.path.exists(output_filepath):
        logger.info(
            "Preprocessing has already occurred for this combination of "
            "timestep, sensor, and region."
        )
        logger.debug("Use `force=True` to rerun this preprocessing step.")
        return output_filepath

    # read in the preprocessed region
    region = read_region(region, location=read_region_location or read_location)
    source = settings.FIRE_SOURCE
    logger.info(
        f"filtering and clustering {t[0]}-{t[1]}-{t[2]} {t[3]}, {source}, {region[0]}"
    )
    if source == "VIIRS":
        dfs = []
        for sat in ["SNPP", "NOAA20", "NOAA21"]:
            try:
                dfs.append(read_preprocessed_input(t, sat=sat, location=read_location))
            except (FileNotFoundError, pd.errors.EmptyDataError) as e:
                logger.info(f"{sat} file or data not available at {t=}: '{str(e)}'")
        if len(dfs) == 0:
            raise ValueError(f"NOAA20, NOAA21, and SNPP files are not available for {t=}")
        else:
            df = pd.concat(dfs, ignore_index=True)
    else:
        df = read_preprocessed_input(t, sat=source, location=read_location)

    # do regional filtering
    shp_Reg = FireIO.get_reg_shp(region[1])
    df = FireIO.AFP_regfilter(df, shp_Reg)

    columns = [
        "Lat",
        "Lon",
        "FRP",
        "Sat",
        "DT",
        "DS",
        "input_filename",
        "datetime",
        "ampm",
        "x",
        "y",
    ]

    if settings.FIRE_NRT == True:
        columns.append("version") # preserve version type with NRT data

    if not df.empty:
        # return selected columns
        df = df[columns]

        # do preliminary clustering using new active fire locations (assign cid to each pixel)
        df = do_clustering(df, settings.CONNECTIVITY_CLUSTER_KM)

        # assign a uuid to each pixel and put it as the first column
        df.insert(0, "uuid", [uuid.uuid4() for _ in range(len(df.index))])
    else:
        # make a dummy DataFrame with the right columns so that we know later that
        # we don't need to do this step again.
        df = pd.DataFrame(columns=["uuid", *columns, "initial_cid"])

    # make nested path if necessary
    os.makedirs(os.path.dirname(output_filepath), exist_ok=True)

    df.to_csv(output_filepath, index=False)

    return output_filepath

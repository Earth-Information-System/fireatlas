""" FireIO
This module include functions used to read and save data
"""

# ------------------------------------------------------------------------------
# %% Read and filter active fire data
# ------------------------------------------------------------------------------

# Try to read a Geopandas file several times. Sometimes, the read fails the
# first time for mysterious reasons.
import s3fs
import os
import pandas as pd
import numpy as np
import geopandas as gpd
import shapely
import rasterio
import pyproj
import fsspec
import pickle
import xarray as xr
import warnings
from shapely.geometry import Point, Polygon
from datetime import datetime, date

from fireatlas.FireLog import logger
from fireatlas.FireTypes import TimeStep
from fireatlas import FireTime, settings


def preprocess_polygon(
    polygon_gdf,
    data_source,
    id_col,
    firename_col=None,
    startdate_col=None,
    enddate_col=None,
    geometry_col=None,
):
    """
    Process a geodataframe to ensure it has the correct column names, date types, and crs.
    Also buffers geometry by 375m to account for VIIRS pixel size.
    Need to specify data source and column used for ID. Each fire dataset has a potentially different ID so need to
    include this to make the ID meaningful.
    I expect a geodataframe with columns for fire name, start date, end date, and geometry.
    Names come from CALFIRE's FRAP fire perimeter dataset.
    If any column names differ from expected, user should specify which columns to use.

    Parameters
    ----------
    polygon_gdf : geopandas.GeoDataFrame
        A geodataframe containing fire perimeters.
    data_source : str, data source e.g. "FRAP", "NIFC WFIGS", etc.
    id_col : str, name of the id column used to record id e.g. "INC_NUM", "IRWIN", etc.
    firename_col : str, column name for fire name if different from expected
    startdate_col : str, column name for start date if different from expected
    enddate_col : str, column name for end date if different from expected
    geometry : str, column name for geometry if different from expected

    Returns
    -------
    geopandas.GeoDataFrame
        A geodataframe with the correct column names, date types, and crs.
    """

    # test if it is a geodataframe, if not return an error message
    if not isinstance(polygon_gdf, gpd.GeoDataFrame):
        raise TypeError("Input should be a GeoDataFrame")

    # buffer polygon by VIIRS res. need to convert to crs in meters
    polygon_m = polygon_gdf.geometry.to_crs(9311)
    polygon_m = polygon_m.buffer(settings.VIIRSbuf)

    # reset geometry and convert to lat/lon
    polygon_gdf = polygon_gdf.set_geometry(polygon_m).to_crs(4326)

    # select the columns that are needed
    colnames = ["FIRE_NAME", "ALARM_DATE", "CONT_DATE", "geometry"]
    updated_colnames = [firename_col, startdate_col, enddate_col, geometry_col]

    # coalesces the updated column names with the original column names, taking the non-None values, select cols
    updated_colnames = [
        new_colname or orig_colname
        for new_colname, orig_colname in zip(updated_colnames, colnames)
    ]
    perim_processed = polygon_gdf[updated_colnames]
    # renames columns to what you want
    perim_processed.columns = colnames

    # add in data source, id name, and fire id
    if id_col in polygon_gdf.columns:
        fire_id = polygon_gdf[id_col]
    else:
        raise ValueError(f"The column {id_col} does not exist in the data")
    perim_processed.insert(1, "FIRE_ID", fire_id)
    perim_processed.insert(2, "ID_NAME", id_col)
    perim_processed.insert(3, "DATA_SOURCE", data_source)

    # finally, make sure the date columns are formatted correctly
    perim_processed = perim_processed.copy()
    date_cols = colnames[1:3]
    for d in date_cols:
        if perim_processed[d].dtype != "datetime64[ns, UTC]":
            perim_processed[d] = pd.to_datetime(perim_processed[d])

    return perim_processed


def gpd_read_file(filename, parquet=False, **kwargs):
    itry = 0
    maxtries = 5
    fun = gpd.read_parquet if parquet else gpd.read_file
    while itry < maxtries:
        try:
            dat = fun(filename, **kwargs)
            return dat
        except Exception as e:
            itry += 1
            print(f"Attempt {itry}/{maxtries} failed.")
            if not itry < maxtries:
                raise e


def viirs_pixel_size(sample, band="i", rtSCAN_ANGLE=False):
    """calculate approximate size of i-band (375-m) or m-band (750-m) VIIRS pixel
        Adapted from L. Giolio's IDL code
    Usage: DS, DT = viirs_pixel_size(200,band='m',rtSCAN_ANGLE=False)

    Parameters
    ----------
    sample : int
        sample number
    band : str, 'i'|'m'
        i (default) or m band
    rtSCAN_ANGLE : bool, default = False
        flag to return scan angle

    Returns
    -------
    DT : float
        length in along-track dimension [km]
    DS : float
        length in along-scan dimension [km]
    """

    # set constants
    earth_radius = 6371.0  # earth radius [km]
    h = 824.0  # SUOMI-NPP orbit altitude [km]
    pt = 0.361  # nadir pixel resolution [km]
    scan_rate = 56.06 / 6304  # scan rate (deg)
    st1, st2, st3 = 1184, 1920, 3200  # sample tier steps
    sb1, sb2, sb3 = 0, 3552, 5024  # sample base for each tier
    ps_const = 0.128  # constant to convert zone number to along-scan deg for 1 pixel

    # adjust constants for m band
    if band == "m":
        pt *= 2
        scan_rate *= 2
        st1 /= 2
        st2 /= 2
        st3 /= 2
        sb1 /= 2
        sb2 /= 2
        sb3 /= 2
        ps_const *= 2

    # derive more constants
    st = pt / h  # along-track deg for 1 pixel[rad]
    r = earth_radius + h  # r for satellite[km]

    # calculate along-scan degrees
    abs_sample = (sample >= st3) * (sample + 1 - st3) + (sample < st3) * (st3 - sample)
    zone1 = abs_sample <= st1
    zone2 = np.logical_and(abs_sample > st1, abs_sample <= st2)
    zone3 = abs_sample > st2
    ps = ps_const * (3 * zone1 + 2 * zone2 + zone3)
    ss = ps / h  # along-scan deg for 1 pixel[rad]
    abs_scan = (
        zone1 * (abs_sample * 3 * scan_rate)
        + zone2 * ((sb2 + ((abs_sample - st1) * 2)) * scan_rate)
        + zone3 * ((sb3 + (abs_sample - st2)) * scan_rate)
    )  # scan angle, in degree
    theta = np.deg2rad(abs_scan)  # scan angle, in radian
    cos_theta = np.cos(theta)

    # calculate pixel size DS and DT
    temp = (earth_radius / r) ** 2.0 - np.sin(theta) ** 2.0
    sqrt_temp = np.sqrt(temp)
    DS = earth_radius * ss * (cos_theta / sqrt_temp - 1.0)
    DT = r * st * (cos_theta - sqrt_temp)

    if rtSCAN_ANGLE == True:
        scan_angle = np.rad2deg(theta)
        return DS, DT, scan_angle
    else:
        return DS, DT


def read_geojson_nv_CA(y0=2012, y1=2019):
    """read non vegetation fire data in California

    Returns
    -------
    gdf : GeoDataFrame
        the points of non vegetation fire location
    """
    fnm = (
        settings.dirextdata + "CA/Calnvf/FCs_nv_" + str(y0) + "-" + str(y1) + ".geojson"
    )
    gdf = gpd_read_file(fnm)

    return gdf


def VNP14IMGML_filepath(t: TimeStep):
    """Filepath for monthly S-NPP VIIRS data

    Parameters
    ----------
    t : tuple, (int,int,int,str)
        the year, month, day and 'AM'|'PM' during the initialization

    Returns
    -------
    filepath : str | None
        Path to input data or None if file does not exist
    """
    year, month = t[0], t[1]

    file_dir = os.path.join(
        settings.dirextdata,
        "VIIRS",
        "VNP14IMGML",
    )

    filepath = os.path.join(file_dir, f"VNP14IMGML.{year}{month:02}.C1.05.txt")
    if not settings.fs.exists(filepath):
        filepath = os.path.join(file_dir, f"VNP14IMGML.{year}{month:02}.C2.01.txt")
    if not settings.fs.exists(filepath):
        print("No data available for file", filepath)
        return

    return filepath


def read_VNP14IMGML(filepath: str):
    """read monthly S-NPP VIIRS data

    Parameters
    ----------
    filepath : str
        Path to input data. Can be local or s3.

    Returns
    -------
    df : pandas.DataFrame
        monthly DataFrame containing standardized columns of VIIRS active fires
    """
    usecols = [
        "YYYYMMDD",
        "HHMM",
        "Lat",
        "Lon",
        "Line",
        "Sample",
        "FRP",
        "Confidence",
        "Type",
        "DNFlag",
    ]

    df = pd.read_csv(
        filepath,
        dtype={"YYYYMMDD": "string", "HHMM": "string"},
        usecols=usecols,
        skipinitialspace=True,
        na_values={"FRP": "*******"},
    )
    df["datetime"] = pd.to_datetime(
        df["YYYYMMDD"] + " " + df["HHMM"], format="%Y%m%d %H%M"
    )
    df["DT"], df["DS"] = viirs_pixel_size(df["Sample"].values)
    df = df.drop(columns=["Sample", "Line"])
    return df


def VNP14IMGTDL_filepath(t: TimeStep):
    """Filepath for NRT S-NPP VIIRS data

    Parameters
    ----------
    t : tuple, (int,int,int,str)
        the year, month, day and 'AM'|'PM' during the initialization

    Returns
    -------
    filepath : str
        Path to input data or None if file does not exist
    """
    d = date(t[0], t[1], t[2])

    filepath = os.path.join(
        settings.dirextdata,
        "VIIRS",
        "VNP14IMGTDL",
        f"SUOMI_VIIRS_C2_Global_VNP14IMGTDL_NRT_{d.strftime('%Y%j')}.txt",
    )
    if not settings.fs.exists(filepath):
        print("No data available for file", filepath)
        return

    return filepath


def read_VNP14IMGTDL(filepath: str):
    """Read daily NRT S-NPP VIIRS fire location data

    Parameters
    ----------
    filepath : str
        Path to input data. Can be local or s3.

    Returns
    -------
    df : pandas.DataFrame
        DataFrame containing standardized columns of daily active fires
    """
    usecols = [
        "latitude",
        "longitude",
        "daynight",
        "frp",
        "confidence",
        "scan",
        "track",
        "acq_date",
        "acq_time",
    ]
    df = pd.read_csv(
        filepath,
        dtype={"acq_date": "string", "acq_time": "string"},
        usecols=usecols,
        skipinitialspace=True,
    )
    df["datetime"] = pd.to_datetime(
        df["acq_date"] + " " + df["acq_time"], format="%Y-%m-%d %H:%M"
    )
    df = df.rename(
        columns={
            "latitude": "Lat",
            "longitude": "Lon",
            "frp": "FRP",
            "scan": "DS",
            "track": "DT",
        }
    )
    return df


def VJ114IMGML_filepath(t: TimeStep):
    """Filepath for monthly NOAA20 VIIRS data

    Parameters
    ----------
    t : tuple, (int,int,int,str)
        the year, month, day and 'AM'|'PM' during the initialization

    Returns
    -------
    filepath : str
        Path to input data or None if file does not exist
    """
    filepath = os.path.join(
        settings.dirextdata,
        "VIIRS",
        "VJ114IMGML",
        str(t[0]),
        f"VJ114IMGML_{t[0]}{t[1]:02}.txt",
    )
    if not settings.fs.exists(filepath):
        print("No data available for file", filepath)
        return

    return filepath


def read_VJ114IMGML(filepath: str):
    """read monthly NOAA20 VIIRS fire location data

    Parameters
    ----------
    filepath : str
        Path to input data. Can be local or s3.

    Returns
    -------
    df : pandas.DataFrame
        monthly DataFrame containing standardized columns of VIIRS active fires
    """
    usecols = [
        "year",
        "month",
        "day",
        "hh",
        "mm",
        "lon",
        "lat",
        "mask",
        "line",
        "sample",
        "frp",
    ]

    df = pd.read_csv(
        filepath,
        dtype={col: "string" for col in ["year", "month", "day", "hh", "mm"]},
        usecols=usecols,
        skipinitialspace=True,
    )
    df["datetime"] = pd.to_datetime(
        df["year"]
        + "-"
        + df["month"]
        + "-"
        + df["day"]
        + " "
        + df["hh"]
        + ":"
        + df["mm"],
        format="%Y-%m-%d %H:%M",
    )
    df = df.rename(
        columns={
            "lat": "Lat",
            "lon": "Lon",
            "frp": "FRP",
            "line": "Line",
            "sample": "Sample",
        }
    )
    df["DT"], df["DS"] = viirs_pixel_size(df["Sample"].values)
    return df


def VJ114IMGTDL_filepath(t: TimeStep):
    """Filepath for daily NRT NOAA20 VIIRS data

    Parameters
    ----------
    t : tuple, (int,int,int,str)
        the year, month, day and 'AM'|'PM' during the initialization

    Returns
    -------
    filepath : str
        Path to input data or None if file does not exist
    """
    d = date(t[0], t[1], t[2])

    filepath = os.path.join(
        settings.dirextdata,
        "VIIRS",
        "VJ114IMGTDL",
        f"J1_VIIRS_C2_Global_VJ114IMGTDL_NRT_{d.strftime('%Y%j')}.txt",
    )
    if not settings.fs.exists(filepath):
        print("No data available for file", filepath)
        return

    return filepath


def read_VJ114IMGTDL(filepath: str):
    """Read daily NRT NOAA20 VIIRS fire location data

    NOTE: this function expects julian date to be encoded in the filename

    Parameters
    ----------
    filepath : str
        Path to input data. Can be local or s3.

    Returns
    -------
    df : pandas.DataFrame
        DataFrame containing standardized columns of daily active fires
    """
    # strip the julian date out of the filepath
    julian_date = filepath.split("_")[-1].split(".")[0]

    # convert to a datetime.date object
    d = datetime.strptime(julian_date, "%Y%j").date()

    df = pd.read_csv(filepath)
    df["acq_date"] = str(d)
    df["datetime"] = pd.to_datetime(
        df["acq_date"] + " " + df["acq_time"], format="%Y-%m-%d %H:%M"
    )
    df = df.rename(
        columns={
            "latitude": "Lat",
            "longitude": "Lon",
            "frp": "FRP",
            "scan": "DS",
            "track": "DT",
        }
    )
    return df


def AFP_regfilter(df, shp_Reg):
    """filter fire pixels using a given shp_Reg

    Parameters
    ----------
    df : pandas DataFrame
        fire pixel data with lat and lon information
    shp_Reg : geometry
        the geometry of the region shape

    Returns
    -------
    df_filtered : pandas DataFrame
        the filtered fire pixels
    """
    # preliminary spatial filter and quality filter
    regext = shp_Reg.bounds
    newfirepixels = df.loc[
        (df["Lat"] >= regext[1])
        & (df["Lat"] <= regext[3])
        & (df["Lon"] >= regext[0])
        & (df["Lon"] <= regext[2])
    ]
    point_data = [Point(xy) for xy in zip(newfirepixels["Lon"], newfirepixels["Lat"])]
    gdf_filtered = gpd.GeoDataFrame(newfirepixels, geometry=point_data, crs=4326)

    # Do detailed filtering (within shp_Reg)
    gdf_filtered = gdf_filtered[gdf_filtered["geometry"].within(shp_Reg)]

    # drop geometry column and project to epsg
    df_filtered = AFP_toprj(gdf_filtered)

    return df_filtered


def AFP_nonstatfilter(df):
    """Function used to filter non-static VIIRS active fire data

    Parameters
    ----------
    df : pandas DataFrame
        fire pixel data

    Returns
    -------
    df_filtered : pandas DataFrame
        the filtered fire pixels
    """

    # filter non-veg fires using a pre-derived mask
    gdf_nv = read_geojson_nv_CA()
    buf_nv = 0.0071  # buffer for each nv point (now set to 0.0071 deg = 0.71 km)
    df_filtered = df[
        (df.geometry.within(gdf_nv.iloc[0].geometry.buffer(0.005)) == False)
    ]

    return df_filtered


def AFP_setampm(df):
    """Function to set ampm column (using local hour) and update the df

    Parameters
    ----------
    df : pandas DataFrame
        fire pixel data, with 'Lon' and 'datetime' column

    Returns
    -------
    df_withampm : pandas DataFrame
        the DataFrame with 'ampm' column
    """
    # calculate local hour using the longitude and datetime column
    localhour = (pd.to_timedelta(df.Lon / 15, unit="hours") + df["datetime"]).dt.hour

    # set am/pm flag based on local hour
    df_withampm = df.assign(
        ampm=np.where(((localhour > 6) & (localhour < 18)), "PM", "AM")
    )

    return df_withampm


def AFP_toprj(gdf):
    """Transforms lat/lon coordinates to projected coordinate system for computations in m
    for global studies probably proj:cea would be a good choice
    for the boreals we use North Pole LAEA (epsg:3571)
    for CA may use WGS 84 / UTM zone 10N (epsg: 32610)
    for US may use US National Atlas Equal Area (epsg: 9311)"""

    gdf = gdf.to_crs(epsg=settings.EPSG_CODE)
    gdf["x"] = gdf.geometry.x
    gdf["y"] = gdf.geometry.y
    df = pd.DataFrame(gdf.drop(columns="geometry"))

    return df


# ------------------------------------------------------------------------------
# %% Read other datasets
# ------------------------------------------------------------------------------


def read_mcd64_pixels(year, ext=[]):
    """Reads all Modis burned area pixels of a year from text file

    Parameters
    ----------
    year: int
    ext: list or tuple, (minx, miny, maxx, maxy) extent for filtering pixels

    Returns
    -------
    df: pd.dataframe, dataframe containing lat, lon and day of burning
        from modis burned area"""

    mcd64dir = os.path.join(settings.dirextdata, "MCD64A1")
    # filelist_viirs = glob.glob(dirpjdata+str(year)+'/Snapshot/*_NFP.txt')
    # df = pd.concat([pd.read_csv(i, dtype = {'lat':np.float64, 'lon':np.float64},
    #                             usecols = ['lat', 'lon'])
    #                 for i in filelist_viirs], ignore_index=True)
    # df['sensor'] = 'viirs'

    filelist_mcd64 = os.path.join(mcd64dir, "ba_centroids_" + str(year) + ".csv")
    df = pd.read_csv(
        filelist_mcd64,
        dtype={"lat": np.float64, "lon": np.float64},
        usecols=["lat", "lon", "doy"],
    )
    # filter by bounding box
    if len(ext) == 4:
        minx, miny, maxx, maxy = ext
        df = df.loc[
            (df["lat"] >= miny)
            & (df["lat"] <= maxy)
            & (df["lon"] >= minx)
            & (df["lon"] <= maxx)
        ]

    # df = pd.concat([df, df2])

    return df


def load_mcd64(year, xoff=0, yoff=0, xsize=None, ysize=None):
    """get annual circumpolar modis burned/unburned tif
    optional: clip using UL pixel and pixel dimensions

    Parameters
    ----------
    year: int
    xoff, yoff: int, UL pixel coordinates for clipping
    xsize,ysize: int, pixels in x and y direction for clipping

    Returns
    -------
    arr: np.array, clipped image"""

    # TODO: remove dead codepath
    import gdal

    mcd64dir = os.path.join(settings.dirextdata, "MCD64A1")
    fnm = os.path.join(mcd64dir, "mcd64_" + str(year) + ".tif")
    ds = gdal.Open(fnm)
    # arr = ds.ReadAsArray(xsize=xsize, ysize=ysize)
    arr = ds.ReadAsArray(xoff=xoff, yoff=yoff, xsize=xsize, ysize=ysize)

    return arr


def get_any_shp(filename):
    """get shapefile of any region given the input file name

    Parameters
    ----------
    filename : str
        the shapefile names saved in the directory dirextdata/shapefiles/
    """
    # find the california shapefile
    dirshape = os.path.join(settings.dirextdata, "shapefiles")
    statefnm = os.path.join(dirshape, filename)

    # read the geometry
    shp = gpd_read_file(statefnm).iloc[0].geometry

    return shp


def get_Cal_shp():
    """get shapefile of California
    !!! this function can be realized using get_reg_shp(); still keep here for convenience...;
        will be deleted or modified later !!!
    """
    # find the california shapefile
    statefnm = os.path.join(settings.dirextdata, "CA", "Calshape", "California.shp")

    # read the geometry
    shp_Cal = gpd_read_file(statefnm).iloc[0].geometry

    return shp_Cal


def get_Cty_shp(ctr):
    """get shapefile of a country

    Parameters
    ----------
    ctr : str
        country name
    """
    ctyfnm = os.path.join(settings.dirextdata, "World", "country.shp")

    gdf_cty = gpd_read_file(ctyfnm)

    if ctr in gdf_cty["CNTRY_NAME"].values:
        g = gdf_cty[gdf_cty.CNTRY_NAME == ctr].iloc[0].geometry
        return g
    else:
        return None


def get_reg_shp(reg):
    """return the shape of a region, given an optional reg input

    Parameters
    ----------
    reg : object
        region definition, one of the following
         - a geometry
         - a four-element list showing the extent of the region [lonmin,latmin,lonmax,latmax]
         - a country name

    Returns
    -------
    shp_Reg : geometry
        the geometry of the region shape
    """

    # read or form shape used for filtering active fires
    if isinstance(reg, shapely.geometry.base.BaseGeometry):
        shp_Reg = reg
    elif isinstance(reg, str):
        shp_Reg = get_Cty_shp(reg)
        if shp_Reg is None:
            print("Please input a valid Country name")
            return None
    elif isinstance(reg, list):
        shp_Reg = Polygon(
            [
                [reg[0], reg[1]],
                [reg[0], reg[3]],
                [reg[2], reg[3]],
                [reg[2], reg[1]],
                [reg[0], reg[1]],
            ]
        )
    else:
        print(
            "Please use geometry, country name (in str), or [lonmin,latmin,lonmax,latmax] list for the parameter region"
        )
        return None

    return shp_Reg


def get_LCT(locs):
    """Get land cover type for active fires

    Parameters
    ----------
    locs : list of lists (nx2)
        lat and lon values for each active fire detection

    Returns
    -------
    vLCT : list of ints
        land cover types for all input active fires
    """
    # read NLCD 500m data
    fnmLCT = os.path.join(settings.dirextdata, "CA", "nlcd_510m.tif")
    dataset = rasterio.open(fnmLCT)
    transformer = pyproj.Transformer.from_crs("epsg:4326", dataset.crs)
    locs_crs_x, locs_crs_y = transformer.transform(
        # NOTE: EPSG 4326 expected coordinate order latitude, longitude, but
        # `locs` is x (longitude), y (latitude). That's why `l[1]`, then `l[0]`
        # here.
        [l[1] for l in locs],
        [l[0] for l in locs],
    )
    locs_crs = list(zip(locs_crs_x, locs_crs_y))
    samps = list(dataset.sample(locs_crs))
    vLCT = [int(s) for s in samps]
    return vLCT


def get_LCT_CONUS(locs):
    """Get land cover type for active fires - CONUS scale.
        This is the same function as get_LCT but with a CONUS wide file.

    Parameters
    ----------
    locs : np.array (nx2)
        lat and lon values for each active fire detection

    Returns
    -------
    vLCT : list of ints
        land cover types for all input active fires
    """
    from fireatlas.preprocess import preprocessed_landcover_filename

    # read NLCD 500m data
    fnmLCT = preprocessed_landcover_filename("nlcd_export_510m_simplified")
    dataset = rasterio.open(fnmLCT)
    vLCT = dataset.sample(locs, indexes=1)
    vLCT = [lc[0] for lc in vLCT]
    return vLCT


def get_LCT_Global(locs):
    """Get land cover type for active fires - CONUS scale.
        This is the same function as get_LCT but with global file.

    Parameters
    ----------
    locs : np.array (nx2)
        lat and lon values for each active fire detection

    Returns
    -------
    vLCT : list of ints
        land cover types for all input active fires
    """
    fnmLCT = os.path.join(settings.dirextdata, "GlobalLC", "global_lc_mosaic.tif")
    dataset = rasterio.open(fnmLCT)

    # previous LC data sources were in a different crs and needed a transform
    # the VIIRS and LC data in this case are both in EPSG:4326, so we can sample directly
    samps = list(dataset.sample(locs))
    vLCT = [s[0] for s in samps]
    return vLCT


def get_LCT_NLCD(locs):
    """Get land cover type from NCLD for multiple locations

    Parameters
    ----------
    locs : np.array (nx2)
        lon and lat values for each active fire detection

    Returns
    -------
    vLCT : list of ints
        land cover types for all input active fires
    """
    # read NLCD 500m data
    fnmLCT = os.path.join(settings.dirextdata, "CA", "nlcd_510m_latlon.tif")
    dataset = rasterio.open(fnmLCT)
    vLCT = dataset.sample(locs, indexes=1)
    vLCT = [lc[0] for lc in vLCT]  # list values

    return vLCT


def get_FM1000(t, lon, lat):
    """Get fm1000 for a point at t

    Parameters
    ----------
    t : datetime date
        date
    lon : float
        longitude value
    lat : float
        latitude value
    Returns
    -------
    FM1000_loc : list of floats
        fm1000 value for all input active fires
    """
    warnings.simplefilter("ignore")

    # read annual fm1000 data
    dirGridMET = os.path.join(settings.dirextdata, "GridMET") + "/"
    fnm = dirGridMET + "fm1000_" + t.strftime("%Y") + ".zarr"
    ds = xr.open_zarr(fnm)
    FM1000_all = ds["dead_fuel_moisture_1000hr"]

    # extract daily data at t
    try:
        FM1000_day = FM1000_all.sel(day=t.strftime("%Y-%m-%d"))
    except:  # if data are not available, use the last available date
        FM1000_day = FM1000_all.isel(day=-1)

    # extract data near the given location
    FM1000_loc = FM1000_day.sel(lon=lon, lat=lat, method="nearest").item()

    return FM1000_loc


# ------------------------------------------------------------------------------
# %% read and load object, gdf and summary related files
# ------------------------------------------------------------------------------
def check_filefolder(fnm):
    """if the folder containing a file does not exist, make it

    Parameters
    ----------
    fnm : str
        file name
    """
    # folder name
    dirnm = os.path.dirname(fnm)
    if dirnm.startswith("s3://"):
        # No concept of "folders" in S3, so no need to create them.
        return None

    # create the folder if needed
    if not os.path.exists(dirnm):
        os.makedirs(dirnm)


def check_folder(dirnm):
    """if the folder does not exist, make it

    Parameters
    ----------
    dirfnm : str
        folder name
    """
    if dirnm.startswith("s3://"):
        # No concept of "folders" in S3, so no need to create them.
        return None

    # create the folder if needed
    if not os.path.exists(dirnm):
        os.makedirs(dirnm)


def correct_dtype(gdf, op=""):
    """correct the datatype for gdfs loaded from geojson files

    Parameters
    ----------
    gdf : gpd DataFrame
        the gdf directly read from the geojson file
    op : str
        the option for the geojson data -
        '': fire perimeter and attributes
        'FL': active fire line
        'NFP': new fire pixels

    Returns
    -------
    gdf : gpd DataFrame
        the gdf with correct datatype as given in FireConsts
    """
    # TODO: remove dead codepath
    from fireatlas.FireConsts import dd

    # explicitly set the attributes data types
    if op == "":
        for v, tp in dd.items():
            gdf[v] = gdf[v].astype(tp)
    else:
        gdf["fireID"] = gdf["fireID"].astype("int")

    return gdf


def get_fobj_fnm(t, regnm, activeonly=False):
    """Return the fire object pickle file name at a time step
    Parameters
    ----------
    t : tuple, (year,month,day,str)
        the day and 'AM'|'PM' during the intialization

    Returns
    ----------
    fnm : str
        pickle file name
    """
    d = date(*t[:-1])
    # fnm = dirpjdata+regnm+'/'+d.strftime('%Y')+'/Serialization/'+d.strftime('%Y%m%d')+t[-1]+'.pkl'
    if activeonly:
        fnm = os.path.join(
            settings.diroutdata,
            regnm,
            d.strftime("%Y"),
            "Serialization",
            d.strftime("%Y%m%d") + t[-1] + ".pkl",
        )
    else:
        fnm = os.path.join(
            settings.diroutdata,
            regnm,
            d.strftime("%Y"),
            "Serialization",
            d.strftime("%Y%m%d") + t[-1] + "_full.pkl",
        )
    return fnm


def check_fobj(t, regnm, activeonly=False):
    """Check if the pickle file storing a daily allfires object exists

    Parameters
    ----------
    t : tuple, (year,month,day,str)
        the day and 'AM'|'PM' during the intialization
    """

    fnm = get_fobj_fnm(t, regnm, activeonly=activeonly)
    return settings.fs.exists(fnm)


def save_fobj(allfires, t, regnm, activeonly=False):
    """Save a daily allfires object to a pickle file

    Parameters
    ----------
    allfires : obj of Allfires class
        daily allfires
    t : tuple, (year,month,day,str)
        the day and 'AM'|'PM' during the intialization
    activeonly : bool
        the flag to save activeonly or full pickle file
    """
    from fireatlas.FireObj import Allfires

    if activeonly:
        # we only want to save active and sleeper fires
        allfires_out = Allfires(t)  # create a new allfires object
        # allfires_out.update_t(t)  # correct the time (previous time step was used in the intialization)
        fids_out = (
            allfires.fids_active + allfires.fids_sleeper
        )  # active fire or sleeper
        id_dict = []
        allfires_index = 0
        for fid in fids_out:
            # allfires_out.fires.append(allfires.fires[fid])
            allfires_out.fires[fid] = allfires.fires[fid]
            id_dict.append((allfires_index, fid))
            allfires_index += 1

        # copy the heritage from the original allfires object
        allfires_out.heritages = allfires.heritages
        allfires_out.id_dict = id_dict
    else:
        allfires_out = allfires

    # get output file name
    fnm = get_fobj_fnm(t, regnm, activeonly=activeonly)

    # check folder
    check_filefolder(fnm)

    # save
    with fsspec.open(fnm, "wb") as f:
        pickle.dump(allfires_out, f)


def load_fobj(t, regnm, activeonly=False):
    """Load a daily allfires object from a pickle file

    Parameters
    ----------
    t : tuple, (year,month,day,str)
        the day and 'AM'|'PM' during the intialization

    Returns
    ----------
    data : obj of Allfires class
        daily allfires
    """
    # get file name
    fnm = get_fobj_fnm(t, regnm, activeonly=activeonly)

    # load data
    with fsspec.open(fnm, "rb") as f:
        data = pickle.load(f)
    return data


# def update_fobj_ia(data,t,regnm):
#     ''' Update the pickle file saving inactive fire objects (one for each year)
#
#     Parameters
#     ----------
#     data : obj of Allfires class
#         daily inactive allfires for updating
#     t : tuple, (year,month,day,str)
#         the day and 'AM'|'PM' during the intialization
#     '''
#     import pickle
#     from FireConsts import diroutdata
#     from datetime import date
#     import os
#
#     # read and update
#     fnm = os.path.join(diroutdata, regnm, str(t[0]),'Serialization', 'Inactive_fobj.pkl')
#     if os.path.exists(fnm):  # if file exist, read and update data
#         with open(fnm,'rb') as f:
#             allfires_ia = pickle.load(f)
#             allfires_ia.fires =  allfires_ia.fires + data.fires
#     else:   # if not exist, use the input data
#         # check/create folder
#         check_filefolder(fnm)
#         allfires_ia = data
#
#     # save data back to the file
#     with open(fnm,'wb') as f:
#         pickle.dump(allfires_ia, f)


def get_gdfobj_fnm(t, regnm, op=""):
    """Return the fire object gpkg file name at a time step
    Parameters
    ----------
    t : tuple, (year,month,day,str)
        the day and 'AM'|'PM' during the intialization
    op : str
        the option for the geojson data -
        '': fire perimeter and attributes
        'FL': active fire line
        'NFP': new fire pixels
    Returns
    ----------
    fnm : str
        gdf file name
    """
    d = date(*t[:-1])
    if op == "":
        fnm = os.path.join(
            settings.diroutdata,
            regnm,
            d.strftime("%Y"),
            "Snapshot",
            d.strftime("%Y%m%d") + t[-1] + ".gpkg",
        )
    else:
        fnm = os.path.join(
            settings.diroutdata,
            regnm,
            d.strftime("%Y"),
            "Snapshot",
            d.strftime("%Y%m%d") + t[-1] + "_" + op + ".gpkg",
        )

    return fnm


def get_gpkgobj_fnm(t, regnm):
    """Return the gpkg snapshot file name at a time step
    Parameters
    ----------
    t : tuple, (year,month,day,str)
        the day and 'AM'|'PM' during the intialization
    Returns
    ----------
    fnm : str
        gpkg file name
    """
    # determine output dir
    d = date(*t[:-1])
    strdir = os.path.join(settings.diroutdata, regnm, d.strftime("%Y"), "Snapshot")

    # get the output file name
    fnm = os.path.join(strdir, d.strftime("%Y%m%d") + t[-1])

    return fnm


def get_gpkgsfs_fnm(t, fid, regnm):
    """Return the gpkg snapshot file name at a time step
    Parameters
    ----------
    t : tuple, (year,month,day,str)
        the day and 'AM'|'PM' during the intialization
    Returns
    ----------
    fnm : str
        gpkg file name
    """
    # determine output dir
    d = date(*t[:-1])
    strdir = os.path.join(settings.diroutdata, regnm, d.strftime("%Y"), "Largefire")

    # get the output file name
    fnm = os.path.join(strdir, "F" + str(int(fid)) + "_" + d.strftime("%Y%m%d") + t[-1])

    return fnm


def get_gpkgsfs_dir(yr, regnm):
    """Return the gpkg snapshot directory name for a year
    Parameters
    ----------
    yr : int
        year
    regnm : str
        region name
    Returns
    ----------
    fnm : str
        gpkg file name
    """
    # determine output dir
    strdir = os.path.join(settings.diroutdata, regnm, str(yr), "Largefire")

    return strdir


def get_NFPlistsfs_fnm(t, fid, regnm):
    """Return the gpkg snapshot file name at a time step
    Parameters
    ----------
    t : tuple, (year,month,day,str)
        the day and 'AM'|'PM' during the intialization
    Returns
    ----------
    fnm : str
        gpkg file name
    """
    # determine output dir
    d = date(*t[:-1])
    strdir = os.path.join(settings.diroutdata, regnm, d.strftime("%Y"), "Largefire")

    # get the output file name
    fnm = os.path.join(
        strdir, "F" + str(int(fid)) + "_" + d.strftime("%Y%m%d") + t[-1] + "_NFP.txt"
    )

    return fnm


def check_gpkgobj(t, regnm):
    """Check if the gpkg file storing a daily allfires attributes exists

    Parameters
    ----------
    t : tuple, (year,month,day,str)
        the day and 'AM'|'PM' during the intialization
    regnm : str
        the name of the region
    """
    # d = date(*t[:-1])
    fnm = get_gpkgobj_fnm(t, regnm)

    return settings.fs.exists(fnm)


def check_gdfobj(t, regnm, op=""):
    """Check if the gpkg file storing a daily allfires attributes exists

    Parameters
    ----------
    t : tuple, (year,month,day,str)
        the day and 'AM'|'PM' during the intialization
    op : str
        the option for the geojson data -
        '': fire perimeter and attributes
        'FL': active fire line
        'NFP': new fire pixels
    """
    d = date(*t[:-1])
    fnm = get_gdfobj_fnm(t, regnm)

    return settings.fs.exists(fnm)


def save_gdfobj(gdf, t, regnm, param="", fid="", op=""):
    """Save geopandas to a gpgk file

    Parameters
    ----------
    gdf : geopandas DataFrame
        daily allfires diagnostic parameters
    t : tuple, (int,int,int,str)
        the year, month, day and 'AM'|'PM'
    param : string
        if empty: save daily allfires diagnostic dataframe
        if 'large': save daily single fire
        if 'ign': save ignition layer
        if 'final': save final perimeters layer
    fid : int,
        fire id (needed only with param = 'large')
    """
    # get file name
    if param == "":
        fnm = get_gdfobj_fnm(t, regnm, op=op)
    elif param == "large":
        fnm = get_gdfobj_sf_fnm(t, fid, regnm, op=op)
    else:

        d = date(*t[:-1])

        # get file name
        fnm = os.path.join(
            settings.diroutdata,
            regnm,
            d.strftime("%Y"),
            "Summary",
            param + d.strftime("%Y") + ".gpkg",
        )

    # check folder
    check_filefolder(fnm)

    # create a new id column used for indexing
    # (id column will be dropped when reading gpkg)
    gdf["id"] = gdf.index
    gdf = gdf.set_index("id")
    if op == "FL":
        gdf["fireID"] = gdf["fireID"].astype(int)  # data types are screwed up in fline

    # save file
    gdf.to_file(fnm, driver="GPKG")


def save_gpkgobj(
    t, regnm, gdf_fperim=None, gdf_fline=None, gdf_nfp=None, gdf_uptonow=None
):
    """Save geopandas to a gpkg fire object file

    Parameters
    ----------
    gdf : geopandas DataFrame
        daily allfires diagnostic parameters
    t : tuple, (int,int,int,str)
        the year, month, day and 'AM'|'PM'
    """
    # date of the time step
    d = date(*t[:-1])

    # get file name
    fnm = get_gpkgobj_fnm(t, regnm)

    # check folder
    check_filefolder(fnm)

    # # create a new id column used for indexing
    # #(id column will be dropped when reading gpkg)
    # gdf["id"] = gdf.index
    # gdf = gdf.set_index('id')
    # if op == 'FL':
    #     gdf['fid'] = gdf['fid'].astype(int) # data types are screwed up in fline

    # save file
    if gdf_fperim is not None:
        gdf_fperim.to_file(f"{fnm}/perimeter.fgb", driver="FlatGeobuf")
        if settings.export_to_veda:
            copy_from_local_to_veda_s3(f"{fnm}/perimeter.fgb", regnm)

    if gdf_fline is not None:
        gdf_fline.to_file(f"{fnm}/fireline.fgb", driver="FlatGeobuf")
        if settings.export_to_veda:
            copy_from_local_to_veda_s3(f"{fnm}/fireline.fgb", regnm)

    if gdf_nfp is not None:
        gdf_nfp.to_file(f"{fnm}/newfirepix.fgb", driver="FlatGeobuf")
        if settings.export_to_veda:
            copy_from_local_to_veda_s3(f"{fnm}/newfirepix.fgb", regnm)

    if gdf_uptonow is not None:
        gdf_uptonow.to_file(f"{fnm}/uptonow.fgb", driver="FlatGeobuf")


def save_gpkgsfs(
    t, fid, regnm, gdf_fperim=None, gdf_fline=None, gdf_nfp=None, gdf_nfplist=None
):
    """Save geopandas to a gpkg fire object file

    Parameters
    ----------
    gdf : geopandas DataFrame
        daily allfires diagnostic parameters
    t : tuple, (int,int,int,str)
        the year, month, day and 'AM'|'PM'
    """
    # date of the time step
    d = date(*t[:-1])

    # get file name
    fnm = get_gpkgsfs_fnm(t, fid, regnm)

    # check folder
    check_filefolder(fnm)

    # # create a new id column used for indexing
    # #(id column will be dropped when reading gpkg)
    # gdf["id"] = gdf.index
    # gdf = gdf.set_index('id')
    # if op == 'FL':
    #     gdf['fid'] = gdf['fid'].astype(int) # data types are screwed up in fline

    # save file
    if gdf_fperim is not None:
        gdf_fperim.to_file(f"{fnm}/perimeter.fgb", driver="FlatGeobuf")

    # if len(gdf_fline) > 0:
    if gdf_fline is not None:
        gdf_fline.to_file(f"{fnm}/fireline.fgb", driver="FlatGeobuf")

    # if len(gdf_NFP) > 0:
    if gdf_nfp is not None:
        gdf_nfp.to_file(f"{fnm}/newfirepix.fgb", driver="FlatGeobuf")

    if gdf_nfplist is not None:
        gdf_nfplist.to_file(f"{fnm}/nfplist.fgb", driver="FlatGeobuf")


def load_gpkgobj(t, regnm, layer="perimeter"):
    """Load geopandas from a gpkg fire object file
    Parameters
    ----------
    t : tuple, (int,int,int,str)
        the year, month, day and 'AM'|'PM'
    layer : str
        the layer name, 'perimeter'|'fireline'|'newfirepix'
    Returns
    -------
    gdf : geopandas DataFrame
        daily allfires diagnostic parameters
    """
    # get file name
    fnm = get_gpkgobj_fnm(t, regnm)
    try:
        gdf = read_gpkg(fnm, layer=layer)
        if gdf is None:
            return gdf
        # set fireID as index column (force to be int type)
        if "t" in gdf.columns:
            gdf["t"] = pd.to_datetime(gdf.t)
        if "t_st" in gdf.columns:
            gdf["t_st"] = pd.to_datetime(gdf.t_st)
        if "t_ed" in gdf.columns:
            gdf["t_ed"] = pd.to_datetime(gdf.t_ed)

        gdf.fireID = gdf.fireID.astype("int")
        gdf = gdf.set_index("fireID")
        return gdf

    except Exception as e:
        print("Encountered the following error:", e)
        gdf = None
        return gdf


def load_gpkgsfs(t, fid, regnm, layer="perimeter"):
    """Load geopandas from a gpkg fire object file

    Parameters
    ----------
    t : tuple, (int,int,int,str)
        the year, month, day and 'AM'|'PM'
    layer : str
        the layer name, 'perimeter'|'fireline'|'newfirepix'
    Returns
    -------
    gdf : geopandas DataFrame
        daily allfires diagnostic parameters
    """
    # get file name
    fnm = get_gpkgsfs_fnm(t, fid, regnm)

    try:
        gdf = read_gpkg(fnm, layer=layer)

        if gdf is not None:
            if "t" in gdf.columns:
                gdf["t"] = gdf["t"].astype("datetime64")
        # set fireID as index column (force to be int type)
        # gdf = gdf.set_index('index')

        return gdf
    except:
        return None


def load_gdfobj(regnm, t="", op=""):
    """Load daily allfires diagnostic dataframe as geopandas gdf

    Parameters
    ----------
    t : tuple, (int,int,int,str)
        the year, month, day and 'AM'|'PM'
    Returns
    ----------
    gdf : geopandas DataFrame
        daily allfires diagnostic parameters
    """
    # get file name
    fnm = get_gdfobj_fnm(t, regnm, op=op)

    # read data as gpd DataFrame
    gdf = gpd_read_file(fnm)

    # with gpkg the id column has to be renamed back to fid
    gdf = gdf.rename(columns={"id": "fid"})

    # correct the datatype
    gdf = correct_dtype(gdf, op=op)

    # set fid as the index
    gdf = gdf.set_index("fireID")

    return gdf


def save_newyearfidmapping(fidmapping, year, regnm):
    """Save the cross year fid mapping tuples"""
    # convert list to dataframe
    df = pd.DataFrame(fidmapping, columns=["oldfid", "newfid"])

    # determine output file name
    strdir = os.path.join(settings.diroutdata, regnm, str(year), "Summary")
    fnmout = os.path.join(strdir, "CrossyrFidmapping_" + str(year) + ".csv")
    check_filefolder(fnmout)

    # save
    df.to_csv(fnmout)  # use df = pd.read_csv(fnmout,index_col=0) to read


def save_FP_txt(df, t, regnm):

    # get filename of new fire pixels product
    fnm = get_gdfobj_fnm(t, regnm, op="NFP")
    fnm = fnm[:-4] + "txt"  # change ending to txt

    # check folder
    check_filefolder(fnm)

    # write out
    if len(df) > 0:
        df.to_csv(fnm)


def save_FPsfs_txt(df, t, fid, regnm):

    # get filename of new fire pixels product
    fnm = get_NFPlistsfs_fnm(t, fid, regnm)

    # check folder
    check_filefolder(fnm)

    # write out
    if len(df) > 0:
        df.to_csv(fnm)


def load_FP_txt(t, regnm):
    # get filename of new fire pixels product
    fnm = get_gdfobj_fnm(t, regnm, op="NFP")
    fnm = fnm[:-4] + "txt"  # change ending to txt

    if settings.fs.exists(fnm):
        df = pd.read_csv(fnm, parse_dates=["datetime"], index_col=0)
        return df


def load_lake_geoms(t, fid, regnm):
    """Load final perimeters as geopandas gdf

    Parameters
    ----------
    t: time tuple,
        needed for the year
    fid : int,
        the fire fid
    Returns
    ----------
    geoms : tuple (geom, geom)
        geometry of all lake perimeters within the fire
    """
    d = date(*t[:-1])

    # get file name
    fnm_lakes = os.path.join(
        settings.diroutdata,
        regnm,
        d.strftime("%Y"),
        "Summary",
        "lakes" + d.strftime("%Y") + ".gpkg",
    )

    # read data and extract target geometry
    gdf = gpd_read_file(fnm_lakes)
    gdf = gdf.set_index("mergid")
    if fid in gdf.index:
        gdf = gdf.loc[fid]
        geom_lakes = gdf["geometry"]
    else:
        geom_lakes = None

    return geom_lakes


def get_gdfobj_sf_fnm(t, fid, regnm, op=""):
    """Return the single fire fire object pickle file name at a time step
    Parameters
    ----------
    t : tuple, (year,month,day,str)
        the day and 'AM'|'PM' during the intialization
    fid : int
        fire id
    Returns
    ----------
    fnm : str
        gdf file name
    """
    d = date(*t[:-1])

    if op == "":
        fnm = os.path.join(
            settings.diroutdata,
            regnm,
            d.strftime("%Y"),
            "Largefire",
            "F" + str(int(fid)) + "_" + d.strftime("%Y%m%d") + t[-1] + ".gpkg",
        )
    else:
        fnm = os.path.join(
            settings.diroutdata,
            regnm,
            d.strftime("%Y"),
            "Largefire",
            "F"
            + str(int(fid))
            + "_"
            + d.strftime("%Y%m%d")
            + t[-1]
            + "_"
            + op
            + ".gpkg",
        )

    return fnm


def get_gdfobj_sf_fnms_year(year, fid, regnm, op=""):
    """Return the single fire fire object pickle file name at a time step
    Parameters
    ----------
    year : int
        year
    fid : int
        fire id
    Returns
    ----------
    fnms : list of str
        geojson file names
    """
    if op == "":
        fnms = settings.fs.glob(
            os.path.join(
                settings.diroutdata,
                regnm,
                str(year),
                "Largefire",
                "F" + str(int(fid)) + "_??????????.gpkg",
            )
        )
    else:
        fnms = settings.fs.glob(
            os.path.join(
                settings.diroutdata,
                regnm,
                str(year),
                "Largefire",
                "F" + str(int(fid)) + "_??????????_" + op + ".gpkg",
            )
        )

    return fnms


def save_gdfobj_sf(gdf, t, fid, regnm, op=""):
    """Save daily single fire allfires diagnostic dataframe to a geojson file

    Parameters
    ----------
    gdf : geopandas DataFrame
        daily allfires diagnostic parameters
    t : tuple, (int,int,int,str)
        the year, month, day and 'AM'|'PM'
    fid : int
        fire id
    """
    # get file name
    fnm = get_gdfobj_sf_fnm(t, fid, regnm, op=op)

    # check folder
    check_filefolder(fnm)

    # save data to file
    gdf.to_file(fnm, driver="GeoJSON")


def load_gdfobj_sf(t, fid, regnm, op=""):
    """Load single fire daily allfires diagnostic dataframe to a geojson file

    Parameters
    ----------
    t : tuple, (int,int,int,str)
        the year, month, day and 'AM'|'PM'
    fid : int
        fire id
    Returns
    ----------
    gdf : geopandas DataFrame
        time series of daily single fire diagnostic parameters
    """
    # get file name
    fnm = get_gdfobj_sf_fnm(t, fid, regnm, op=op)

    # read data as gpd DataFrame
    gdf = gpd_read_file(fnm)

    # parse time step to Datetime and set as index
    gdf["index"] = pd.to_datetime(gdf["index"])
    gdf = gdf.set_index("index")

    return gdf


def get_summary_fnm(t, regnm):
    """Return the fire summary file name at year end
    Parameters
    ----------
    t : tuple, (year,month,day,str)
        the day and 'AM'|'PM' during the intialization
    Returns
    ----------
    fnm : str
        summary netcdf file name
    """
    d = date(*t[:-1])
    fnm = os.path.join(
        settings.diroutdata,
        regnm,
        d.strftime("%Y"),
        "Summary",
        "fsummary_" + d.strftime("%Y%m%d") + t[-1] + ".nc",
    )

    # check folder
    check_filefolder(fnm)

    return fnm


def get_summary_fnm_lt(t, regnm):
    """Return the latest time step before current time when summary file exists
    Parameters
    ----------
    t : tuple, (year,month,day,str)
        the day and 'AM'|'PM' during the intialization

    Returns
    ----------
    pt : tuple, (year, month, day, pmpm)
        the lastest time step with summary file
    """
    # if there's no summary file for this year, return the first time step of the year
    fnms = settings.fs.glob(
        os.path.join(settings.diroutdata, regnm, str(t[0]), "Summary", "fsummary_*.nc")
    )
    if len(fnms) == 0:
        return None

    # if there is, find the nearest one
    endloop = False
    pt = FireTime.t_nb(t, nb="previous")
    while endloop == False:
        if settings.fs.exists(get_summary_fnm(pt, regnm)):
            return pt
        else:
            pt = FireTime.t_nb(pt, nb="previous")
            if pt[0] != t[0]:
                return None


def check_summary(t, regnm):
    """Check if the summary file storing a daily allfires object exists

    Parameters
    ----------
    t : tuple, (year,month,day,str)
        the day and 'AM'|'PM' during the intialization
    """
    # get file name
    fnm = get_summary_fnm(t, regnm)

    # return if it's present
    return settings.fs.exists(fnm)


def save_summary(ds, t, regnm):
    """Save summary info as of t a netcdf file

    Parameters
    ----------
    ds : xarray dataset
        year end summary dataset
    t : tuple, (int,int,int,str)
        the year, month, day and 'AM'|'PM'
    """
    # get file name
    fnm = get_summary_fnm(t, regnm)

    # check folder
    check_filefolder(fnm)

    # save netcdf file
    ds.to_netcdf(fnm)


def load_summary(t, regnm):
    """Load summary info from a netcdf file at t

    Parameters
    ----------
    t : tuple, (int,int,int,str)
        the year, month, day and 'AM'|'PM'

    Returns
    -------
    ds : xarray dataset
        summary dataset
    """
    # get file name
    fnm = get_summary_fnm(t, regnm)

    # read data as xarray Dataset
    ds = xr.open_dataset(fnm)

    return ds


def save_summarycsv(df, year, regnm, op="heritage"):
    """save summary csv files

    Parameters
    ----------
    df : pandas DataFrame
        the data
    year : int
        year
    op : str
        option, 'heritage'|'large'
    """
    fnm = os.path.join(
        settings.diroutdata,
        regnm,
        str(year),
        "Summary",
        "Flist_" + op + "_" + str(year) + ".csv",
    )
    check_filefolder(fnm)

    df.to_csv(fnm)


def read_summarycsv(year, regnm, op="heritage"):
    """read summary csv files

    Parameters
    ----------
    year : int
        year
    op : str
        option, 'heritage'|'large'

    Returns
    -------
    df : pandas DataFrame
        the data
    """
    fnm = os.path.join(
        settings.diroutdata,
        regnm,
        str(year),
        "Summary",
        "Flist_" + op + "_" + str(year) + ".csv",
    )
    df = pd.read_csv(fnm, index_col=0)
    return df


def get_lts_VNP14IMGTDL(year=None):
    if year == None:
        year = date.today().year

    dirFC = os.path.join(settings.dirextdata, "VNP14IMGTDL") + "/"
    fnms = settings.fs.glob(
        os.path.join(
            dirFC, "SUOMI_VIIRS_C2_Global_VNP14IMGTDL_NRT_" + str(year) + "*.txt"
        )
    )
    fnms.sort()
    DOY_lts = int(os.path.splitext(os.path.basename(fnms[-1]))[0][-3:])

    return DOY_lts


def get_lts_serialization(regnm, year=None):
    """get the time with lastest pkl data"""
    if year == None:
        year = date.today().year

    fnms = settings.fs.glob(
        os.path.join(settings.diroutdata, regnm, str(year), "Serialization", "*.pkl")
    )

    if len(fnms) > 0:
        fnms.sort()
        fnm_lts = os.path.basename(fnms[-1])

        lts = [int(fnm_lts[0:4]), int(fnm_lts[4:6]), int(fnm_lts[6:8]), fnm_lts[8:10]]
    else:
        lts = None

    return lts


def read_gpkg(fnm, layer="perimeter"):
    """read gpkg data
    op: 'perimeter', 'fireline', 'newfirepix'
    """
    try:
        gdf = gpd_read_file(fnm, layer=layer)
        return gdf
    except Exception as e:
        print(f"Failed to read GPKG file {fnm} with error: {str(e)}")
        return None


# %% other functions related to read/write


def save2gtif(arr, outfile, cols, rows, geotrans, proj):
    """write out a geotiff"""
    # TODO: remove dead codepath
    import gdal

    driver = gdal.GetDriverByName("GTiff")
    ds = driver.Create(outfile, cols, rows, 1, gdal.GDT_Byte)
    ds.SetGeoTransform(geotrans)
    ds.SetProjection(proj)
    band = ds.GetRasterBand(1)
    band.WriteArray(arr)


def geo_to_polar(lon_arr, lat_arr):
    """transform lists of geographic lat lon coordinates to polar LAEA grid (projected)"""

    proj4str = "epsg:3571"
    p_modis_grid = pyproj.Proj(proj4str)

    x_arr, y_arr = p_modis_grid(lon_arr, lat_arr)
    x_arr = np.array(x_arr)
    y_arr = np.array(y_arr)

    return x_arr, y_arr


def polar_to_geo(x_arr, y_arr):
    """transform lists of geographic lat lon coordinates to polar LAEA grid (projected)"""
    proj4str = "epsg:3571"
    p_modis_grid = pyproj.Proj(proj4str)

    lon_arr, lat_arr = p_modis_grid(x_arr, y_arr, inverse=True)
    lon_arr = np.array(lon_arr)
    lat_arr = np.array(lat_arr)

    return lon_arr, lat_arr


def world2Pixel(gt, Xgeo, Ygeo):
    """Uses a geomatrix (gdal.GetGeoTransform()) to calculate the pixel
    location of a geospatial coordinate"""

    gt = list(gt)
    Xpx = np.rint((Xgeo - gt[0]) / gt[1]).astype(int)
    Ypx = np.rint((Ygeo - gt[3]) / gt[5]).astype(int)

    return (Xpx, Ypx)


def pixel2World(gt, Xpixel, Ypixel):
    """Uses a gdal geomatrix (gdal.GetGeoTransform()) to calculate the
    geospatial coordinate of a pixel location"""

    Xgeo = gt[0] + Xpixel * gt[1] + Ypixel * gt[2]
    Ygeo = gt[3] + Xpixel * gt[4] + Ypixel * gt[5]

    return (Xgeo, Ygeo)


def copy_from_local_to_veda_s3(local_filepath: str, regnm: str, fs: s3fs.S3FileSystem | None = None):
    from_maap_s3_path = local_filepath.replace(settings.LOCAL_PATH, settings.S3_PATH)
    if fs is None:
        fs = s3fs.S3FileSystem()
    filename = os.path.basename(from_maap_s3_path)
    filename_no_ext = os.path.splitext(filename)[0]

    if "fireline.fgb" == filename:
        select_cols = ["fireID", "mergeid", "t", "primarykey", "region", "geometry"]

    elif "newfirepix.fgb" == filename:
        select_cols = ["fireID", "mergeid", "t", "primarykey", "region", "geometry"]

    elif "perimeter.fgb" == filename:
        select_cols = [
            "fireID",
            "n_pixels",
            "n_newpixels",
            "farea",
            "fperim",
            "flinelen",
            "duration",
            "pixden",
            "meanFRP",
            "isactive",
            "t",
            "primarykey",
            "region",
            "geom_counts",
            "low_confidence_grouping",
            "geometry",
        ]

    # for combine_largefire
    elif "lf_permimeter" in from_maap_s3_path:
        select_cols = [
            "n_pixels",
            "n_newpixels",
            "farea",
            "fperim",
            "flinelen",
            "duration",
            "pixden",
            "meanFRP",
            "t",
            "fireID",
            "primarykey",
            "region",
            "geom_counts",
            "low_confidence_grouping",
            "geometry",
        ]
    elif "lf_fireline" in from_maap_s3_path:
        select_cols = ["fireID", "t", "primarykey", "region", "geometry"]
    elif "lf_newfirepix" in from_maap_s3_path:
        select_cols = ["fireID", "t", "primarykey", "region", "geometry"]

    else:
        return

    new_region_name = regnm.lower().replace("_", "")

    if "lf_" in from_maap_s3_path:
        new_key_layer_name = f"{filename_no_ext}_{new_region_name}.gpkg"
        local_tmp_filepath = f"/tmp/{new_key_layer_name}"
        to_veda_s3_path = f"EIS/FEDSoutput/LFArchive/{new_key_layer_name}"
    else:
        new_key_layer_name = f"snapshot_{filename_no_ext}_nrt_{new_region_name}.gpkg"
        local_tmp_filepath = f"/tmp/{new_key_layer_name}"
        to_veda_s3_path = f"EIS/FEDSoutput/Snapshot/{new_key_layer_name}"

    gdf = gpd.read_file(local_filepath)

    # rename columns for snapshots
    # but get rid of any possible duplicate columns first
    try:
        gdf = gdf.drop(columns=["t"])
    except:
        pass
    # all LF and Snapshot outpus should have `t_ed` at this point
    gdf = gdf.rename(columns={"t_ed": "t"})

    # fiona has a bug where it cannot write GPKG files to s3 even though FileGeobuf work fine
    # so to work around this issue we just write them locally to /tmp first
    gdf[select_cols].to_file(local_tmp_filepath, driver="GPKG")
    fs.put_file(local_tmp_filepath, f"s3://veda-data-store-staging/{to_veda_s3_path}")


def copy_from_local_to_s3(filepath: str, fs: s3fs.S3FileSystem, **tags):
    """Copy from local to s3 adding any specified tags

    Some default tags will be added from the environment and specified FireConsts
    """
    dst = filepath.replace(settings.LOCAL_PATH, settings.S3_PATH)

    fs.put_file(filepath, dst)

    tags = {}
    settings_to_include_in_tags = [
        "EPSG_CODE",
        "remove_static_sources",
        "FTYP_OPT",
    ]

    # TODO: wait until JPL changes bucket policies for DPS workers to allow this
    # default_tags = {
    #     "processedBy": os.environ.get("CHE_WORKSPACE_NAMESPACE", None) or os.environ.get("JUPYTERHUB_USER", None),
    #     **settings.model_dump(include=settings_to_include_in_tags),
    # }
    # tags = {k: v for k, v in {**default_tags, **tags}.items() if v is not None}
    #
    # # there are some limitations on tags in s3:
    # #   - objects can only have 10 tags
    # #   - tag keys have a utf-16 limit of 128
    # #   - tag values have a utf-16 limit of 256
    # tags = dict([(str(k)[:64], str(v)[:128]) for k, v in tags.items()][-10:])
    #
    # s3.put_tags(dst, tags)

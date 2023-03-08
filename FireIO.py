""" FireIO
This module include functions used to read and save data
"""

# ------------------------------------------------------------------------------
#%% Read and filter active fire data
# ------------------------------------------------------------------------------

# Try to read a Geopandas file several times. Sometimes, the read fails the
# first time for mysterious reasons.
import boto3
import re
import time
from FireLog import logger


def gpd_read_file(filename, parquet=False, **kwargs):
    import geopandas as gpd
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

def os_path_exists(filename):
    """Alternative to os.path.exists that also works with S3 paths."""
    if filename.startswith("s3://"):
        import s3fs
        s3 = s3fs.S3FileSystem(anon=False)
        return s3.exists(filename)
    else:
        import os
        return os.path.exists(filename)

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
    import numpy as np

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
    """ read non vegetation fire data in California

    Returns
    -------
    gdf : GeoDataFrame
        the points of non vegetation fire location
    """
    import geopandas as gpd
    from FireConsts import dirextdata

    fnm = dirextdata + "CA/Calnvf/FCs_nv_" + str(y0) + "-" + str(y1) + ".geojson"
    gdf = gpd_read_file(fnm)

    return gdf


def read_VNP14IMGML(yr, mo, ver="C1.05"):
    """ read monthly VNP14IMGML data

    Parameters
    ----------
    yr : int
        year
    mo : int
        month

    Returns
    -------
    df : dataframe
        the monthly dataframe of VIIRS active fires
    """
    from FireConsts import dirextdata
    import os
    import pandas as pd

    # set monthly file name
    fnmFC = os.path.join(
        dirextdata,
        "VIIRS","VNP14IMGML",
        "VNP14IMGML." + str(yr) + str(mo).zfill(2) + "." + ver + ".txt",)
    
    # read
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

    
    if os_path_exists(fnmFC):
        df = pd.read_csv(
            fnmFC,
            parse_dates=[["YYYYMMDD", "HHMM"]],
            usecols=usecols,
            skipinitialspace=True,
        )
        
         # sometimes the FRP value is '*******' and cause incorrect dtype, need to correct this
        df = df.replace('*******',0)
        df['FRP'] = df['FRP'].astype('float')
        
        df["DT"], df["DS"] = viirs_pixel_size(df["Sample"].values)
        df = df.drop(columns=["Sample", "Line"])
        return df
    else:
        print('No data available for file',fnmFC)
        return None


def read_VNP14IMGTDL(t):
    """ Read daily NRT S-NPP VIIRS fire location data

    Parameters
    ----------
    t : tuple, (int,int,int,str)
        the year, month, day and 'AM'|'PM' during the intialization

    Returns
    -------
    vlist : list
        (lat,lon) tuple of all daily active fires
    """
    from FireConsts import dirextdata

    import pandas as pd
    import os
    from datetime import date

    d = date(*t[:-1])

    # derive monthly file name with path
    dirFC = os.path.join(dirextdata,'VIIRS', "VNP14IMGTDL") + "/"
    fnmFC = os.path.join(
        dirFC, "SUOMI_VIIRS_C2_Global_VNP14IMGTDL_NRT_" + d.strftime("%Y%j") + ".txt"
    )

    # read and extract
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
    if os_path_exists(fnmFC):
        # read data
        df = pd.read_csv(
            fnmFC,
            parse_dates=[["acq_date", "acq_time"]],
            usecols=usecols,
            skipinitialspace=True,
        )
        df = df.rename(
            columns={
                "latitude": "Lat",
                "longitude": "Lon",
                "frp": "FRP",
                "scan": "DS",
                "track": "DT",
                "acq_date_acq_time": "YYYYMMDD_HHMM",
            }
        )
        return df
    else:
        return None


def read_VJ114IMGML(yr, mo):
    """ read monthly VNP14IMGML data

    Parameters
    ----------
    yr : int
        year
    mo : int
        month

    Returns
    -------
    df : dataframe
        the monthly dataframe of VIIRS active fires
    """
    from FireConsts import dirextdata
    import os
    import pandas as pd

    # set monthly file name
    fnmFC = os.path.join(
        dirextdata,
        "VJ114IMGML",
        str(yr),
        "VJ114IMGML_" + str(yr) + str(mo).zfill(2) + ".txt",
    )

    # read
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

    def parser(yr, mo, dy, h, m):
        return pd.to_datetime(
            yr + "-" + mo + "-" + dy + " " + h + ":" + m, format="%Y-%m-%d %H:%M"
        )

    if os_path_exists(fnmFC):
        # df = pd.read_csv(fnmFC,usecols=usecols,skipinitialspace=True)
        df = pd.read_csv(
            fnmFC,
            parse_dates={"YYYYMMDD_HHMM": ["year", "month", "day", "hh", "mm"]},
            date_parser=parser,
            usecols=usecols,
            skipinitialspace=True,
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
    else:
        return None


def read_VJ114IMGTDL(t):
    """ Read daily NRT NOAA20 VIIRS fire location data

    Parameters
    ----------
    t : tuple, (int,int,int,str)
        the year, month, day and 'AM'|'PM' during the intialization

    Returns
    -------
    vlist : list
        (lat,lon) tuple of all daily active fires
    """
    from FireConsts import dirextdata

    import pandas as pd
    import os
    from datetime import date

    d = date(*t[:-1])
    # derive monthly file name with path
    dirFC = os.path.join(dirextdata,'VIIRS', "VJ114IMGTDL") + "/"
    fnmFC = os.path.join(
        dirFC, "J1_VIIRS_C2_Global_VJ114IMGTDL_NRT_" + d.strftime("%Y%j") + ".txt"
    )

    # read and extract
    # usecols = ['latitude','longitude','daynight','frp','confidence']
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
    if os_path_exists(fnmFC):
        # read data
        #df = pd.read_csv(
        #    fnmFC,
        #    parse_dates=[["acq_date", "acq_time"]],
        #    usecols=usecols,
        #    skipinitialspace=True,
        #)
        df = pd.read_csv(fnmFC)
        df['acq_date'] = str(d)
        df["acq_date_acq_time"] = pd.to_datetime(df['acq_date'] +' '+ df['acq_time'])
        #print(df['acq_date_acq_time'])
        df = df.rename(
            columns={
                "latitude": "Lat",
                "longitude": "Lon",
                "frp": "FRP",
                "scan": "DS",
                "track": "DT",
                "acq_date_acq_time": "YYYYMMDD_HHMM",
            }
        )
        return df
    else:
        return None


def AFP_regfilter(df, shp_Reg, strlat="Lat", strlon="Lon"):
    """ filter fire pixels using a given shp_Reg

    Parameters
    ----------
    df : pandas DataFrame
        fire pixel data with lat and lon information
    shp_Reg : geometry
        the geometry of the region shape
    strlat : str
        the column name for the latitude
    strlon : str
        the column name for the longitude

    Returns
    -------
    df_filtered : pandas DataFrame
        the filtered fire pixels
    """
    from FireConsts import remove_static_sources_bool
    from shapely.geometry import Point
    import geopandas as gpd
    import pandas as pd

    # preliminary spatial filter and quality filter
    regext = shp_Reg.bounds
    newfirepixels = df.loc[
        (df[strlat] >= regext[1])
        & (df[strlat] <= regext[3])
        & (df[strlon] >= regext[0])
        & (df[strlon] <= regext[2])
    ]
    point_data = [Point(xy) for xy in zip(newfirepixels[strlon], newfirepixels[strlat])]
    gdf_filtered = gpd.GeoDataFrame(newfirepixels, geometry=point_data, crs=4326)

    # if shp_Reg is not a rectangle, do detailed filtering (within shp_Reg)
    import math

    if ((not math.isclose(shp_Reg.minimum_rotated_rectangle.area, shp_Reg.area)) | remove_static_sources_bool):
        gdf_filtered = gdf_filtered[gdf_filtered["geometry"].within(shp_Reg)]

    # drop geometry column
    from FireConsts import epsg
    
    #print('current proj:',gdf_filtered.crs)
    
    df_filtered = AFP_toprj(gdf_filtered, epsg=epsg)
    
    #print('converting to proj:',epsg)
    # df_filtered = pd.DataFrame(gdf_filtered.drop(columns='geometry'))

    return df_filtered


def AFP_nonstatfilter(df):
    """ Function used to filter non-static VIIRS active fire data

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
    """ Function to set ampm column (using local hour) and update the df

    Parameters
    ----------
    df : pandas DataFrame
        fire pixel data, with 'Lon' and 'YYYYMMDD_HHMM' column

    Returns
    -------
    df_withampm : pandas DataFrame
        the DataFrame with 'ampm' column
    """
    import pandas as pd
    import numpy as np

    # calculate local hour using the longitude and YYYYMMDD_HHMM column
    localhour = (
        pd.to_timedelta(df.Lon / 15, unit="hours") + df["YYYYMMDD_HHMM"]
    ).dt.hour

    # set am/pm flag based on local hour
    df_withampm = df.assign(
        ampm=np.where(((localhour > 6) & (localhour < 18)), "PM", "AM")
    )

    return df_withampm


def AFP_toprj(gdf, epsg=32610):
    """Transforms lat/lon coordinates to projected coordinate system for computations in m
    for global studies probably proj:cea would be a good choice
    for the boreals we use North Pole LAEA (epsg:3571)
    for CA may use WGS 84 / UTM zone 10N (epsg: 32610)
    for US may use US National Atlas Equal Area (epsg: 9311)"""

    import pandas as pd

    gdf = gdf.to_crs(epsg=epsg)
    gdf["x"] = gdf.geometry.x
    gdf["y"] = gdf.geometry.y
    df = pd.DataFrame(gdf.drop(columns="geometry"))

    return df


def save_AFPtmp(df, d, head="VNP14IMGML.", tail=".pkl"):
    """ Function used to save regional filtered data to a temporary directory
    """
    from FireConsts import dirtmpdata
    import os
    import pickle
    import fsspec

    # temporary file name
    fnm_tmp = os.path.join(dirtmpdata, head + d.strftime("%Y%m") + tail)
    check_filefolder(fnm_tmp)

    # save
    with fsspec.open(fnm_tmp, "wb") as f:
        pickle.dump(df, f)


def load_AFPtmp(d, head="VNP14IMGML.", tail=".pkl"):
    """ Function used to load regional filtered data to a temporary directory
    """
    from FireConsts import dirtmpdata
    import os
    import fsspec
    import pickle

    # temporary file name
    fnm_tmp = os.path.join(dirtmpdata, head + d.strftime("%Y%m") + tail)

    # load and return
    if os_path_exists(fnm_tmp):
        with fsspec.open(fnm_tmp, "rb") as f:
            df = pickle.load(f)
            return df
    else:
        print(f"Temporary file not found: {fnm_tmp}.")
        return None


def read_AFPVIIRS(
    t,
    region,
    sat="SNPP",
    cols=["Lat", "Lon", "FRP", "Sat", "DT", "DS", "YYYYMMDD_HHMM", "ampm", "x", "y"],
):
    """ Read half-daily fire location for a region from monthly NOAA20 VIIRS fire data.
        For the monthly VIIRS data, in order to reduce repeat calculation, we
        do the spatial and quality filtering once and save the filtered data to
        files in a temporary directory. For most time steps, we read the pre-saved
        VIIRS data.

    Parameters
    ----------
    t : tuple, (int,int,int,str)
        the year, month, day and 'AM'|'PM' during the intialization
    region : 2-element tuple [regnm, regshp].
        Regnm is the name of the region; Regshp can be
         - a geometry
         - a four-element list showing the extent of the region [lonmin,latmin,lonmax,latmax]
         - a country name
    cols : list
        the column names of the output dataframe
    Returns
    -------
    df_AFP : pandas DataFrame
        list of active fire pixels (filtered) from NOAA20
    """

    from datetime import date

    # read from pre-saved file
    d = date(*t[:-1])
    if sat == "SNPP":
        sathead = "VNP14IMGML"
    elif sat == "NOAA20":
        sathead = "VJ114IMGML"
    else:
        print("please set SNPP or NOAA20 for sat")

    # load from pre-saved file
    df = load_AFPtmp(d, head=region[0] + "_" + sathead + ".")

    # if no pre-saved file, read from original VNP14IMGML file and do/save spatial filtering
    if df is None:
        # read or form shape used for filtering active fires
        shp_Reg = get_reg_shp(region[1])
        
        # read VJ114IMGML monthly file
        if sat == "SNPP":
            df = read_VNP14IMGML(t[0], t[1])
            df = df.loc[df["Type"] == 0]  # type filtering
        elif sat == "NOAA20":
            df = read_VJ114IMGML(t[0], t[1])
            df = df.loc[df["mask"] >= 7]
        else:
            print("please set SNPP or NOAA20 for sat")

        # do regional filtering
        df = AFP_regfilter(df, shp_Reg)
        # set ampm
        df = AFP_setampm(df)

        # save to temporary file
        save_AFPtmp(df, d, head=region[0] + "_" + sathead + ".")

    # extract active pixels at current time step  (day and ampm filter)
    df = df.loc[(df["YYYYMMDD_HHMM"].dt.day == t[2]) & (df["ampm"] == t[-1])]

    # add the satellite information
    df["Sat"] = sat

    # return selected columns; need to change column names first
    df_AFP = df[cols]

    return df_AFP


def read_AFPVIIRSNRT(
    t,
    region,
    sat="SNPP",
    cols=["Lat", "Lon", "FRP", "Sat", "DT", "DS", "YYYYMMDD_HHMM", "ampm", "x", "y"],):
    
    """ Read half-daily active fire pixels in a region from daily NRT NOAA20 VIIRS
    fire location data

    Parameters
    ----------
    t : tuple, (int,int,int,str)
        the year, month, day and 'AM'|'PM' during the intialization
    region : 2-element tuple [regnm, regshp].
        Regnm is the name of the region; Regshp can be
         - a geometry
         - a four-element list showing the extent of the region [lonmin,latmin,lonmax,latmax]
         - a country name
    cols : list
        the column names of the output dataframe
    Returns
    -------
    df_AFP : pandas DataFrame
        list of active fire pixels (filtered) from SNPPNRT
    """
    from datetime import date

    # read from VNP14IMGTDL file
    if sat == "SNPP":
        df = read_VNP14IMGTDL(t)
    elif sat == "NOAA20":
        df = read_VJ114IMGTDL(t)
    else:
        print("please set SNPP or NOAA20 for sat")

    if df is None:
        return None

    # extract active pixes at current time step
    df = AFP_setampm(df)
    df = df.loc[(df["ampm"] == t[-1])]

    # do regional filtering
    shp_Reg = get_reg_shp(region[1])
    df = AFP_regfilter(df, shp_Reg)

    # do non-static filtering; now only available at CA
    # df = AFP_nonstatfilter(df)

    # add the satellite information
    df["Sat"] = sat

    # return selected columns; need to change column names first
    df_AFP = df[cols]

    return df_AFP


def read_AFP(t, src="SNPP", nrt=False, region=None):
    """ The wrapper function used to read and extract half-daily active fire
    pixels in a region.

    Parameters
    ----------
    t : tuple, (int,int,int,str)
        the year, month, day and 'AM'|'PM' during the intialization
    src : str
        'SNPP' | 'NOAA20' | 'VIIRS' | 'BAMOD'
    nrt : bool
        if set to True, read NRT data instead of monthly data
    region : 2-element tuple [regnm, regshp].
        Regnm is the name of the region; Regshp can be
         - a geometry
         - a four-element list showing the extent of the region [lonmin,latmin,lonmax,latmax]
         - a country name
    Returns
    -------
    vlist : pandas DataFrame
        extracted active fire pixel data used for fire tracking
    """
    import pandas as pd

    if src == "VIIRS":
        if nrt:
            vlist_SNPP = read_AFPVIIRSNRT(t, region, sat="SNPP")
            vlist_NOAA20 = read_AFPVIIRSNRT(t, region, sat="NOAA20")
            vlist = pd.concat([vlist_NOAA20, vlist_SNPP], ignore_index=True)
        else:
            vlist_SNPP = read_AFPVIIRS(t, region, sat="SNPP")
            vlist_NOAA20 = read_AFPVIIRS(t, region, sat="NOAA20")
            vlist = pd.concat([vlist_NOAA20, vlist_SNPP], ignore_index=True)
    elif src == "SNPP":
        if nrt:
            vlist = read_AFPVIIRSNRT(t, region, sat="SNPP")
        else:
            vlist = read_AFPVIIRS(t, region, sat="SNPP")
    elif src == "NOAA20":
        if nrt:
            vlist = read_AFPVIIRSNRT(t, region, sat="NOAA20")
        else:
            vlist = read_AFPVIIRS(t, region, sat="NOAA20")
    elif src == "BAMOD":
        vlist = read_BAMOD(t, region)
    else:
        print("Please set src to SNPP, NOAA20, or BAMOD")
        return None
    
    if vlist is None: 
        print("No data available for this source for time",t) 
        return
    
    return vlist.reset_index(drop=True)


# ------------------------------------------------------------------------------
#%% Read other datasets
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
        from modis burned area """

    from FireConsts import dirextdata

    # import glob
    import numpy as np
    import pandas as pd
    import os

    mcd64dir = os.path.join(dirextdata, "MCD64A1")
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

    from FireConsts import dirextdata
    import gdal
    import os

    mcd64dir = os.path.join(dirextdata, "MCD64A1")
    fnm = os.path.join(mcd64dir, "mcd64_" + str(year) + ".tif")
    ds = gdal.Open(fnm)
    # arr = ds.ReadAsArray(xsize=xsize, ysize=ysize)
    arr = ds.ReadAsArray(xoff=xoff, yoff=yoff, xsize=xsize, ysize=ysize)

    return arr


def get_any_shp(filename):
    """ get shapefile of any region given the input file name

    Parameters
    ----------
    filename : str
        the shapefile names saved in the directory dirextdata/shapefiles/
    """
    from FireConsts import dirextdata
    import geopandas as gpd
    import os

    # find the california shapefile
    dirshape = os.path.join(dirextdata, "shapefiles")
    statefnm = os.path.join(dirshape, filename)

    # read the geometry
    shp = gpd_read_file(statefnm).iloc[0].geometry

    return shp


def get_Cal_shp():
    """ get shapefile of California
    !!! this function can be realized using get_reg_shp(); still keep here for convenience...;
        will be deleted or modified later !!!
    """
    from FireConsts import dirextdata
    import geopandas as gpd
    import os

    # find the california shapefile
    statefnm = os.path.join(dirextdata, "CA", "Calshape", "California.shp")

    # read the geometry
    shp_Cal = gpd_read_file(statefnm).iloc[0].geometry

    return shp_Cal


def get_Cty_shp(ctr):
    """ get shapefile of a country

    Parameters
    ----------
    ctr : str
        country name
    """
    import geopandas as gpd
    import os
    from FireConsts import dirextdata

    ctyfnm = os.path.join(dirextdata, "World", "country.shp")

    gdf_cty = gpd_read_file(ctyfnm)

    if ctr in gdf_cty["CNTRY_NAME"].values:
        g = gdf_cty[gdf_cty.CNTRY_NAME == ctr].iloc[0].geometry
        return g
    else:
        return None


def get_reg_shp(reg):
    """ return the shape of a region, given an optional reg input

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
    from shapely.geometry import Point, Polygon
    import shapely

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
    """ Get land cover type for active fires

    Parameters
    ----------
    locs : list of lists (nx2)
        lat and lon values for each active fire detection

    Returns
    -------
    vLCT : list of ints
        land cover types for all input active fires
    """
    from FireConsts import dirextdata

    import rasterio
    import pyproj
    import os

    # read NLCD 500m data
    fnmLCT = os.path.join(dirextdata, "CA", "nlcd_510m.tif")
    dataset = rasterio.open(fnmLCT)
    transformer = pyproj.Transformer.from_crs("epsg:4326", dataset.crs)
    locs_crs_x, locs_crs_y = transformer.transform(
        # NOTE: EPSG 4326 expected coordinate order latitude, longitude, but
        # `locs` is x (longitude), y (latitude). That's why `l[1]`, then `l[0]`
        # here.
        [l[1] for l in locs],
        [l[0] for l in locs]
    )
    locs_crs = list(zip(locs_crs_x, locs_crs_y))
    samps = list(dataset.sample(locs_crs))
    vLCT = [int(s) for s in samps]
    return vLCT

def get_LCT_CONUS(locs):
    """ Get land cover type for active fires - CONUS scale.
        This is the same function as get_LCT but with a CONUS wide file.
    
    Parameters
    ----------
    locs : list of lists (nx2)
        lat and lon values for each active fire detection

    Returns
    -------
    vLCT : list of ints
        land cover types for all input active fires
    """
    from FireConsts import dirextdata

    import rasterio
    import pyproj
    import os

    # read NLCD 500m data
    fnmLCT = os.path.join(dirextdata, "NLCD", "nlcd_export_510m_simplified.tif")
    dataset = rasterio.open(fnmLCT)
    transformer = pyproj.Transformer.from_crs("epsg:4326", dataset.crs)
    locs_crs_x, locs_crs_y = transformer.transform(
        # NOTE: EPSG 4326 expected coordinate order latitude, longitude, but
        # `locs` is x (longitude), y (latitude). That's why `l[1]`, then `l[0]`
        # here.
        [l[1] for l in locs],
        [l[0] for l in locs]
    )
    locs_crs = list(zip(locs_crs_x, locs_crs_y))
    samps = list(dataset.sample(locs_crs))
    vLCT = [int(s) for s in samps]
    return vLCT

def get_LCT_NLCD(locs):
    """ Get land cover type from NCLD for multiple locations

    Parameters
    ----------
    locs : list of lists (nx2)
        lon and lat values for each active fire detection

    Returns
    -------
    vLCT : list of ints
        land cover types for all input active fires
    """
    from FireConsts import dirextdata
    import rasterio

    # from osgeo import gdal
    # import pyproj
    import os

    # read NLCD 500m data
    fnmLCT = os.path.join(dirextdata, "CA", "nlcd_510m_latlon.tif")
    dataset = rasterio.open(fnmLCT)
    # locs = [(x,y) for (y,x) in locs]
    vLCT = dataset.sample(locs, indexes=1)
    vLCT = [lc[0] for lc in vLCT]  # list values

    return vLCT


def get_FM1000(t, loc):
    """ Get fm1000 for a point at t

    Parameters
    ----------
    t : datetime date
        date
    loc : tuple
        (lat,lon) value

    Returns
    -------
    FM1000_loc : list of floats
        fm1000 value for all input active fires
    """
    from FireConsts import dirextdata

    import xarray as xr
    import os

    import warnings

    warnings.simplefilter("ignore")

    # read annual fm1000 data
    dirGridMET = os.path.join(dirextdata, "GridMET") + "/"
    fnm = dirGridMET + "fm1000_" + t.strftime("%Y") + ".zarr"
    ds = xr.open_zarr(fnm)
    FM1000_all = ds["dead_fuel_moisture_1000hr"]

    # extract daily data at t
    try:
        FM1000_day = FM1000_all.sel(day=t.strftime("%Y-%m-%d"))
    except:  # if data are not available, use the last available date
        FM1000_day = FM1000_all.isel(day=-1)

    # extract data near the given location
    FM1000_loc = FM1000_day.sel(lon=loc[0], lat=loc[1], method="nearest").item()

    return FM1000_loc


# def get_stFM1000(fhull,locs,t):
#     ''' get FM1000 for a fire at a time
#
#     Parameters
#     ----------
#     fhull : geometry
#         the hull of the fire
#     locs : list of lists (nx2)
#         lat and lon values for each active fire detection
#     t : tuple, (int,int,int,str)
#         the year, month, day and 'AM'|'PM'
#
#     Returns
#     -------
#     FM1000 : list of floats
#         fm1000 value for all input active fires
#     '''
#     import FireClustering
#
#     from datetime import date
#
#     # centroid (lat,lon) is defined as that of hull (faster) or locs
#     if fhull is not None:
#         cent = (fhull.centroid.y,fhull.centroid.x)
#     else:
#         cent = FireClustering.cal_centroid(locs)
#
#     # datetime date of the time tuple
#     d_st = date(*t[:-1])
#
#     # call get_FM1000 function to extract fm1000 for the centroid at a time
#     FM1000 = get_FM1000(d_st, cent)
#
#     return FM1000

# ------------------------------------------------------------------------------
#%% read and load object, gdf and summary related files
# ------------------------------------------------------------------------------
def check_filefolder(fnm):
    """ if the folder containing a file does not exist, make it

    Parameters
    ----------
    fnm : str
        file name
    """
    import os

    # folder name
    dirnm = os.path.dirname(fnm)
    if dirnm.startswith("s3://"):
        # No concept of "folders" in S3, so no need to create them.
        return None

    # create the folder if needed
    if not os.path.exists(dirnm):
        os.makedirs(dirnm)


def check_folder(dirnm):
    """ if the folder does not exist, make it

    Parameters
    ----------
    dirfnm : str
        folder name
    """
    import os

    if dirnm.startswith("s3://"):
        # No concept of "folders" in S3, so no need to create them.
        return None

    # create the folder if needed
    if not os.path.exists(dirnm):
        os.makedirs(dirnm)


def correct_dtype(gdf, op=""):
    """ correct the datatype for gdfs loaded from geojson files

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
    from FireConsts import dd

    # explicitly set the attributes data types
    if op == "":
        for v, tp in dd.items():
            gdf[v] = gdf[v].astype(tp)
    else:
        gdf["fireID"] = gdf["fireID"].astype("int")

    return gdf


def get_fobj_fnm(t, regnm, activeonly=False):
    """ Return the fire object pickle file name at a time step
    Parameters
    ----------
    t : tuple, (year,month,day,str)
        the day and 'AM'|'PM' during the intialization

    Returns
    ----------
    fnm : str
        pickle file name
    """
    from FireConsts import diroutdata
    from datetime import date
    import os

    d = date(*t[:-1])
    # fnm = dirpjdata+regnm+'/'+d.strftime('%Y')+'/Serialization/'+d.strftime('%Y%m%d')+t[-1]+'.pkl'
    if activeonly:
        fnm = os.path.join(
            diroutdata,
            regnm,
            d.strftime("%Y"),
            "Serialization",
            d.strftime("%Y%m%d") + t[-1] + ".pkl",
        )
    else:
        fnm = os.path.join(
            diroutdata,
            regnm,
            d.strftime("%Y"),
            "Serialization",
            d.strftime("%Y%m%d") + t[-1] + "_full.pkl",
        )
    return fnm


def check_fobj(t, regnm, activeonly=False):
    """ Check if the pickle file storing a daily allfires object exists

    Parameters
    ----------
    t : tuple, (year,month,day,str)
        the day and 'AM'|'PM' during the intialization
    """
    import os

    fnm = get_fobj_fnm(t, regnm, activeonly=activeonly)
    return os_path_exists(fnm)


def save_fobj(allfires, t, regnm, activeonly=False):
    """ Save a daily allfires object to a pickle file

    Parameters
    ----------
    allfires : obj of Allfires class
        daily allfires
    t : tuple, (year,month,day,str)
        the day and 'AM'|'PM' during the intialization
    activeonly : bool
        the flag to save activeonly or full pickle file
    """

    import pickle
    import os
    import fsspec
    from FireObj import Allfires

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
    """ Load a daily allfires object from a pickle file

    Parameters
    ----------
    t : tuple, (year,month,day,str)
        the day and 'AM'|'PM' during the intialization

    Returns
    ----------
    data : obj of Allfires class
        daily allfires
    """
    import pickle
    import fsspec

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
    """ Return the fire object gpkg file name at a time step
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
    from FireConsts import diroutdata
    from datetime import date
    import os

    d = date(*t[:-1])
    if op == "":
        fnm = os.path.join(
            diroutdata,
            regnm,
            d.strftime("%Y"),
            "Snapshot",
            d.strftime("%Y%m%d") + t[-1] + ".gpkg",
        )
    else:
        fnm = os.path.join(
            diroutdata,
            regnm,
            d.strftime("%Y"),
            "Snapshot",
            d.strftime("%Y%m%d") + t[-1] + "_" + op + ".gpkg",
        )

    return fnm


def get_gpkgobj_fnm(t, regnm):
    """ Return the gpkg snapshot file name at a time step
    Parameters
    ----------
    t : tuple, (year,month,day,str)
        the day and 'AM'|'PM' during the intialization
    Returns
    ----------
    fnm : str
        gpkg file name
    """
    from FireConsts import diroutdata
    from datetime import date
    import FireIO
    import os

    # determine output dir
    d = date(*t[:-1])
    strdir = os.path.join(diroutdata, regnm, d.strftime("%Y"), "Snapshot")

    # get the output file name
    fnm = os.path.join(strdir, d.strftime("%Y%m%d") + t[-1])

    return fnm


def get_gpkgsfs_fnm(t, fid, regnm):
    """ Return the gpkg snapshot file name at a time step
    Parameters
    ----------
    t : tuple, (year,month,day,str)
        the day and 'AM'|'PM' during the intialization
    Returns
    ----------
    fnm : str
        gpkg file name
    """
    from FireConsts import diroutdata
    from datetime import date
    import FireIO
    import os

    # determine output dir
    d = date(*t[:-1])
    strdir = os.path.join(diroutdata, regnm, d.strftime("%Y"), "Largefire")

    # get the output file name
    fnm = os.path.join(
        strdir, "F" + str(int(fid)) + "_" + d.strftime("%Y%m%d") + t[-1]
    )

    return fnm

def get_gpkgsfs_dir(yr, regnm):
    """ Return the gpkg snapshot directory name for a year
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
    from FireConsts import diroutdata
    import os

    # determine output dir
    strdir = os.path.join(diroutdata, regnm, str(yr), "Largefire")

    return strdir

def get_NFPlistsfs_fnm(t, fid, regnm):
    """ Return the gpkg snapshot file name at a time step
    Parameters
    ----------
    t : tuple, (year,month,day,str)
        the day and 'AM'|'PM' during the intialization
    Returns
    ----------
    fnm : str
        gpkg file name
    """
    from FireConsts import diroutdata
    from datetime import date
    import FireIO
    import os

    # determine output dir
    d = date(*t[:-1])
    strdir = os.path.join(diroutdata, regnm, d.strftime("%Y"), "Largefire")

    # get the output file name
    fnm = os.path.join(
        strdir, "F" + str(int(fid)) + "_" + d.strftime("%Y%m%d") + t[-1] + "_NFP.txt"
    )

    return fnm


def check_gpkgobj(t, regnm):
    """ Check if the gpkg file storing a daily allfires attributes exists

    Parameters
    ----------
    t : tuple, (year,month,day,str)
        the day and 'AM'|'PM' during the intialization
    regnm : str
        the name of the region
    """
    from datetime import date
    import os

    # d = date(*t[:-1])
    fnm = get_gpkgobj_fnm(t, regnm)

    return os_path_exists(fnm)


def check_gdfobj(t, regnm, op=""):
    """ Check if the gpkg file storing a daily allfires attributes exists

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
    from datetime import date
    import os

    d = date(*t[:-1])
    fnm = get_gdfobj_fnm(t, regnm)

    return os_path_exists(fnm)


def save_gdfobj(gdf, t, regnm, param="", fid="", op=""):
    """ Save geopandas to a gpgk file

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
    import os

    # get file name
    if param == "":
        fnm = get_gdfobj_fnm(t, regnm, op=op)
    elif param == "large":
        fnm = get_gdfobj_sf_fnm(t, fid, regnm, op=op)
    else:
        from datetime import date
        from FireConsts import diroutdata

        d = date(*t[:-1])

        # get file name
        fnm = os.path.join(
            diroutdata,
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
    """ Save geopandas to a gpkg fire object file

    Parameters
    ----------
    gdf : geopandas DataFrame
        daily allfires diagnostic parameters
    t : tuple, (int,int,int,str)
        the year, month, day and 'AM'|'PM'
    """
    from datetime import date

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
        copy_from_maap_to_veda_s3(f"{fnm}/perimeter.fgb")

    if gdf_fline is not None:
        gdf_fline.to_file(f"{fnm}/fireline.fgb", driver="FlatGeobuf")
        copy_from_maap_to_veda_s3(f"{fnm}/fireline.fgb")

    if gdf_nfp is not None:
        gdf_nfp.to_file(f"{fnm}/newfirepix.fgb", driver="FlatGeobuf")
        copy_from_maap_to_veda_s3(f"{fnm}/newfirepix.fgb")

    if gdf_uptonow is not None:
        gdf_uptonow.to_file(f"{fnm}/uptonow.fgb", driver="FlatGeobuf")


def save_gpkgsfs(
    t, fid, regnm, gdf_fperim=None, gdf_fline=None, gdf_nfp=None, gdf_nfplist=None
):
    """ Save geopandas to a gpkg fire object file

    Parameters
    ----------
    gdf : geopandas DataFrame
        daily allfires diagnostic parameters
    t : tuple, (int,int,int,str)
        the year, month, day and 'AM'|'PM'
    """
    from datetime import date

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
        copy_from_maap_to_veda_s3(f"{fnm}/perimeter.fgb")

    # if len(gdf_fline) > 0:
    if gdf_fline is not None:
        gdf_fline.to_file(f"{fnm}/fireline.fgb", driver="FlatGeobuf")
        copy_from_maap_to_veda_s3(f"{fnm}/fireline.fgb")

    # if len(gdf_NFP) > 0:
    if gdf_nfp is not None:
        gdf_nfp.to_file(f"{fnm}/newfirepix.fgb", driver="FlatGeobuf")
        copy_from_maap_to_veda_s3(f"{fnm}/newfirepix.fgb")

    if gdf_nfplist is not None:
        gdf_nfplist.to_file(f"{fnm}/nfplist.fgb", driver="FlatGeobuf")
        copy_from_maap_to_veda_s3(f"{fnm}/nfplist.fgb")

def load_gpkgobj(t, regnm, layer="perimeter"):
    """ Load geopandas from a gpkg fire object file
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
    import os
    import pandas as pd
    from time import sleep
    
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
        print('Encountered the following error:',e)
        gdf = None
        return gdf
    
def load_gpkgsfs(t, fid, regnm, layer="perimeter"):
    """ Load geopandas from a gpkg fire object file

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
    import os

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
    """ Load daily allfires diagnostic dataframe as geopandas gdf

    Parameters
    ----------
    t : tuple, (int,int,int,str)
        the year, month, day and 'AM'|'PM'
    Returns
    ----------
    gdf : geopandas DataFrame
        daily allfires diagnostic parameters
    """
    import geopandas as gpd

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
    """ Save the cross year fid mapping tuples
    """
    from FireConsts import diroutdata
    import os
    import pandas as pd

    # convert list to dataframe
    df = pd.DataFrame(fidmapping, columns=["oldfid", "newfid"])

    # determine output file name
    strdir = os.path.join(diroutdata, regnm, str(year), "Summary")
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
    import os
    import pandas as pd

    # get filename of new fire pixels product
    fnm = get_gdfobj_fnm(t, regnm, op="NFP")
    fnm = fnm[:-4] + "txt"  # change ending to txt

    if os_path_exists(fnm):
        df = pd.read_csv(fnm, parse_dates=["datetime"], index_col=0)
        return df


def load_lake_geoms(t, fid, regnm):
    """ Load final perimeters as geopandas gdf

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
    import geopandas as gpd
    from datetime import date
    from FireConsts import diroutdata
    import os

    d = date(*t[:-1])

    # get file name
    fnm_lakes = os.path.join(
        diroutdata,
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
    """ Return the single fire fire object pickle file name at a time step
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
    from FireConsts import diroutdata
    from datetime import date
    import os

    d = date(*t[:-1])

    if op == "":
        fnm = os.path.join(
            diroutdata,
            regnm,
            d.strftime("%Y"),
            "Largefire",
            "F" + str(int(fid)) + "_" + d.strftime("%Y%m%d") + t[-1] + ".gpkg",
        )
    else:
        fnm = os.path.join(
            diroutdata,
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
    """ Return the single fire fire object pickle file name at a time step
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
    from FireConsts import diroutdata
    from datetime import date
    from glob import glob
    import os

    if op == "":
        fnms = glob(
            os.path.join(
                diroutdata,
                regnm,
                str(year),
                "Largefire",
                "F" + str(int(fid)) + "_??????????.gpkg",
            )
        )
    else:
        fnms = glob(
            os.path.join(
                diroutdata,
                regnm,
                str(year),
                "Largefire",
                "F" + str(int(fid)) + "_??????????_" + op + ".gpkg",
            )
        )

    return fnms


def save_gdfobj_sf(gdf, t, fid, regnm, op=""):
    """ Save daily single fire allfires diagnostic dataframe to a geojson file

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
    """ Load single fire daily allfires diagnostic dataframe to a geojson file

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
    import geopandas as gpd
    import pandas as pd

    # get file name
    fnm = get_gdfobj_sf_fnm(t, fid, regnm, op=op)

    # read data as gpd DataFrame
    gdf = gpd_read_file(fnm)

    # parse time step to Datetime and set as index
    gdf["index"] = pd.to_datetime(gdf["index"])
    gdf = gdf.set_index("index")

    return gdf


def get_summary_fnm(t, regnm):
    """ Return the fire summary file name at year end
    Parameters
    ----------
    t : tuple, (year,month,day,str)
        the day and 'AM'|'PM' during the intialization
    Returns
    ----------
    fnm : str
        summary netcdf file name
    """
    from FireConsts import diroutdata
    from datetime import date
    import os

    d = date(*t[:-1])
    fnm = os.path.join(
        diroutdata,
        regnm,
        d.strftime("%Y"),
        "Summary",
        "fsummary_" + d.strftime("%Y%m%d") + t[-1] + ".nc",
    )

    # check folder
    check_filefolder(fnm)

    return fnm


def get_summary_fnm_lt(t, regnm):
    """ Return the latest time step before current time when summary file exists
    Parameters
    ----------
    t : tuple, (year,month,day,str)
        the day and 'AM'|'PM' during the intialization

    Returns
    ----------
    pt : tuple, (year, month, day, pmpm)
        the lastest time step with summary file
    """
    from FireConsts import diroutdata
    from datetime import date
    import os
    from glob import glob
    import FireObj

    # if there's no summary file for this year, return the first time step of the year
    fnms = glob(os.path.join(diroutdata, regnm, str(t[0]), "Summary", "fsummary_*.nc"))
    if len(fnms) == 0:
        return None

    # if there is, find the nearest one
    endloop = False
    pt = FireTime.t_nb(t, nb="previous")
    while endloop == False:
        if os_path_exists(get_summary_fnm(pt, regnm)):
            return pt
        else:
            pt = FireTime.t_nb(pt, nb="previous")
            if pt[0] != t[0]:
                return None


def check_summary(t, regnm):
    """ Check if the summary file storing a daily allfires object exists

    Parameters
    ----------
    t : tuple, (year,month,day,str)
        the day and 'AM'|'PM' during the intialization
    """
    import os

    # get file name
    fnm = get_summary_fnm(t, regnm)

    # return if it's present
    return os_path_exists(fnm)


def save_summary(ds, t, regnm):
    """ Save summary info as of t a netcdf file

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
    """ Load summary info from a netcdf file at t

    Parameters
    ----------
    t : tuple, (int,int,int,str)
        the year, month, day and 'AM'|'PM'

    Returns
    -------
    ds : xarray dataset
        summary dataset
    """
    import xarray as xr

    # get file name
    fnm = get_summary_fnm(t, regnm)

    # read data as xarray Dataset
    ds = xr.open_dataset(fnm)

    return ds


def save_summarycsv(df, year, regnm, op="heritage"):
    """ save summary csv files

    Parameters
    ----------
    df : pandas DataFrame
        the data
    year : int
        year
    op : str
        option, 'heritage'|'large'
    """
    from FireConsts import diroutdata
    import os

    fnm = os.path.join(
        diroutdata,
        regnm,
        str(year),
        "Summary",
        "Flist_" + op + "_" + str(year) + ".csv",
    )
    check_filefolder(fnm)

    df.to_csv(fnm)


def read_summarycsv(year, regnm, op="heritage"):
    """ read summary csv files

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
    from FireConsts import diroutdata
    import pandas as pd
    import os

    fnm = os.path.join(
        diroutdata,
        regnm,
        str(year),
        "Summary",
        "Flist_" + op + "_" + str(year) + ".csv",
    )
    df = pd.read_csv(fnm, index_col=0)
    return df


def get_lts_VNP14IMGTDL(year=None):
    from FireConsts import dirextdata
    from datetime import date
    from glob import glob
    import os

    if year == None:
        year = date.today().year

    dirFC = os.path.join(dirextdata, "VNP14IMGTDL") + "/"
    fnms = glob(
        os.path.join(
            dirFC, "SUOMI_VIIRS_C2_Global_VNP14IMGTDL_NRT_" + str(year) + "*.txt"
        )
    )
    fnms.sort()
    DOY_lts = int(os.path.splitext(os.path.basename(fnms[-1]))[0][-3:])

    return DOY_lts


def get_lts_serialization(regnm, year=None):
    """ get the time with lastest pkl data
    """
    from FireConsts import diroutdata
    from datetime import date
    from glob import glob
    import os

    if year == None:
        year = date.today().year

    if diroutdata.startswith("s3://"):
        # Can't use glob for S3. Use s3.ls instead.
        import s3fs
        s3 = s3fs.S3FileSystem(anon=False)
        s3path = os.path.join(diroutdata, regnm, str(year), "Serialization")
        fnms = [f for f in s3.ls(s3path) if f.endswith(".pkl")]
    else:
        fnms = glob(os.path.join(diroutdata, regnm, str(year), "Serialization", "*.pkl"))

    if len(fnms) > 0:
        fnms.sort()
        fnm_lts = os.path.basename(fnms[-1])

        lts = [int(fnm_lts[0:4]), int(fnm_lts[4:6]), int(fnm_lts[6:8]), fnm_lts[8:10]]
    else:
        lts = None

    return lts


def read_gpkg(fnm, layer="perimeter"):
    """ read gpkg data
    op: 'perimeter', 'fireline', 'newfirepix'
    """
    import geopandas as gpd
    import pandas as pd

    try:
        gdf = gpd_read_file(fnm, layer=layer)
        return gdf
    except Exception as e:
        print(f"Failed to read GPKG file {fnm} with error: {str(e)}")
        return None


#%% other functions related to read/write


def save2gtif(arr, outfile, cols, rows, geotrans, proj):
    """write out a geotiff"""

    import gdal

    driver = gdal.GetDriverByName("GTiff")
    ds = driver.Create(outfile, cols, rows, 1, gdal.GDT_Byte)
    ds.SetGeoTransform(geotrans)
    ds.SetProjection(proj)
    band = ds.GetRasterBand(1)
    band.WriteArray(arr)


def geo_to_polar(lon_arr, lat_arr):
    """transform lists of geographic lat lon coordinates to polar LAEA grid (projected)"""

    import numpy as np
    import pyproj

    proj4str = "epsg:3571"
    p_modis_grid = pyproj.Proj(proj4str)

    x_arr, y_arr = p_modis_grid(lon_arr, lat_arr)
    x_arr = np.array(x_arr)
    y_arr = np.array(y_arr)

    return x_arr, y_arr


def polar_to_geo(x_arr, y_arr):
    """transform lists of geographic lat lon coordinates to polar LAEA grid (projected)"""

    import numpy as np
    import pyproj

    proj4str = "epsg:3571"
    p_modis_grid = pyproj.Proj(proj4str)

    lon_arr, lat_arr = p_modis_grid(x_arr, y_arr, inverse=True)
    lon_arr = np.array(lon_arr)
    lat_arr = np.array(lat_arr)

    return lon_arr, lat_arr


def world2Pixel(gt, Xgeo, Ygeo):
    """ Uses a geomatrix (gdal.GetGeoTransform()) to calculate the pixel
    location of a geospatial coordinate"""

    import numpy as np

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


def copy_from_maap_to_veda_s3(from_maap_s3_path):
    s3_client = boto3.client('s3')

    # the geopandas writes should block until finished writing to remote s3
    # but who knows what AWS does after that point to "finalize" persisted state
    # and have the file be "ready" for a read, so add a little buffer
    # TODO: if all writes were local this wouldn't be a problem
    time.sleep(15)

    if "Largefire" in from_maap_s3_path:
        try:
            fname_regex = r"^s3://maap.*?(/Largefire/)(?P<fid>F[0-9_a-zA-Z]+)/(?P<fname>fireline.fgb|perimeter.fgb|newfirepix.fgb|nfplist.fgb)$"
            # note that `destination_dict` should resemble this output with a match if the URL was a perimeter file:
            # {'fid': 'F1013_20230104PM', 'fname': 'perimeter.fgb'}
            destination_dict = re.compile(fname_regex).match(from_maap_s3_path).groupdict()
        except AttributeError:
            logger.error(f"[ NO REGEX MATCH FOUND ]: for file f{from_maap_s3_path}")
            return

        s3_client.copy_object(
            CopySource=from_maap_s3_path,  # full bucket path
            Bucket='veda-data-store-staging',  # Destination bucket
            Key=f"EIS/FEDSoutput/Largefire/{destination_dict['fid']}/{destination_dict['fname']}"  # Destination path/filename
        )

    elif "Snapshot" in from_maap_s3_path:
        try:
            fname_regex = r"^s3://maap.*(?P<fname>fireline.fgb|perimeter.fgb|newfirepix.fgb)$"
            # note that `destination_fname` should resemble this output with a match if the URL was a perimeter file:
            # {'fname': 'perimeter.fgb'}
            destination_fname = re.compile(fname_regex).match(from_maap_s3_path).groupdict()['fname']
        except AttributeError:
            logger.error(f"[ NO REGEX MATCH FOUND ]: for file f{from_maap_s3_path}")
            return

        s3_client.copy_object(
            CopySource=from_maap_s3_path,  # full bucket path
            Bucket='veda-data-store-staging',  # destination bucket
            Key=f'EIS/FEDSoutput/Snapshot/{destination_fname}'  # destination path/filename
        )
    else:
        logger.error(f"[ NO S3 COPY EXPORTED ]: for file f{from_maap_s3_path}")


import os
import uuid
import pandas as pd
from shapely import to_geojson, from_geojson

from FireConsts import CONNECTIVITY_CLUSTER_KM, dirextdata, dirprpdata
from FireClustering import do_clustering
from FireLog import logger
from FireIO import read_VNP14IMGTDL, read_VJ114IMGTDL, read_VNP14IMGML, read_VJ114IMGML, get_reg_shp, AFP_regfilter, AFP_setampm
from FireMain import maybe_remove_static_sources
from utils import timed

# INPUT_DIR = "/projects/shared-buckets/gsfc_landslides/FEDSinput/VIIRS/"
# OUTPUT_DIR = "/projects/shared-buckets/jsignell/processed/"

INPUT_DIR = dirextdata
OUTPUT_DIR = dirprpdata


@timed
def preprocess_region(region):
    region = maybe_remove_static_sources(region, INPUT_DIR)
    with open(f"{OUTPUT_DIR}{region[0]}.json", "w") as f:
        f.write(to_geojson(region[1], indent=2))


@timed
def read_region(region):
    with open(f"{OUTPUT_DIR}{region[0]}.json", "r") as f:
        shape = from_geojson(f.read())
    return (region[0], shape)


@timed
def preprocess_landcover(filename="nlcd_export_510m_simplified.tif"):
    import rasterio
    from rasterio.warp import calculate_default_transform, reproject, Resampling

    fnmLCT = os.path.join(dirextdata, "NLCD", filename)
    dst_crs = f'EPSG:4326'

    with rasterio.open(fnmLCT) as src:
        transform, width, height = calculate_default_transform(
            src.crs, dst_crs, src.width, src.height, *src.bounds)
        kwargs = src.meta.copy()
        kwargs.update({
            'crs': dst_crs,
            'transform': transform,
            'width': width,
            'height': height
        })

        with rasterio.open(fnmLCT.replace(".tif", "_latlon.tif"), 'w', **kwargs) as dst:
            for i in range(1, src.count + 1):
                reproject(
                    source=rasterio.band(src, i),
                    destination=rasterio.band(dst, i),
                    src_transform=src.transform,
                    src_crs=src.crs,
                    dst_transform=transform,
                    dst_crs=dst_crs,
                    resampling=Resampling.nearest)


def preprocessed_filename(t, sat, *, region=None, suffix="", preprocessed_data_dir=OUTPUT_DIR):
    return os.path.join(
        preprocessed_data_dir,
        sat,
        *([] if region is None else [region[0]]),
        f"{t[0]}{t[1]:02}{t[2]:02}_{t[3]}{suffix}.txt"
    )


@timed
def preprocess_NRT_file(t, sat):
    logger.info(f"preprocessing NRT file for {t[0]}-{t[1]}-{t[2]}, {sat}")
    
    if sat == "SNPP":
        df = read_VNP14IMGTDL(t, input_data_dir=INPUT_DIR)
    elif sat == "NOAA20":
        df = read_VJ114IMGTDL(t, input_data_dir=INPUT_DIR)
    else:
        raise ValueError("please set SNPP or NOAA20 for sat")

    # set ampm
    df = AFP_setampm(df)
    
    # add the satellite information
    df["Sat"] = sat

    # return selected columns
    df = df[["Lat", "Lon", "FRP", "Sat", "DT", "DS", "YYYYMMDD_HHMM", "ampm"]]

    for ampm in ["AM", "PM"]:
        time_filtered_df = df.loc[df["ampm"] == ampm]

        output_filepath = preprocessed_filename((*t[:3], ampm), sat)

        # make nested path if necessary
        os.makedirs(os.path.dirname(output_filepath), exist_ok=True)
        
        # save active pixels at this time step (day and ampm filter)    
        time_filtered_df.to_csv(output_filepath, index=False)


@timed
def preprocess_monthly_file(t, sat):
    """Preprocess a monthly file for a given satellite"""
    logger.info(f"preprocessing monthly file for {t[0]}-{t[1]}, {sat}")

    if sat == "SNPP":
        df = read_VNP14IMGML(t, input_data_dir=INPUT_DIR)
        df = df.loc[df["Type"] == 0]  # type filtering
    elif sat == "NOAA20":
        df = read_VJ114IMGML(t, input_data_dir=INPUT_DIR)
        df = df.loc[df["mask"] >= 7]
    else:
        raise ValueError("please set SNPP or NOAA20 for sat")

    # set ampm
    df = AFP_setampm(df)

    # add the satellite information
    df["Sat"] = sat

    # return selected columns
    df = df[["Lat", "Lon", "FRP", "Sat", "DT", "DS", "YYYYMMDD_HHMM", "ampm"]]

    days = df["YYYYMMDD_HHMM"].dt.date.unique()
    for day in days:
        for ampm in ["AM", "PM"]:
            time_filtered_df = df.loc[(df["YYYYMMDD_HHMM"].dt.date == day) & (df["ampm"] == ampm)]

            output_filepath = preprocessed_filename((day.year, day.month, day.day, ampm), sat)
            
            # make nested path if necessary
            os.makedirs(os.path.dirname(output_filepath), exist_ok=True)

            # save active pixels at this time step (day and ampm filter)
            time_filtered_df.to_csv(output_filepath, index=False)


@timed
def read_preprocessed(t, sat, *, region=None):
    output_filepath = preprocessed_filename(t, sat ,region=region)
    return pd.read_csv(output_filepath)


@timed
def preprocess_region_t(t, sat, region):
    logger.info(f"filtering and clustering {t[0]}-{t[1]}-{t[2]} {t[3]}, {sat}, {region[0]}")
    if sat == "VIIRS":
        df = pd.concat(
            [
                read_preprocessed(t, sat="SNPP"),
                read_preprocessed(t, sat="NOAA20"),
            ], ignore_index=True
        )
    else:
        df = read_preprocessed(t, sat=sat)

    # do regional filtering
    shp_Reg = get_reg_shp(region[1])
    df = AFP_regfilter(df, shp_Reg)

    # return selected columns
    df = df[["Lat", "Lon", "FRP", "Sat", "DT", "DS", "YYYYMMDD_HHMM", "ampm", "x", "y"]]

    # do preliminary clustering using new active fire locations (assign cid to each pixel)
    df = do_clustering(df, CONNECTIVITY_CLUSTER_KM)

    # assign a uuid to each pixel
    df['uuid'] = [uuid.uuid4() for _ in range(len(df.index))]

    output_filepath = preprocessed_filename(t, sat, region=region)
    
    # make nested path if necessary
    os.makedirs(os.path.dirname(output_filepath), exist_ok=True)

    df.to_csv(output_filepath, index=False)

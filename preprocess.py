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


INPUT_DIR = FireConsts.dirextdata
OUTPUT_DIR = FireConsts.dirprpdata


def preprocessed_region_filename(region: Region):
    return os.path.join(OUTPUT_DIR, f"{region[0]}.json")


@timed
def preprocess_region(region: Region):
    from FireMain import maybe_remove_static_sources

    output_filepath = preprocessed_region_filename(region)
    
    # make path if necessary
    os.makedirs(os.path.dirname(output_filepath), exist_ok=True)

    region = maybe_remove_static_sources(region, INPUT_DIR)

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


@timed
def preprocess_NRT_file(t: TimeStep, sat: Literal["NOAA20", "SNPP"]):
    logger.info(f"preprocessing NRT file for {t[0]}-{t[1]}-{t[2]}, {sat}")

    if sat == "SNPP":
        df = FireIO.read_VNP14IMGTDL(t, INPUT_DIR)
    elif sat == "NOAA20":
        df = FireIO.read_VJ114IMGTDL(t, INPUT_DIR)
    else:
        raise ValueError("please set SNPP or NOAA20 for sat")

    # set ampm
    df = FireIO.AFP_setampm(df)

    # add the satellite information
    df["Sat"] = sat

    # return selected columns
    df = df[["Lat", "Lon", "FRP", "Sat", "DT", "DS", "YYYYMMDD_HHMM", "ampm"]]

    filtered_df_paths = []
    for ampm in ["AM", "PM"]:
        time_filtered_df = df.loc[df["ampm"] == ampm]

        output_filepath = preprocessed_filename((*t[:3], ampm), sat)

        # make nested path if necessary
        os.makedirs(os.path.dirname(output_filepath), exist_ok=True)

        # save active pixels at this time step (day and ampm filter)
        time_filtered_df.to_csv(output_filepath, index=False)

        filtered_df_paths.append(output_filepath)

    return filtered_df_paths


@timed
def preprocess_monthly_file(t: TimeStep, sat: Literal["NOAA20", "SNPP"]):
    """Preprocess a monthly file for a given satellite"""
    logger.info(f"preprocessing monthly file for {t[0]}-{t[1]}, {sat}")

    if sat == "SNPP":
        df = FireIO.read_VNP14IMGML(t, INPUT_DIR)
        df = df.loc[df["Type"] == 0]  # type filtering
    elif sat == "NOAA20":
        df = FireIO.read_VJ114IMGML(t, INPUT_DIR)
        df = df.loc[df["mask"] >= 7]
    else:
        raise ValueError("please set SNPP or NOAA20 for sat")

    # set ampm
    df = FireIO.AFP_setampm(df)

    # add the satellite information
    df["Sat"] = sat

    # return selected columns
    df = df[["Lat", "Lon", "FRP", "Sat", "DT", "DS", "YYYYMMDD_HHMM", "ampm"]]

    days = df["YYYYMMDD_HHMM"].dt.date.unique()
    for day in days:
        for ampm in ["AM", "PM"]:
            time_filtered_df = df.loc[
                (df["YYYYMMDD_HHMM"].dt.date == day) & (df["ampm"] == ampm)
            ]

            output_filepath = preprocessed_filename(
                (day.year, day.month, day.day, ampm), sat
            )

            # make nested path if necessary
            os.makedirs(os.path.dirname(output_filepath), exist_ok=True)

            # save active pixels at this time step (day and ampm filter)
            time_filtered_df.to_csv(output_filepath, index=False)


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

    # if regional output already exists, exit early so we don't reprocess
    output_filepath = preprocessed_filename(t, sensor, region=region)
    if os.path.exists(output_filepath):
        return output_filepath

    # do regional filtering
    shp_Reg = FireIO.get_reg_shp(region[1])
    df = FireIO.AFP_regfilter(df, shp_Reg)

    if df.empty:
        return

    # return selected columns
    df = df[["Lat", "Lon", "FRP", "Sat", "DT", "DS", "YYYYMMDD_HHMM", "ampm", "x", "y"]]

    # do preliminary clustering using new active fire locations (assign cid to each pixel)
    df = FireClustering.do_clustering(df, FireConsts.CONNECTIVITY_CLUSTER_KM)

    # assign a uuid to each pixel
    df["uuid"] = [uuid.uuid4() for _ in range(len(df.index))]

    # make nested path if necessary
    os.makedirs(os.path.dirname(output_filepath), exist_ok=True)

    df.to_csv(output_filepath, index=False)

    return output_filepath

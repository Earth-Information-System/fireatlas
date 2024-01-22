import os
import pandas as pd
from shapely import to_geojson, from_geojson

from FireConsts import CONNECTIVITY_CLUSTER_KM
from FireClustering import do_clustering
from FireIO import read_VNP14IMGTDL, read_VJ114IMGTDL, read_VNP14IMGML, read_VJ114IMGML, get_reg_shp, AFP_regfilter, AFP_setampm
from FireMain import maybe_remove_static_sources

# INPUT_DIR = "/projects/shared-buckets/gsfc_landslides/FEDSinput/VIIRS/"
# OUTPUT_DIR = "/projects/shared-buckets/jsignell/processed/"

INPUT_DIR = "data/"
OUTPUT_DIR = "processed/"

sat = "SNPP"
epsg = 9311 
region = ('CONUS', [-126.401171875,24.071240929282325,
                    -61.36210937500001,49.40003415463647])
t = (2022, 8, 19)


def preprocess_region(region):
    region = maybe_remove_static_sources(region, INPUT_DIR)
    with open(f"{OUTPUT_DIR}{region[0]}.json", "w") as f:
        f.write(to_geojson(region[1], indent=2))


def read_region(region):
    with open(f"{OUTPUT_DIR}{region[0]}.json", "r") as f:
        shape = from_geojson(f.read())
    return (region[0], shape)


def preprocessed_NRT_filename(t, region, sat, preprocessed_data_dir):
    return os.path.join(
        preprocessed_data_dir,
        sat,
        region[0],
        f"NRT_{t[0]}{t[1]:02}{t[2]:02}_{t[3]}.txt"
    )


def preprocessed_filename(t, region, sat, preprocessed_data_dir):
    return os.path.join(
        preprocessed_data_dir,
        sat,
        region[0],
        f"{t[0]}{t[1]:02}{t[2]:02}_{t[3]}.txt"
    )


def preprocess_NRT_file(t, region, sat):
    print(f"preprocessing NRT file for {t[0]}-{t[1]}-{t[2]}, {region[0]}, {sat}")
    
    if sat == "SNPP":
        df = read_VNP14IMGTDL(t, input_data_dir=INPUT_DIR)
    elif sat == "NOAA20":
        df = read_VJ114IMGTDL(t, input_data_dir=INPUT_DIR)
    else:
        print("please set SNPP or NOAA20 for sat")

    # do regional filtering
    shp_Reg = get_reg_shp(region[1])
    df = AFP_regfilter(df, shp_Reg)

    # set ampm
    df = AFP_setampm(df)
    
    # add the satellite information
    df["Sat"] = sat

    # return selected columns
    df = df[["Lat", "Lon", "FRP", "Sat", "DT", "DS", "YYYYMMDD_HHMM", "ampm", "x", "y"]]

    for ampm in ["AM", "PM"]:
        time_filtered_df = df.loc[df["ampm"] == ampm]

        # do preliminary clustering using new active fire locations (assign cid to each pixel)
        time_filtered_df = do_clustering(time_filtered_df, CONNECTIVITY_CLUSTER_KM)

        output_filepath = preprocessed_NRT_filename((*t[:3], ampm), region, sat, OUTPUT_DIR)

        # make nested path if necessary
        os.makedirs(os.path.dirname(output_filepath), exist_ok=True)
        
        # save active pixels at this time step (day and ampm filter)    
        time_filtered_df.to_csv(output_filepath, index=False)


def preprocess_monthly_file(t, region, sat):
    """Preprocess a monthly file for a given region and satellite"""
    print(f"preprocessing monthly file for {t[0]}-{t[1]}, {region[0]}, {sat}")

    df = None
    if sat in ["SNPP", "VIIRS"]:
        df = read_VNP14IMGML(t, input_data_dir=INPUT_DIR)
        df = df.loc[df["Type"] == 0]  # type filtering
        df_SNPP = df
    if sat in ["NOAA20", "VIIRS"]:
        df = read_VJ114IMGML(t, input_data_dir=INPUT_DIR)
        df = df.loc[df["mask"] >= 7]
        df_NOAA20 = df
    if sat == "VIIRS":
        df = pd.concat([df_SNPP, df_NOAA20], ignore_index=True)
    if df is None:
        print("please set SNPP or NOAA20 for sat")
    
    print("filtering...")
    # do regional filtering
    shp_Reg = get_reg_shp(region[1])
    df = AFP_regfilter(df, shp_Reg)

    print("setting ap and pm...")
    # set ampm
    df = AFP_setampm(df)

    # add the satellite information
    df["Sat"] = sat

    # return selected columns
    df = df[["Lat", "Lon", "FRP", "Sat", "DT", "DS", "YYYYMMDD_HHMM", "ampm", "x", "y"]]

    days = df["YYYYMMDD_HHMM"].dt.date.unique()
    for day in days:
        for ampm in ["AM", "PM"]:
            time_filtered_df = df.loc[(df["YYYYMMDD_HHMM"].dt.date == day) & (df["ampm"] == ampm)]

            print("clustering...")
            # do preliminary clustering using new active fire locations (assign cid to each pixel)
            time_filtered_df = do_clustering(time_filtered_df, CONNECTIVITY_CLUSTER_KM)

            output_filepath = preprocessed_filename((day.year, day.month, day.day, ampm), region, sat, OUTPUT_DIR)
            print(f"writing {output_filepath}")
            
            # make nested path if necessary
            os.makedirs(os.path.dirname(output_filepath), exist_ok=True)

            # save active pixels at this time step (day and ampm filter)
            time_filtered_df.to_csv(output_filepath, index=False)


def read_preprocessed(t, region, sat):
    output_filepath = preprocessed_filename(t, region, sat, OUTPUT_DIR)
    return pd.read_csv(output_filepath)
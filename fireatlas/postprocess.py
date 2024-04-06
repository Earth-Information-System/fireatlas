import os

import fsspec
import numpy as np
import pandas as pd
import geopandas as gpd

from shapely.ops import unary_union
import warnings

warnings.filterwarnings("ignore", "GeoSeries.notna", UserWarning)

from fireatlas.utils import timed
from fireatlas.FireTypes import Region, TimeStep, Location
from fireatlas.FireTime import t2dt, t_generator
from fireatlas.FireGpkg_sfs import getdd as singlefire_getdd
from fireatlas.FireGpkg import getdd as snapshot_getdd
from fireatlas import settings


def all_dir(tst: TimeStep, region: Region, location: Location = None):
    return os.path.join(settings.get_path(location), settings.OUTPUT_DIR, region[0], str(tst[0]))


def get_t_of_last_allfires_run(tst: TimeStep, ted: TimeStep, region: Region, location: Location = None):
    """Look at the files in a given location and figure out the t of the last
    allfires run. 
    
    Returns
    -------
    t: TimeStep
       latest t within range for which there are allfires and allpixels files
    """
    fs = fsspec.filesystem(location or settings.READ_LOCATION)
    all_filenames = {
        os.path.basename(f).split(".")[0]
        for f in [
            *fs.glob(os.path.join(all_dir(tst, region, location=location), "allpixels*")),
            *fs.glob(os.path.join(all_dir(tst, region, location=location), "allfires*")),
        ]
    }

    for t in list(t_generator(tst, ted))[::-1]:
        t_string = f"{t[0]}{t[1]:02}{t[2]:02}_{t[3]}"
        if all_filenames.issuperset({f"allfires_{t_string}", f"allpixels_{t_string}"}):
            return t


def allpixels_filepath(
    tst: TimeStep,
    ted: TimeStep,
    region: Region,
    location: Location = None,
):
    filename = f"allpixels_{ted[0]}{ted[1]:02}{ted[2]:02}_{ted[3]}.csv"
    return os.path.join(all_dir(tst, region, location), filename)


@timed
def save_allpixels(allpixels, tst: TimeStep, ted: TimeStep, region: Region):
    output_filepath = allpixels_filepath(tst, ted, region, location="local")

    # make path if necessary
    os.makedirs(os.path.dirname(output_filepath), exist_ok=True)

    allpixels[allpixels.t <= t2dt(ted)].to_csv(output_filepath)
    return output_filepath


@timed
def read_allpixels(
    tst: TimeStep,
    ted: TimeStep,
    region: Region,
    location: Location = None,
):
    filepath = allpixels_filepath(tst, ted, region, location=location)

    return pd.read_csv(filepath, index_col="uuid", parse_dates=["t"])


def allfires_filepath(
    tst: TimeStep,
    ted: TimeStep,
    region: Region,
    location: Location = None,
):
    filename = f"allfires_{ted[0]}{ted[1]:02}{ted[2]:02}_{ted[3]}.parq"
    return os.path.join(all_dir(tst, region, location), filename)


@timed
def save_allfires_gdf(allfires_gdf, tst: TimeStep, ted: TimeStep, region: Region):
    output_filepath = allfires_filepath(tst, ted, region, location="local")

    # make path if necessary
    os.makedirs(os.path.dirname(output_filepath), exist_ok=True)

    allfires_gdf.to_parquet(output_filepath)
    return output_filepath


@timed
def read_allfires_gdf(
    tst: TimeStep,
    ted: TimeStep,
    region: Region,
    location: Location = None,
):
    filepath = allfires_filepath(tst, ted, region, location=location)

    return gpd.read_parquet(filepath)


def snapshot_folder(
    region: Region,
    tst: TimeStep,
    ted: TimeStep,
    location: Location = None,
):
    return os.path.join(
        all_dir(tst, region, location=location),
        "Snapshot",
        f"{ted[0]}{ted[1]:02}{ted[2]:02}{ted[3]}",
    )


@timed
def save_snapshot_layers(allfires_gdf_t, region: Region, tst: TimeStep, ted: TimeStep):
    output_dir = snapshot_folder(region, tst, ted, location="local")
    os.makedirs(output_dir, exist_ok=True)

    dt = t2dt(ted)

    for layer in ["perimeter", "fireline", "newfirepix"]:
        columns = [col for col in snapshot_getdd(layer)]
        data = allfires_gdf_t[[*columns, "invalid", "fireID"]].copy()

        if layer == "perimeter":
            data["geometry"] = allfires_gdf_t["hull"]
        if layer == "newfirepix":
            data["geometry"] = allfires_gdf_t["nfp"]
        if layer == "fireline":
            data["geometry"] = allfires_gdf_t["fline"]

        data = data.set_geometry("geometry")
        data = data[data.geometry.notna() & ~data.geometry.is_empty]

        if layer == "perimeter":
            # figure out the fire state given current t
            data["isignition"] = dt == data["t_st"]
            data["t_inactive"] = (dt - data["t_ed"]).dt.days

            data["isactive"] = ~data["invalid"] & (data["t_inactive"] <= settings.maxoffdays)
            data["isdead"] = ~data["invalid"] & (data["t_inactive"] > settings.limoffdays)
            data["mayreactivate"] = (
                ~data["invalid"]
                & (settings.maxoffdays < data["t_inactive"])
                & (data["t_inactive"] <= settings.limoffdays)
            )

            # map booleans to integers
            for col in ["isignition", "isactive", "isdead", "mayreactivate"]:
                data[col] = data[col].astype(int)

            # apply filter flag
            data["geom_counts"] = (
                data[["fireID", "geometry"]]
                .explode(index_parts=True)
                .groupby(["fireID"])
                .nunique()["geometry"]
            )  # count number of polygons
            data["low_confidence_grouping"] = np.where(
                data["geom_counts"] > 5, 1, 0
            )  # if more than 5 geometries are present, flag it

        data["region"] = str(region[0])

        # primary key is: region + fireID + 12hr slice
        data["primarykey"] = (
            data["region"] + "|" + data["fireID"].astype(str) + "|" + dt.isoformat()
        )

        # drop the columns we don't actually need
        data = data.drop(columns=["invalid", "fireID"])

        data.to_file(os.path.join(output_dir, f"{layer}.fgb"), driver="FlatGeobuf")


@timed
def save_snapshots(allfires_gdf, region, tst, ted):
    gdf = allfires_gdf.reset_index()

    for t in t_generator(tst, ted):
        dt = t2dt(t)
        data = gdf[gdf.t <= dt].drop_duplicates("fireID", keep="last")
        save_snapshot_layers(data, region, tst, t)


@timed
def find_largefires(allfires_gdf):
    gdf = allfires_gdf.reset_index()

    last_seen = gdf.drop_duplicates("fireID", keep="last")
    last_large = last_seen[
        (last_seen.farea > settings.LARGEFIRE_FAREA) & (last_seen.invalid == False)
    ]
    return last_large.fireID.values


def largefire_folder(region: Region, fid, tst: TimeStep, location: Location = None):
    return os.path.join(
        all_dir(tst, region, location=location),
        "Largefire",
        str(fid)
    )


@timed
def save_fire_nplist(allpixels_fid, region, fid, tst):
    output_dir = largefire_folder(region, fid, tst, location="local")
    os.makedirs(output_dir, exist_ok=True)

    data = allpixels_fid[
        ["x", "y", "FRP", "DS", "DT", "ampm", "datetime", "Sat"]
    ].copy()
    data.columns = ["x", "y", "frp", "DS", "DT", "ampm", "datetime", "sat"]
    data["geometry"] = gpd.points_from_xy(data.x, data.y)
    data = data.set_geometry("geometry")

    data.to_file(os.path.join(output_dir, "nfplist.fgb"), driver="FlatGeobuf")


@timed
def save_large_fires_nplist(allpixels, region, large_fires, tst):
    for fid in large_fires:
        data = allpixels[allpixels["fid"] == fid]
        save_fire_nplist(data, region, fid, tst)


@timed
def save_fire_layers(allfires_gdf_fid, region, fid, tst):

    output_dir = largefire_folder(region, fid, tst, location="local")
    os.makedirs(output_dir, exist_ok=True)

    for layer in ["perimeter", "fireline", "newfirepix"]:
        columns = [col for col in singlefire_getdd(layer)]
        data = allfires_gdf_fid[columns].copy()
        if layer == "perimeter":
            data["geometry"] = allfires_gdf_fid["hull"]
        elif layer == "newfirepix":
            data["geometry"] = allfires_gdf_fid["nfp"]
        elif layer == "fireline":
            data["geometry"] = allfires_gdf_fid["fline"]

        data = data.set_geometry("geometry")
        data = data[data.geometry.notna() & ~data.geometry.is_empty]

        data.to_file(os.path.join(output_dir, f"{layer}.fgb"), driver="FlatGeobuf")


@timed
def merge_rows(allfires_gdf_fid):
    """For a subset of allfires data containing only one fire, merge any
    rows that have the same `t`
    """
    output = allfires_gdf_fid.drop_duplicates(subset=["t"]).set_index("t").copy()

    # clean up any merges that are needed
    for dt, rows in allfires_gdf_fid[allfires_gdf_fid.t.duplicated(False)].groupby("t"):
        # first get the weighted sums for pixden and meanFRP
        pixweight = (rows["pixden"] * rows["farea"]).sum()
        FRPweight = (rows["meanFRP"] * rows["n_pixels"]).sum()

        for col in ["n_pixels", "n_newpixels", "farea", "fperim", "flinelen"]:
            output.loc[dt, col] = rows[col].sum()

        output.loc[dt, "t_st"] = rows["t_st"].min()
        output.loc[dt, "pixden"] = pixweight / output.loc[dt, "farea"]
        output.loc[dt, "meanFRP"] = FRPweight / output.loc[dt, "n_pixels"]

        dissolved = rows.dissolve()
        for col in ["hull", "fline", "nfp"]:
            output.loc[dt, col] = dissolved[col].item()

    return output.reset_index()


@timed
def save_large_fires_layers(allfires_gdf, region, large_fires, tst):
    gdf = allfires_gdf.reset_index()

    merge_needed = (gdf.mergeid != gdf.fireID) & (gdf.invalid == False)
    print(f"{merge_needed.sum()} rows that potentially need a merge")

    # we'll set the "fireID" to "mergeid" in those spots
    gdf.loc[merge_needed, "fireID"] = gdf.loc[merge_needed, "mergeid"]

    for fid, data in gdf[gdf["fireID"].isin(large_fires)].groupby("fireID"):
        # merge any rows that have the same t
        if data.t.duplicated().any():
            data = merge_rows(data)

        save_fire_layers(data, region, fid, tst)


def individual_fires_path(tst, ted, region, location: Location = None):
    return os.path.join(
        all_dir(tst, region, location=location),
        f"mergedDailyFires_{ted[0]}{ted[1]:02}{ted[2]:02}_{ted[3]}.fgb"
    )


def cumunion(x):
    for i in range(1, len(x)):
        x[i] = unary_union([x[i - 1], x[i]])
    return x


@timed
def save_individual_fire(allfires_gdf, tst, ted, region):
    """save daily perimeters for an individual fire. This function ignores fids
    and key to this function is unioning the hulls from previous days, regardless of
    the original fireID/mergeid.
    """

    allfires_gdf = allfires_gdf.reset_index()

    # summarize file by t. make sure you dissolve the hull
    merged_t = (
        allfires_gdf.set_geometry("hull")
        .dissolve(
            by="t",
            aggfunc={
                "meanFRP": lambda x: (
                    allfires_gdf.loc[x.index, "meanFRP"]
                    * allfires_gdf.loc[x.index, "n_newpixels"]
                ).sum()
                / allfires_gdf.loc[x.index, "n_newpixels"].sum(),
                "n_newpixels": "sum",
                "duration": "max",
            },
        )
        .reset_index()
    )

    # calculate cumulative sum of n_newpixels
    merged_t["n_pixels"] = merged_t["n_newpixels"].cumsum()

    # combine the hulls from previous days
    merged_t["hull"] = cumunion(merged_t["hull"].tolist())

    # do the rest of the calculations
    merged_t["farea"] = merged_t["hull"].area / 10**6  # convert to km^2
    merged_t["pixden"] = merged_t["n_pixels"] / merged_t["farea"]

    # reorder columns
    col_order = [
        "t",
        "duration",
        "n_pixels",
        "n_newpixels",
        "meanFRP",
        "pixden",
        "farea",
        "hull",
    ]
    merged_t = merged_t[col_order]

    # get file output name
    output_filepath = individual_fires_path(tst, ted, region, location="local")
    # make path if necessary
    os.makedirs(os.path.dirname(output_filepath), exist_ok=True)

    # save
    merged_t.to_file(output_filepath, driver="FlatGeobuf")

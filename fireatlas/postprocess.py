import os

import datetime
from typing import Literal

import fsspec
import numpy as np
import pandas as pd
import geopandas as gpd

from shapely.ops import unary_union
import warnings

warnings.filterwarnings("ignore", "GeoSeries.notna", UserWarning)
warnings.filterwarnings("ignore", "Large object * detected in task graph", UserWarning)
warnings.filterwarnings("ignore", "Sending large graph", UserWarning)


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
    fs = fsspec.filesystem(location or settings.READ_LOCATION, use_listings_cache=False)
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
    df = pd.read_csv(filepath, index_col="uuid")
    for col in ["t", "datetime", "ext_until"]:
        df[col] = pd.to_datetime(df[col], format='ISO8601')

    return df

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

def create_snapshot_data(
        allfires_gdf,
        layer: Literal["perimeter", "fireline", "newfirepix"],
        region: Region,
        ted: datetime.datetime
):
    columns = [col for col in snapshot_getdd(layer)]
    data = allfires_gdf[[*columns, "invalid", "fireID"]].copy()

    if layer == "perimeter":
        data["geometry"] = allfires_gdf["hull"]
    if layer == "newfirepix":
        data["geometry"] = allfires_gdf["nfp"]
    if layer == "fireline":
        data["geometry"] = allfires_gdf["fline"]

    data = data.set_geometry("geometry")
    data = data[data.geometry.notna() & ~data.geometry.is_empty]

    if layer == "perimeter":
        # figure out the fire state given current t
        data["isignition"] = ted == data["t_st"]
        data["t_inactive"] = (ted - data["t_ed"]).dt.days

        data['invalid'] = data['invalid'].fillna(False)
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

    data['t'] = data['t'].dt.strftime('%Y-%m-%dT%H:%M:%S')

    # primary key is: region + fireID + 12hr slice
    data["primarykey"] = (
            data["region"] + "|" + data["fireID"].astype(str) + "|" + data['t']
    )

    # drop the columns we don't actually need
    data = data.drop(columns=["invalid"])
    return data


@timed
def save_snapshot_layers(allfires_gdf_t, region: Region, tst: TimeStep, ted: TimeStep):
    output_dir = snapshot_folder(region, tst, ted, location="local")
    os.makedirs(output_dir, exist_ok=True)

    dt = t2dt(ted)

    for layer in ["perimeter", "fireline", "newfirepix"]:
        data = create_snapshot_data(allfires_gdf_t, layer, region, dt)
        if layer == "perimeter":
            # only include perimeters what are active or may reactivate
            data = data[(data['isactive'] == 1) | (data['mayreactivate'] == 1)]
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
    # don't include fires that later get merged in.
    return last_large["mergeid"].unique()


def largefire_folder(region: Region, fid, tst: TimeStep, location: Location = None):
    return os.path.join(
        all_dir(tst, region, location=location),
        "Largefire",
        str(fid)
    )


def combined_largefire_folder(region: Region, tst: TimeStep, ted: TimeStep, location: Location = None):
    return os.path.join(
        all_dir(tst, region, location=location),
        "CombinedLargefire",
        f"{ted[0]}{ted[1]:02}{ted[2]:02}{ted[3]}",
    )


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
def save_combined_large_fire_layers(allfires_gdf, tst: TimeStep, ted: TimeStep, region: Region):
    output_dir = combined_largefire_folder(region, tst, ted, location="local")
    os.makedirs(output_dir, exist_ok=True)

    dt = t2dt(ted)
    for layer in ["perimeter", "fireline", "newfirepix"]:
        data = create_snapshot_data(allfires_gdf, layer, region, dt)
        data.to_file(os.path.join(output_dir, f"lf_{layer}.fgb"), driver="FlatGeobuf")


@timed
def fill_activefire_rows(allfires_gdf, ted):
    """For a subset of allfires data, add any rows where a fire is still
    active, but is not burning.
    """
    dd = singlefire_getdd("all")
    
    if allfires_gdf.index.names == ['fireID', 't']:
        gdf = allfires_gdf.reset_index()
    else:
        gdf = allfires_gdf

    dt = t2dt(ted)
    all_new_rows = []
    for fid, allfires_gdf_fid in gdf.groupby("fireID"):
        d = allfires_gdf_fid.set_index("t")
        
        # set last timestep so we fill after the fire has stopped burning
        if not d.iloc[-1]["invalid"]:
            last_t = d.index[-1]
            if last_t != dt:
                last_t = min(last_t + datetime.timedelta(days=settings.maxoffdays), dt)
                d.loc[last_t] = None

        ffilled = d.resample("12H").ffill(limit=settings.limoffdays*2).dropna(how="all")

        # get all the rows that are new
        new_rows = ffilled[~ffilled.index.isin(d.index)]

        # set values that should not be forward filled.
        new_rows.loc[:,["n_newpixels", "meanFRP", "nfp"]] = 0, None, None
        
        all_new_rows.append(new_rows.reset_index())
            
    output = pd.concat([gdf, *all_new_rows]).sort_values(["t", "fireID"])
    for k, tp in dd.items():
        output[k] = output[k].astype(tp)
    return output


def merge_rows(allfires_gdf_fid, fid: int | str):
    """For a subset of allfires data containing only one fire, merge any
    rows that have the same `t`
    """
    output = allfires_gdf_fid.dissolve(
        by="t",
        aggfunc={
            "meanFRP": lambda x: (
                None 
                if (n_newpixels := allfires_gdf_fid.loc[x.index, "n_newpixels"]).sum() == 0 
                else (x * n_newpixels).sum() / n_newpixels.sum()
            ),
            "n_newpixels": "sum",
            "fline": unary_union,  # still questioning if this will give the right outcome
            "nfp": unary_union,
            "t_st": "min",
            "t_ed": "max",
        },
    )
    t_diff = output["t_ed"] - output.index.min()
    output["duration"] = t_diff.dt.seconds / 24 / 3600 + t_diff.dt.days
    output["n_pixels"] = output.n_newpixels.cumsum()
    output["farea"] = output.hull.area / 1e6  # km2
    output["fperim"] = output.hull.length / 1e3  # km
    output["flinelen"] = output.fline.apply(lambda x: 0 if x is None else x.length) / 1e3  # km
    output["pixden"] = output.n_pixels / output.farea
    output["fireID"] = fid
    output["mergeid"] = fid
    output["invalid"] = False
    output["ftype"] = allfires_gdf_fid.ftype.mode()[0]
    return output.reset_index()


@timed
def save_large_fires_layers(allfires_gdf, region, large_fires, tst, ted, client=None):
    """will save individual large fire artifacts as well as

    a combined perimeter, newpixel and fireline artifact of all large fires
    """
    gdf = allfires_gdf.reset_index()

    gdf = gdf[gdf["fireID"].isin(large_fires) | gdf["mergeid"].isin(large_fires)]
    
    # forward fill any timesteps that are mising
    gdf = fill_activefire_rows(gdf, ted)
    
    merge_needed = (gdf.mergeid != gdf.fireID) & (gdf.invalid == False)
    print(f"{merge_needed.sum()} rows that potentially need a merge")

    # we'll set the "fireID" to "mergeid" in those spots
    gdf.loc[merge_needed, "fireID"] = gdf.loc[merge_needed, "mergeid"]

    def merge_and_save_fire(data, fid):
        # merge any rows that have the same t
        if data.t.duplicated().any():
            data = merge_rows(data, fid)

        # save off single large fire artifacts
        save_fire_layers(data, region, int(fid), tst)
        
        # accumulate each fid for combined large fires
        return data

    futures = []
    processed_gdfs = []
    for fid, data in gdf[gdf["fireID"].isin(large_fires)].groupby("fireID"):
        if client:
            futures.append(client.submit(merge_and_save_fire, data, fid))
        else:
            processed_gdfs.append(merge_and_save_fire(data, fid))
    if futures:
        processed_gdfs = client.gather(futures)

    # save off all large fire artifacts
    if len(processed_gdfs) != 0:
        all_gdfs = gpd.GeoDataFrame(pd.concat(processed_gdfs, ignore_index=True))
        save_combined_large_fire_layers(all_gdfs, tst, ted, region)


@timed
def save_individual_fire(allfires_gdf, tst, ted, region):
    """save daily perimeters for an individual fire. This function ignores fids
    and key to this function is unioning the hulls from previous days, regardless of
    the original fireID/mergeid.
    """
    gdf = allfires_gdf.reset_index()

    # forward fill any timesteps that are mising
    gdf = fill_activefire_rows(gdf, ted)

    # merge any rows that need merging across t
    data = merge_rows(gdf, fid=region[0])

    # save fire layers - use region name as fid
    save_fire_layers(data, region, region[0], tst)

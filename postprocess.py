import os

import geopandas as gpd
import pandas as pd
from FireTypes import Region, TimeStep
import FireConsts
from utils import timed


def allpixels_filepath(t: TimeStep, region: Region):
    filename = f"allpixels_{t[0]}{t[1]:02}{t[2]:02}_{t[3]}.csv"
    return os.path.join(FireConsts.diroutdata, region[0], filename)


@timed
def save_allpixels(allpixels, ted: TimeStep, region: Region):
    output_filepath = allpixels_filepath(ted, region)

    # make path if necessary
    os.makedirs(os.path.dirname(output_filepath), exist_ok=True)

    allpixels.to_csv(output_filepath)


@timed
def read_allpixels(ted: TimeStep, region: Region):
    filepath = allpixels_filepath(ted, region)

    return pd.read_csv(filepath, index_col="uuid", parse_dates=["t"])


def allfires_filepath(t: TimeStep, region: Region):
    filename = f"allfires_{t[0]}{t[1]:02}{t[2]:02}_{t[3]}.parq"
    return os.path.join(FireConsts.diroutdata, region[0], filename)


@timed
def save_allfires_gdf(allfires_gdf, ted: TimeStep, region: Region):
    output_filepath = allfires_filepath(ted, region)

    # make path if necessary
    os.makedirs(os.path.dirname(output_filepath), exist_ok=True)

    allfires_gdf.to_parquet(output_filepath)


@timed
def read_allfires_gdf(ted: TimeStep, region: Region):
    filepath = allfires_filepath(ted, region)

    return gpd.read_parquet(filepath)


@timed
def find_largefires(allfires_gdf):
    gdf = allfires_gdf.reset_index()

    last_seen = gdf.drop_duplicates("fireID", keep="last")
    last_large = last_seen[
        (last_seen.farea > FireConsts.LARGEFIRE_FAREA) & (last_seen.invalid == False)
    ]
    return last_large.fireID.values


def fire_folder(region: Region, fid):
    return os.path.join(FireConsts.diroutdata, region[0], "fires", str(fid))


@timed
def save_fire_nplist(allpixels_fid, region, fid):

    output_dir = fire_folder(region, fid)
    os.makedirs(output_dir, exist_ok=True)

    subset = allpixels_fid[["x", "y", "FRP", "DS", "DT", "ampm", 'YYYYMMDD_HHMM', "Sat"]].copy()
    subset.columns = ["x", "y", "frp", "DS", "DT", "ampm", 'datetime', "sat"]
    subset["geometry"] = gpd.points_from_xy(subset.x, subset.y)
    subset = subset.set_geometry("geometry")
    
    subset.to_file(os.path.join(output_dir, "nfplist.fgb"), driver="FlatGeobuf")


@timed
def save_large_fires_nplist(allpixels, region, large_fires):
    for fid in large_fires:
        data = allpixels[allpixels["fid"] == fid]
        save_fire_nplist(data, region, fid)


@timed
def save_fire_layers(allfires_gdf_fid, region, fid):
    from FireGpkg_sfs import getdd

    output_dir = fire_folder(region, fid)
    os.makedirs(output_dir, exist_ok=True)

    for layer in ["perimeter", "fireline", "newfirepix"]:
        columns = [col for col in getdd(layer)]
        subset = allfires_gdf_fid[columns].copy()
        if layer == "perimeter":
            subset["geometry"] = allfires_gdf_fid["hull"]
        elif layer == "newfirepix":
            subset["geometry"] = allfires_gdf_fid["nfp"]
        elif layer == "fireline":
            subset["geometry"] = allfires_gdf_fid["fline"]
            subset = subset.dropna(subset=["geometry"])
        subset = subset.set_geometry("geometry")
        
        subset.to_file(os.path.join(output_dir, f"{layer}.fgb"), driver="FlatGeobuf")


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
def save_large_fires_layers(allfires_gdf, region, large_fires):
    gdf = allfires_gdf.reset_index()

    merge_needed = (gdf.mergeid != gdf.fireID) & (gdf.invalid == False)
    print(f"{merge_needed.sum()} rows that potentially need a merge")

    # we'll set the "fireID" to "mergeid" in those spots
    gdf.loc[merge_needed, "fireID"] = gdf.loc[merge_needed, "mergeid"]

    for fid, data in gdf[gdf["fireID"].isin(large_fires)].groupby("fireID"):
        # merge any rows that have the same t
        if data.t.duplicated().any():
            data = merge_rows(data)
        
        save_fire_layers(data, region, fid)

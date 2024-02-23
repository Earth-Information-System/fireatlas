import os

import numpy as np
import pandas as pd
import geopandas as gpd

from FireTypes import Region, TimeStep
import FireConsts, FireTime
from utils import timed

import warnings; warnings.filterwarnings('ignore', 'GeoSeries.notna', UserWarning)

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


def snapshot_folder(region: Region, t: TimeStep):
    return os.path.join(FireConsts.diroutdata, region[0], "Snapshot", f"{t[0]}{t[1]:02}{t[2]:02}{t[3]}")


@timed
def save_snapshot_layers(allfires_gdf_t, region: Region, t: TimeStep):
    from FireGpkg import getdd

    output_dir = snapshot_folder(region, t)
    os.makedirs(output_dir, exist_ok=True)

    dt = FireTime.t2dt(t)

    for layer in ["perimeter", "fireline", "newfirepix"]:
        columns = [col for col in getdd(layer)]
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
            
            data["isactive"] = ~data["invalid"] & (data["t_inactive"] <= FireConsts.maxoffdays)
            data["isdead"] = ~data["invalid"] & (data["t_inactive"] > FireConsts.limoffdays)
            data["mayreactivate"] = ~data["invalid"] & (FireConsts.maxoffdays < data["t_inactive"]) & (data["t_inactive"] <= FireConsts.limoffdays)

            # map booleans to integers
            for col in ["isignition", "isactive", "isdead", "mayreactivate"]:
                data[col] = data[col].astype(int)

            # apply filter flag
            data['geom_counts'] = data[["fireID", "geometry"]].explode(index_parts=True).groupby(['fireID']).nunique()["geometry"] # count number of polygons
            data['low_confidence_grouping'] = np.where(data['geom_counts']>5, 1, 0) # if more than 5 geometries are present, flag it

        
        data["region"] = str(region[0])
        
        # primary key is: region + fireID + 12hr slice 
        data['primarykey'] = data['region'] + '|' + data["fireID"].astype(str) + '|' + dt.isoformat()

        # drop the columns we don't actually need
        data = data.drop(columns=["invalid", "fireID"])

        data.to_file(os.path.join(output_dir, f"{layer}.fgb"), driver="FlatGeobuf")



@timed
def save_snapshots(allfires_gdf, region, tst, ted):
    gdf = allfires_gdf.reset_index()

    for t in FireTime.t_generator(tst, ted):
        dt = FireTime.t2dt(t)
        data = gdf[gdf.t <= dt].drop_duplicates("fireID", keep="last")
        save_snapshot_layers(data, region, t)


@timed
def find_largefires(allfires_gdf):
    gdf = allfires_gdf.reset_index()

    last_seen = gdf.drop_duplicates("fireID", keep="last")
    last_large = last_seen[
        (last_seen.farea > FireConsts.LARGEFIRE_FAREA) & (last_seen.invalid == False)
    ]
    return last_large.fireID.values


def fire_folder(region: Region, fid):
    return os.path.join(FireConsts.diroutdata, region[0], "Largefire", str(fid))


@timed
def save_fire_nplist(allpixels_fid, region, fid):

    output_dir = fire_folder(region, fid)
    os.makedirs(output_dir, exist_ok=True)

    data = allpixels_fid[["x", "y", "FRP", "DS", "DT", "ampm", 'YYYYMMDD_HHMM', "Sat"]].copy()
    data.columns = ["x", "y", "frp", "DS", "DT", "ampm", 'datetime', "sat"]
    data["geometry"] = gpd.points_from_xy(data.x, data.y)
    data = data.set_geometry("geometry")
    
    data.to_file(os.path.join(output_dir, "nfplist.fgb"), driver="FlatGeobuf")


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

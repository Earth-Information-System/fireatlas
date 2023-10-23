import argparse
import pandas as pd
import geopandas as gpd
import numpy as np
import s3fs
import boto3
import re
import time
from pathlib import Path
from FireConsts import diroutdata
from dask.distributed import Client
from FireLog import logger
from FireIO import copy_from_maap_to_veda_s3

LAYERS = ["newfirepix", "fireline", "perimeter"]
MAX_WORKERS = 14
LOCAL_DIR_OUTPUT_PREFIX_PATH = "/tmp/{}/LargeFire_Outputs"


def mkdir_dash_p(parent_output_path):
    """named after linux bash `mkdir -p`

    this function will create all parent folders
    for a path if they don't already exist and
    if they do exist, it will gracefully ignore them

    Examples:
        input: `mkdir_dash_p('/tmp/foo/bar/doobie/README.txt')`
        output: all nested parent directories of the file are created "/tmp/foo/bar/doobie"

    :param parent_output_path:
    :return:
    """
    path = Path(parent_output_path)
    path.parent.mkdir(parents=True, exist_ok=True)


def merge_lf_years(
    parent_years_folder_input_path,
    maap_output_folder_path,
    folder_name,
    layers=LAYERS,
):
    """using `glob` and `concat` merge all large fire layers across years and write back to MAAP s3

    :param years_range: a list of year ints
    :param parent_folder_path: local dir path to the folder that houses all output years from `combine_per
    :param layers: a list of layer strings
    :param folder_name: alias for region
    :return:
    """
    for layer in layers:
        folder = Path(parent_years_folder_input_path)
        logger.info(f"[ PARENT ]: years folder path: {folder}")
        flatgeobufs_by_layer_and_year = list(folder.glob(f"*/lf_{layer}.fgb"))
        logger.info(f"[ CHILD ]: fgb(s) by year: {flatgeobufs_by_layer_and_year}")
        gpd_by_year = [gpd.read_file(fgb) for fgb in flatgeobufs_by_layer_and_year]
        logger.info(f"[ GPD ]: frames by year: {gpd_by_year}")
        gdf = pd.concat(gpd_by_year).pipe(gpd.GeoDataFrame)

        maap_s3_layer_path = f"{maap_output_folder_path}/lf_{layer}_archive.fgb"
        if IS_NRT_RUN:
            maap_s3_layer_path = f"{maap_output_folder_path}/lf_{layer}_nrt.fgb"
        gdf.to_file(
            maap_s3_layer_path,
            driver="FlatGeobuf",
        )
        if IS_PRODUCTION_RUN:
            copy_from_maap_to_veda_s3(maap_s3_layer_path, folder_name)


def load_lf(lf_id, file_path, layer="nfplist", drop_duplicate_geometries=False):
    """find the large fire pickled file by id

    :param lf_id: integer
    :param file_path: s3 MAAP path to inputs
    :param layer: str
    :param drop_duplicate_geometries: bool
    :return: pandas.DataFrame
    """
    try:
        logger.info(f"[ READ FILE ]: {file_path}/{layer}")
        gdf = gpd.read_file(file_path, layer=layer)
    except Exception as e:
        logger.exception(e)
        return
    gdf["ID"] = lf_id

    if (drop_duplicate_geometries == True) and (layer != "nfplist"):
        gdf.sort_values("t", ascending=True, inplace=True)
        gdf = gdf.loc[gdf["geometry"].drop_duplicates(keep="first").index]
    return gdf


def write_lf_layers_by_year(
    year, 
    s3_maap_input_path, 
    LOCAL_DIR_OUTPUT_PREFIX_PATH,
    layers=LAYERS
):
    """ for each layer write out the most recent lf layer

    :param year: int
    :param s3_maap_input_path: the s3 MAAP path
    :param LOCAL_DIR_OUTPUT_PREFIX_PATH: the local dir output prefix
    :param layers: a list of layer strings
    :return: None
    """
    s3 = s3fs.S3FileSystem(anon=False)

    # load in NRT Largefire data
    lf_files = [f for f in s3.ls(s3_maap_input_path)]

    lf_files.sort()
    lf_ids = list(
        set([file.split("Largefire/")[1].split("_")[0] for file in lf_files])
    )  # unique lf ids

    # each largefire id has a file for each timestep which has entire evolution up to that point.
    # the latest one has the most up-to-date info for that fire
    largefire_dict = dict.fromkeys(lf_ids)
    for lf_id in lf_ids:
        most_recent_file = (
            "s3://" + [file for file in lf_files if lf_id in file][-1]
        )  # most recent file is last!
        largefire_dict[lf_id] = most_recent_file

    # NOTE: there's another opportunity to speed up code here for each year process
    # if we spin up a couple more dask workers per year and scatter each layer reads/writes
    for layer in layers:
        all_lf_per_layer = pd.concat(
            [
                load_lf(lf_id, file_path, layer=f"{layer}")
                for lf_id, file_path in largefire_dict.items()
            ],
            ignore_index=True,
        )

        #### APPLY GEOMETRY FLAG ####

        all_lf_per_layer.set_index(['ID','t'],inplace=True) # set index for each day's observation of each fire ID
        multip = all_lf_per_layer[all_lf_per_layer['geometry'].geom_type=='MultiPolygon'].index # get index of observations with multipolygons
        all_lf_per_layer.loc[multip,'max_geom_counts'] = all_lf_per_layer.loc[multip]['geometry'].explode(index_parts=True).groupby(['ID','t']).nunique() # count num polygons
        all_lf_per_layer['max_geom_counts'].fillna(1,inplace=True) # fill the rest in as 1s bc only one polygon geometry
        all_lf_per_layer['low_confidence_grouping'] = np.where(all_lf_per_layer['max_geom_counts']>5, 1, 0) # set the flag
        all_lf_per_layer.reset_index(inplace=True) # reset index for saving

        #############################

        layer_output_fgb_path = f"{LOCAL_DIR_OUTPUT_PREFIX_PATH}/{year}/lf_{layer}.fgb"
        # create all parent directories for local output (if they don't exist already)
        mkdir_dash_p(layer_output_fgb_path)
        all_lf_per_layer.to_file(
            layer_output_fgb_path,
            driver="FlatGeobuf",
        )


def main(
    years_range: list, 
    folder_name: str, 
    in_parallel: bool = False
):
    """
    :param years_range: a list of year integers
    :param in_parallel: bool
    :return: None
    """

    if not in_parallel:
        for year in years_range:
            s3_maap_input_path = f"{diroutdata}{folder_name}/{year}/Largefire/"
            write_lf_layers_by_year(
                year, 
                s3_maap_input_path, 
                LOCAL_DIR_OUTPUT_PREFIX_PATH.format(folder_name)
            )
        merge_lf_years(
            LOCAL_DIR_OUTPUT_PREFIX_PATH,
            f"{diroutdata}{LOCAL_DIR_OUTPUT_PREFIX_PATH.format(folder_name)}/LargeFire_Outputs/merged",
            folder_name
        )
        return

    #############################################
    # in parallel via dask
    #############################################
    # this assumes we are running on our largest worker instance where we have 16 CPU
    # so we limit (using modulo) to 14 workers at most but default to using `len(years_range)` workers

    dask_client = Client(n_workers=(len(years_range) % MAX_WORKERS))
    logger.info(f"workers = {len(dask_client.cluster.workers)}")
    tstart = time.time()
    # set up work items
    futures = [
        dask_client.submit(
            write_lf_layers_by_year,
            year,
            f"{diroutdata}{folder_name}/{year}/Largefire/",
            LOCAL_DIR_OUTPUT_PREFIX_PATH,
        )
        for year in years_range
    ]
    # join children and wait to finish
    dask_client.gather(futures)
    logger.info(
        f"workers after dask_client.gather = {len(dask_client.cluster.workers)}"
    )
    tend = time.time()
    logger.info(f'"write_lf_layers_by_year" in parallel: {(tend - tstart) / 60} minutes')
    dask_client.restart()
    merge_lf_years(
        LOCAL_DIR_OUTPUT_PREFIX_PATH,
        f"{diroutdata}{folder_name}/LargeFire_Outputs/merged",
        folder_name
    )


if __name__ == "__main__":
    """
    Examples:
        # run a single year
        python3 combine_largefire.py -s 2023 -e 2023 
        
        # run multiple years in parallel
        python3 combine_largefire.py -s 2018 -e 2022 -p
        
        # run multiple years in parallel in PRODUCTION mode, by default this is the LF archive (not NRT)
        # PRODUCTION mode means outputs get `aws s3 cp` to s3://veda-data-store-staging
        python3 combine_largefire.py -s 2018 -e 2022 -p -x
        
        # run NRT current year in PRODUCTION mode
        # PRODUCTION mode means outputs get `aws s3 cp` to s3://veda-data-store-staging
        python3 combine_largefire.py -s 2023 -e 20203 -x --nrt
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-s", "--start-year", required=True, type=int, help="start year int"
    )
    parser.add_argument(
        "-e", "--end-year", required=True, type=int, help="end year int"
    )
    parser.add_argument(
        "-p",
        "--parallel",
        action="store_true",
        help="turn on dask processing years in parallel",
    )
    parser.add_argument(
        "-x",
        "--production-run",
        action="store_true",
        help="creates a flag/trap for us to know this is running as the PRODUCTION job to avoid overwrite VEDA s3 data",
    )
    parser.add_argument(
        "-n",
        "--nrt",
        action="store_true",
        help="creates a flag/trap to know if this is LF archive (all years) or LF NRT (current year)",
    )
    # parser argument for maap s3 bucket folder name mapped to folder_name
    parser.add_argument(
        "--folder-name",
        required=True,
        type=str,
        default="CONUS_NRT_DPS",
        help="maap s3 bucket folder name",
    )

    args = parser.parse_args()

    # set global flag/trap to protect VEDA s3 copy
    global IS_PRODUCTION_RUN
    IS_PRODUCTION_RUN = args.production_run

    # set global flag/trap for NRT
    global IS_NRT_RUN
    IS_NRT_RUN = args.nrt

    # validate `years_range` construction
    years_range = list(range(args.start_year, args.end_year + 1))
    if years_range[0] != args.start_year or years_range[-1] != args.end_year:
        raise ValueError(
            f"[ ERROR ]: the range='{years_range}' doesn't start or end with inputs='{args.start_year}-{args.end_year}'"
        )

    start = time.time()
    logger.info(f"Running algo with year range: '{years_range}'")
    main(years_range, args.folder_name, in_parallel=args.parallel)
    total = time.time() - start
    logger.info(f"Total runtime is {str(total / 60)} minutes")

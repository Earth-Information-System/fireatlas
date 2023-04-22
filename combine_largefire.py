import argparse
import pandas as pd
import geopandas as gpd
import s3fs
import boto3
import re
import time
from pathlib import Path
from FireConsts import diroutdata
from dask.distributed import Client
from FireLog import logger


def mkdir_dash_p(parent_output_path):
    """named after linux bash `mkdir -p`

    this function will create all parent folders
    for a path if they don't already exist and
    if they do exist, it will gracefully leave then alone

    Examples:
        input: /tmp/foo/bar/doobie/README.txt
        output: all nested parent directories of the file are created "/tmp/foo/bar/doobie"

    :param parent_output_path:
    :return:
    """
    path = Path(parent_output_path)
    path.parent.mkdir(parents=True, exist_ok=True)


def copy_from_maap_to_veda_s3(from_maap_s3_path):
    s3_client = boto3.client("s3")

    if "LargeFire" in from_maap_s3_path:
        try:
            fname_regex = r"^s3://maap.*?(/LargeFire_Outputs/)merged/(?P<fname>lf_fireline.fgb|lf_perimeter.fgb|lf_newfirepix.fgb|nfplist.fgb)$"
            # note that `destination_dict` should resemble this output with a match if the URL was a perimeter file:
            # {'fname': 'lf_perimeter.fgb'}
            destination_dict = (
                re.compile(fname_regex).match(from_maap_s3_path).groupdict()
            )
        except AttributeError:
            logger.error(f"[ NO REGEX MATCH FOUND ]: for file f{from_maap_s3_path}")
            return

        from_maap_s3_path = from_maap_s3_path.replace("s3://", "")
        s3_client.copy_object(
            CopySource=from_maap_s3_path,  # full bucket path
            Bucket="veda-data-store-staging",  # Destination bucket
            Key=f"EIS/FEDSoutput/LFArchive/{destination_dict['']}/{destination_dict['fname']}"
            # Destination path/filename
        )
    else:
        logger.error(f"[ NO S3 COPY EXPORTED ]: for file f{from_maap_s3_path}")


def merge_df_years(
    parent_years_folder_input_path,
    maap_output_folder_path,
    layers=["nlplist.fgb", "newfirepix.fgb", "fireline.fgb", "perimeter.fgb"],
):
    """
    :param years_range:
    :param parent_folder_path:
    :param layers:
    :return:
    """
    for layer in layers:
        folder = Path(parent_years_folder_input_path)
        logger.info(f"[ PARENT ]: years folder path: {folder}")
        flatgeobufs_by_layer_and_year = list(folder.glob(f"*/lf_{layer}"))
        logger.info(f"[ CHILD ]: fgb(s) by year: {flatgeobufs_by_layer_and_year}")
        gpd_by_year = [gpd.read_file(fgb) for fgb in flatgeobufs_by_layer_and_year]
        logger.info(f"[ GPD ]: frames by year: {gpd_by_year}")
        gdf = pd.concat(gpd_by_year).pipe(gpd.GeoDataFrame)

        maap_s3_layer_path = f"{maap_output_folder_path}/lf_{layer}"
        mkdir_dash_p(maap_s3_layer_path)
        gdf.to_file(
            maap_s3_layer_path,
            driver="FlatGeobuf",
        )
        copy_from_maap_to_veda_s3(maap_s3_layer_path)


def load_lf(lf_id, file_path, layer="nfplist", drop_duplicate_geometries=False):
    """
    :param lf_id: integer
    :param file_path:
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


def combine_by_year(year, s3_maap_input_path, local_dir_output_prefix_path):
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

    all_lf_nfplist = pd.concat(
        [
            load_lf(lf_id, file_path, layer="nfplist")
            for lf_id, file_path in largefire_dict.items()
        ],
        ignore_index=True,
    )
    all_lf_firelines = pd.concat(
        [
            load_lf(lf_id, file_path, layer="fireline")
            for lf_id, file_path in largefire_dict.items()
        ],
        ignore_index=True,
    )
    all_lf_perimeters = pd.concat(
        [
            load_lf(lf_id, file_path, layer="perimeter")
            for lf_id, file_path in largefire_dict.items()
        ],
        ignore_index=True,
    )
    all_lf_newfirepix = pd.concat(
        [
            load_lf(lf_id, file_path, layer="newfirepix")
            for lf_id, file_path in largefire_dict.items()
        ],
        ignore_index=True,
    )

    nfplist_maap_fgb_path = f"{local_dir_output_prefix_path}/{year}/lf_nfplist.fgb"
    # create all parent directories for local output (if they don't exist already)
    # we only need to do this once for the first layer since subsequent
    # layers will have the same parent dirs
    mkdir_dash_p(nfplist_maap_fgb_path)
    all_lf_nfplist.to_file(
        nfplist_maap_fgb_path,
        driver="FlatGeobuf",
    )

    fireline_maap_fgb_path = f"{local_dir_output_prefix_path}/{year}/lf_fireline.fgb"
    all_lf_firelines.to_file(
        fireline_maap_fgb_path,
        driver="FlatGeobuf",
    )

    perimeter_maap_fgb_path = f"{local_dir_output_prefix_path}/{year}/lf_perimeter.fgb"
    all_lf_perimeters.to_file(
        perimeter_maap_fgb_path,
        driver="FlatGeobuf",
    )

    newfirepix_maap_fgb_path = (
        f"{local_dir_output_prefix_path}/{year}/lf_newfirepix.fgb"
    )
    all_lf_newfirepix.to_file(
        newfirepix_maap_fgb_path,
        driver="FlatGeobuf",
    )


def main(years_range, in_parallel=False):
    """
    :param years_range: a list of year integers
    :param in_parallel: bool
    :return: None
    """
    local_dir_output_prefix_path = "/tmp/CONUS_NRT_DPS/LargeFire_Outputs"
    if not in_parallel:
        for year in years_range:
            s3_maap_input_path = f"{diroutdata}CONUS_NRT_DPS/{year}/Largefire/"
            combine_by_year(year, s3_maap_input_path, local_dir_output_prefix_path)
        merge_df_years(
            local_dir_output_prefix_path,
            f"{diroutdata}CONUS_NRT_DPS/LargeFire_Outputs/merged",
        )
        return

    #############################################
    # in parallel via dask
    #############################################
    # this assumes we are running on our largest worker instance where we have 16 CPU
    # so we limit (using modulo) to 14 workers at most but default to using `len(years_range)` workers
    max_workers = 14
    dask_client = Client(n_workers=(len(years_range) % max_workers))
    logger.info(f"workers = {len(dask_client.cluster.workers)}")
    tstart = time.time()
    # set up work items
    futures = [
        dask_client.submit(
            combine_by_year,
            year,
            f"{diroutdata}CONUS_NRT_DPS/{year}/Largefire/",
            local_dir_output_prefix_path,
        )
        for year in years_range
    ]
    # join children and wait to finish
    dask_client.gather(futures)
    logger.info(
        f"workers after dask_client.gather = {len(dask_client.cluster.workers)}"
    )
    tend = time.time()
    logger.info(f'"combine_by_year" in parallel: {(tend - tstart) / 60} minutes')
    dask_client.restart()
    merge_df_years(
        local_dir_output_prefix_path,
        f"{diroutdata}CONUS_NRT_DPS/LargeFire_Outputs/merged",
    )


if __name__ == "__main__":
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
    args = parser.parse_args()

    # validate `years_range` construction
    years_range = list(range(args.start_year, args.end_year + 1))
    if years_range[0] != args.start_year or years_range[-1] != args.end_year:
        raise ValueError(
            f"[ ERROR ]: the range='{years_range}' doesn't start or end with inputs='{args.start_year}-{args.end_year}'"
        )

    start = time.time()
    logger.info(f"Running algo with year range: '{years_range}'")
    main(years_range, in_parallel=args.parallel)
    total = time.time() - start
    logger.info(f"Total runtime is {str(total / 60)} minutes")

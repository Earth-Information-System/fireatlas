import argparse
import pandas as pd
import geopandas as gpd
import s3fs
import time
from FireConsts import diroutdata
from dask.distributed import Client
from FireLog import logger


def load_lf(lf_id, file_path, layer="nfplist", drop_duplicate_geometries=False):
    """
    :param lf_id: integer
    :param file_path:
    :param layer: str
    :param drop_duplicate_geometries: bool
    :return: pandas.DataFrame
    """
    try:
        logger.info(f"[ READ FILE ]: filepath={file_path} for layer={layer}")
        gdf = gpd.read_file(file_path, layer=layer)
    except Exception as e:
        logger.exception(e)
        return
    gdf["ID"] = lf_id

    if (drop_duplicate_geometries == True) and (layer != "nfplist"):
        gdf.sort_values("t", ascending=True, inplace=True)
        gdf = gdf.loc[gdf["geometry"].drop_duplicates(keep="first").index]
    return gdf


def combine_by_year(year, s3_maap_input_path, s3_maap_output_prefix_path):
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

    all_lf_nfplist.to_file(
        f"{s3_maap_output_prefix_path}/{year}/lf_nfplist_{year}.fgb",
        driver="FlatGeobuf",
    )
    all_lf_firelines.to_file(
        f"{s3_maap_output_prefix_path}/{year}/lf_firelines_{year}.fgb",
        driver="FlatGeobuf",
    )
    all_lf_perimeters.to_file(
        f"{s3_maap_output_prefix_path}/{year}/lf_perimeters_{year}.fgb",
        driver="FlatGeobuf",
    )
    all_lf_newfirepix.to_file(
        f"{s3_maap_output_prefix_path}/{year}/lf_newfirepix_{year}.fgb",
        driver="FlatGeobuf",
    )


def main(years_range, in_parallel=False):
    """
    :param years_range: a list of year integers
    :param in_parallel: bool
    :return: None
    """
    if not in_parallel:
        for year in years_range:
            s3_maap_input_path = f"{diroutdata}CONUS_NRT_DPS/{year}/Largefire/"
            s3_maap_output_prefix_path = f"{diroutdata}CONUS_NRT_DPS/LargeFire_Outputs"
            combine_by_year(year, s3_maap_input_path, s3_maap_output_prefix_path)
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
            f"{diroutdata}CONUS_NRT_DPS/LargeFire_Outputs",
        )
        for year in years_range
    ]
    # join children and wait to finish
    dask_client.gather(futures)
    logger.info(
        f"workers after dask_client.gather = {len(dask_client.cluster.workers)}"
    )
    tend = time.time()
    logger.info(f'"combine_by_year" in parallel: {(tend - tstart) / 60.} minutes')
    dask_client.restart()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Combine Largefire Archive")
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
    logger.info("Total runtime is", str(total / 60), "minutes")

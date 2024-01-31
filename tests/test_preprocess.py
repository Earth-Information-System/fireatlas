import os
import FireConsts
import FireClustering
import FireMain
import FireIO
import pathlib
import preprocess
import pytest
import pandas as pd
from datetime import datetime
from FireTypes import TimeStep
from shapely.geometry import Polygon
from unittest.mock import MagicMock

try:
    from shapely import to_geojson, from_geojson
except ImportError:
    from shapely.geometry import mapping as to_geojson, shape as from_geojson

from unittest.mock import call


def test_preprocess_region(tmpdir, monkeypatch):
    # arrange
    monkeypatch.setattr(FireConsts, "dirprpdata", tmpdir)
    monkeypatch.setattr(preprocess, "OUTPUT_DIR", tmpdir)
    test_region = ["Test123", Polygon([(0, 0), (0, 1), (1, 1), (1, 0), (0, 0)])]
    monkeypatch.setattr(FireMain, "maybe_remove_static_sources", lambda x, y: x)

    # act
    preprocess.preprocess_region(test_region)

    # assert
    with open(f"{tmpdir}{test_region[0]}.json", "r") as f:
        assert test_region[1] == from_geojson(f.read())


def test_read_region(tmpdir, monkeypatch):
    # arrange
    monkeypatch.setattr(FireConsts, "dirprpdata", tmpdir)
    monkeypatch.setattr(preprocess, "OUTPUT_DIR", tmpdir)
    expected_region = ("Test123", Polygon([(0, 0), (0, 1), (1, 1), (1, 0), (0, 0)]))
    with open(f"{tmpdir}{expected_region[0]}.json", "w") as f:
        f.write(to_geojson(expected_region[1]))

    # act
    actual_region = preprocess.read_region(expected_region)

    # assert
    assert expected_region == actual_region


def test_preprocess_landcover(tmpdir, mock_rasterio, monkeypatch):
    # arrange
    monkeypatch.setattr(FireConsts, "dirextdata", tmpdir)
    monkeypatch.setattr(preprocess, "INPUT_DIR", tmpdir)

    tmpdir = tmpdir.join("NLCD")
    tmpdir.mkdir()

    file_basename = "nlcd_export_510m_simplified"
    file_path = os.path.join(FireConsts.dirextdata, "NLCD", f"{file_basename}.tif")
    with open(file_path, "w") as f:
        f.write("")

    # act
    preprocess.preprocess_landcover()

    # assert
    actual_first_call_args = mock_rasterio.open.call_args_list[0]
    assert actual_first_call_args == call(file_path)

    actual_second_call_args = mock_rasterio.open.call_args_list[1]
    assert actual_second_call_args == call(
        os.path.join(FireConsts.dirextdata, "NLCD", f"{file_basename}_latlon.tif"), "w"
    )

    mock_rasterio.warp.calculate_default_transform.assert_called()

    # TODO: we need a fixture that creates a TIFF
    # mock_rasterio.warp.reproject.assert_called()


@pytest.mark.parametrize(
    "timestep, sat",
    [
        ((2023, 11, 9, "AM"), "NOAA20"),
        (
            (2023, 11, 9, "AM"),
            "SNPP",
        ),
    ],
)
def test_preprocess_NRT_file(timestep: TimeStep, sat: str, monkeypatch, test_data_dir):
    monkeypatch.setattr(FireConsts, "dirextdata", test_data_dir)
    monkeypatch.setattr(preprocess, "INPUT_DIR", test_data_dir)

    if sat == "SNPP":
        df_filtered_paths = preprocess.preprocess_NRT_file(timestep, sat)
    else:
        df_filtered_paths = preprocess.preprocess_NRT_file(timestep, sat)
    assert len(df_filtered_paths) == 2

    # TODO: more assertions on the filtered CSVs


@pytest.mark.parametrize(
    "region, region_shape_to_filter, output_should_already_exist",
    [
        (
            ["Testing123", Polygon([(0, 0), (0, 2), (2, 2), (2, 0), (0, 0)])],
            Polygon([(0, 0), (0, 1), (1, 1), (1, 0), (0, 0)]),
            False,
        ),  # preprocessed_nrt_snpp_tmpfile fixture should only have one pixel that fits in this region shape
        (
            ["Testing123", Polygon([(0, 0), (0, 2), (2, 2), (2, 0), (0, 0)])],
            Polygon([(0, 0), (0, 3), (3, 3), (3, 0), (0, 0)]),
            True,
        ),  # preprocessed_nrt_snpp_tmpfile fixture should have two pixels that fix in this region shape
    ],
)
def test_preprocess_region_t(
    inputdirs,
    preprocessed_nrt_snpp_tmpfile,
    monkeypatch,
    region,
    region_shape_to_filter,
    output_should_already_exist,
):
    # arrange
    monkeypatch.setattr(FireConsts, "dirextdata", str(inputdirs))
    monkeypatch.setattr(FireConsts, "dirprpdata", str(inputdirs / "processed"))

    # return a preprocessed SNPP file fixture with a couple pixels
    input_df = pd.read_csv(preprocessed_nrt_snpp_tmpfile)

    monkeypatch.setattr(preprocess, "read_preprocessed", lambda x, sat=None: input_df)

    monkeypatch.setattr(FireIO, "get_reg_shp", lambda x: region_shape_to_filter)

    # we should have separate unit tests that really exercise `FireClustering.do_clustering`
    # and `FireClustering.compute_all_spatial_distances` so just return the input DataFrame for now
    monkeypatch.setattr(FireClustering, "do_clustering", lambda x, y: input_df)

    # TODO: handle output already exists branching
    # monkeypatch.setattr(
    #     preprocess, "preprocessed_filename", lambda x, y, region=None: MagicMock()
    # )
    # monkeypatch.setattr(os, "makedirs", lambda x, exist_ok=None: MagicMock())

    # act
    # NOTE: we pass a fake sensor here so we can easily mock the else branch of and handle a single `read_preprocessed`
    outfile_df_path = preprocess.preprocess_region_t(
        (2023, 11, 9, "AM"), "TESTING123", region
    )

    # assert
    assert os.path.exists(outfile_df_path)
    assert len(pd.read_csv(outfile_df_path)) > 0

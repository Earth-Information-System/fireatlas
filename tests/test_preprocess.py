import os
import pytest
import pandas as pd
from shapely.geometry import Polygon
from unittest.mock import MagicMock
from pathlib import Path

from fireatlas import preprocess
from fireatlas import FireClustering
from fireatlas import FireMain
from fireatlas import FireIO
from fireatlas.FireTypes import TimeStep, Region
from fireatlas import settings

try:
    from shapely import to_geojson, from_geojson
except ImportError:
    from shapely.geometry import mapping as to_geojson, shape as from_geojson

from unittest.mock import call


@pytest.mark.parametrize(
    "region, location, new_s3_path, expected_filepath",
    [
        (
            ["TestDefault", None],
            "s3",
            None,
            f"{settings.S3_PATH}/FEDSpreprocessed/TestDefault/TestDefault.json",
        ),
        (
            ["TestOverride", None],
            "s3",
            "s3://big-bucket/small-path/whatever",
            f"s3://big-bucket/small-path/whatever/FEDSpreprocessed/TestOverride/TestOverride.json",
        ),
    ],
)
def test_preprocessed_region_filename_s3(monkeypatch, region: Region,
                                         location: str, new_s3_path: str,
                                         expected_filepath: str):
    # arrange
    if new_s3_path:
        monkeypatch.setattr(settings, "S3_PATH", new_s3_path)

    # act
    actual_filepath = preprocess.preprocessed_region_filename(region, location)

    # assert
    assert actual_filepath == expected_filepath


@pytest.mark.parametrize(
    "region, new_dir_path, expected_filepath",
    [
        (
            ["TestDefault", None],
            None,
            f"{settings.LOCAL_PATH}/FEDSpreprocessed/TestDefault/TestDefault.json",
        ),
        (
            ["TestOverride", None],
            "big-data/whatever",
            f"big-data/whatever/FEDSpreprocessed/TestOverride/TestOverride.json",
        ),
    ],
)
def test_preprocessed_region_filename_local(monkeypatch, region: Region,
                                         new_dir_path: str,
                                         expected_filepath: str):
    # arrange
    if new_dir_path:
        monkeypatch.setattr(settings, "LOCAL_PATH", new_dir_path)

    # act
    actual_filepath = preprocess.preprocessed_region_filename(region, location="local")

    # assert
    assert actual_filepath == expected_filepath


def test_preprocess_region(tmpdir, monkeypatch):
    # arrange
    monkeypatch.setattr(settings, "LOCAL_PATH", str(tmpdir))
    monkeypatch.setattr(FireMain, "maybe_remove_static_sources", lambda region: region)
    test_region = ["Test123", Polygon([(0, 0), (0, 1), (1, 1), (1, 0), (0, 0)])]

    # act
    output_filepath = preprocess.preprocess_region(test_region)

    # assert
    with open(output_filepath, "r") as f:
        assert test_region[1] == from_geojson(f.read())

    assert str(tmpdir) in output_filepath


def test_read_region(tmpdir, monkeypatch):
    # arrange
    monkeypatch.setattr(settings, "S3_PATH", str(tmpdir))
    expected_region = ("Test123", Polygon([(0, 0), (0, 1), (1, 1), (1, 0), (0, 0)]))

    data_dir = tmpdir / settings.PREPROCESSED_DIR
    data_dir.mkdir()

    expected_dir = tmpdir / settings.PREPROCESSED_DIR / expected_region[0]
    expected_dir.mkdir()
    
    with open(expected_dir / f"{expected_region[0]}.json", "w") as f:
        f.write(to_geojson(expected_region[1]))

    # act
    actual_region = preprocess.read_region(expected_region)

    # assert
    assert expected_region == actual_region


def test_preprocess_landcover(tmpdir, mock_rasterio, monkeypatch):
    # arrange
    monkeypatch.setattr(settings, "LOCAL_PATH", str(tmpdir))
    monkeypatch.setattr(settings, "S3_PATH", str(tmpdir))

    data_dir = tmpdir / settings.INPUT_DIR
    data_dir.mkdir()
    test_dir = data_dir / "NLCD"
    test_dir.mkdir()

    file_basename = "nlcd_export_510m_simplified"
    file_path = os.path.join(settings.dirextdata, "NLCD", f"{file_basename}.tif")
    with open(file_path, "w") as f:
        f.write("")

    # act
    preprocess.preprocess_landcover()

    # assert
    actual_first_call_args = mock_rasterio.open.call_args_list[0]
    assert actual_first_call_args == call(file_path)

    actual_second_call_args = mock_rasterio.open.call_args_list[1]
    assert actual_second_call_args == call(
        os.path.join(settings.LOCAL_PATH, settings.PREPROCESSED_DIR, f"{file_basename}_latlon.tif"), "w"
    )

    mock_rasterio.warp.calculate_default_transform.assert_called()

    # TODO: we need a fixture that creates a TIFF
    # mock_rasterio.warp.reproject.assert_called()


@pytest.mark.parametrize(
    "timestep, sat",
    [
        ((2023, 11, 9, "AM"), "NOAA20"),
        ((2023, 11, 9, "AM"), "SNPP"),
    ],
)
def test_preprocess_NRT_file(timestep: TimeStep, sat: str, monkeypatch, test_data_dir):
    monkeypatch.setattr(settings, "READ_LOCATION", "local")
    monkeypatch.setattr(settings, "LOCAL_PATH", test_data_dir)

    if sat == "SNPP":
        df_filtered_paths = preprocess.preprocess_NRT_file(timestep, sat=sat)
    else:
        df_filtered_paths = preprocess.preprocess_NRT_file(timestep, sat=sat)
    assert len(df_filtered_paths) == 2

    # TODO: more assertions on the filtered CSVs


@pytest.mark.parametrize(
    "region, region_shape_to_filter, output_should_already_exist",
    [
        (
            ["TestRegion", Polygon([(0, 0), (0, 2), (2, 2), (2, 0), (0, 0)])],
            Polygon([(0, 0), (0, 1), (1, 1), (1, 0), (0, 0)]),
            False,
        ),  # preprocessed_nrt_snpp_tmpfile fixture should only have one pixel that fits in this region shape
        (
            ["TestRegion", Polygon([(0, 0), (0, 2), (2, 2), (2, 0), (0, 0)])],
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
    # return a preprocessed SNPP file fixture with a couple pixels
    input_df = pd.read_csv(preprocessed_nrt_snpp_tmpfile)

    monkeypatch.setattr(preprocess, "read_preprocessed_input", lambda t, sat=None, location=None: input_df)
    monkeypatch.setattr(preprocess, "read_region", lambda x, location=None: region)
    monkeypatch.setattr(FireIO, "get_reg_shp", lambda x: region_shape_to_filter)

    if output_should_already_exist:
        # mock `FireClustering.do_clustering` so we can assert it is not called
        mock_do_clustering = MagicMock()
        monkeypatch.setattr(FireClustering, "do_clustering", mock_do_clustering)

        monkeypatch.setattr(
            preprocess,
            "preprocessed_filename",
            lambda x: "force_exit.txt",
        )

        monkeypatch.setattr(os.path, "exists", lambda x: True)

        mock_do_clustering.assert_not_called()

    else:
        # we should have separate unit tests that really exercise `FireClustering.do_clustering`
        # and `FireClustering.compute_all_spatial_distances` so just return the input DataFrame for now
        monkeypatch.setattr(FireClustering, "do_clustering", lambda x, y: input_df)
        monkeypatch.setattr(settings, "FIRE_SOURCE", "TESTING123")
        
        # we pass a fake fire source so we can trigger the `else` branch for a single `read_preprocessed` mock
        outfile_df_path = preprocess.preprocess_region_t(
            (2023, 11, 9, "AM"), region
        )

        # assert
        assert outfile_df_path
        assert os.path.exists(outfile_df_path)
        assert len(pd.read_csv(outfile_df_path)) > 0
        os.remove(outfile_df_path)



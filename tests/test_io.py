import os
from unittest.mock import MagicMock

import pytest
import geopandas as gpd
import pandas as pd
from shapely.geometry import Point, Polygon

import fireatlas
from fireatlas import FireIO


@pytest.mark.parametrize(
    "filename, parquet, expect_exception",
    [
        ("valid.geojson", False, None),  # Valid GeoJSON file
        ("invalid.geojson", False, Exception),  # Max retries exceeded
    ],
)
def test_gpd_read_geojson(geojson_file, filename, parquet, expect_exception):

    if expect_exception:
        with pytest.raises(Exception) as exc_info:
            FireIO.gpd_read_file(filename, parquet=parquet)
            assert str(exc_info.value) == ""
    else:
        data = FireIO.gpd_read_file(geojson_file, parquet=parquet)
        assert isinstance(data, gpd.GeoDataFrame)


@pytest.mark.parametrize(
    "filename, parquet, expect_exception",
    [
        ("valid.csv", False, None),  # Valid GeoJSON file
        ("invalid.csv", False, Exception),  # Max retries exceeded
    ],
)
def test_gpd_read_csv(csv_file, filename, parquet, expect_exception):

    if expect_exception:
        with pytest.raises(Exception) as exc_info:
            FireIO.gpd_read_file(filename, parquet=parquet)
            assert str(exc_info.value) == ""
    else:
        data = FireIO.gpd_read_file(csv_file, parquet=parquet)
        assert isinstance(data, gpd.GeoDataFrame)


@pytest.mark.parametrize(
    "filename, parquet, expect_exception",
    [
        ("valid.csv", False, None),  # Valid GeoJSON file
        ("invalid.csv", False, Exception),  # Max retries exceeded
    ],
)
def test_gpd_read_static_source(
    static_source_file_real, filename, parquet, expect_exception
):

    if expect_exception:
        with pytest.raises(Exception) as exc_info:
            FireIO.gpd_read_file(filename, parquet=parquet)
            assert str(exc_info.value) == ""
    else:
        data = FireIO.gpd_read_file(static_source_file_real, parquet=parquet)
        assert isinstance(data, gpd.GeoDataFrame)


@pytest.mark.parametrize(
    "input_value, expected_geometry",
    [
        (Point(0, 0), Point(0, 0)),  # Test with a Shapely Point geometry
        (
            Polygon([(0, 0), (1, 0), (1, 1), (0, 1)]),
            Polygon([(0, 0), (1, 0), (1, 1), (0, 1)]),
        ),  # Test with a Shapely Polygon geometry
        ("United Flakes", None),  # Test with an invalid country name
        (
            [0, 0, 1, 1],
            Polygon([(0, 0), (0, 1), (1, 1), (1, 0), (0, 0)]),
        ),  # Test with a list representing an extent
        (123, None),  # Test with an invalid input type
    ],
)
def test_get_reg_shp(input_value, expected_geometry, monkeypatch):
    monkeypatch.setattr(FireIO, "get_Cty_shp", lambda reg: None)
    result = FireIO.get_reg_shp(input_value)
    assert result == expected_geometry


@pytest.mark.parametrize(
    "region_shape_to_filter, all_pixels_inside",
    [
        (Polygon([(0, 0), (0, 1), (1, 1), (1, 0), (0, 0)]), False),
        (Polygon([(0, 0), (0, 3), (3, 3), (3, 0), (0, 0)]), True),
    ],
)
def test_afp_regfilter(
    preprocessed_nrt_snpp_tmpfile, region_shape_to_filter, all_pixels_inside
):
    # arrange
    input_df = pd.read_csv(preprocessed_nrt_snpp_tmpfile)

    # act
    df_expected_count = len(input_df)
    df_filtered = FireIO.AFP_regfilter(input_df, region_shape_to_filter)

    # assert
    if all_pixels_inside:
        assert df_expected_count == len(df_filtered)
    else:
        assert df_expected_count > len(df_filtered)


@pytest.mark.parametrize(
    "local_settings_dir, s3_settings_dir, filename, fs_mock",
    [
        ("/data", "s3://magic-mike", "firelines.fgb", MagicMock()),
        ("/data", "s3://magic-mike", "allfires.parq", MagicMock()),
    ],
)
def test_copy_from_local_to_s3(
    tmp_settings_context_manager,
    local_settings_dir,
    s3_settings_dir,
    filename,
    fs_mock
):
    # arrange
    expected_local_filepath = os.path.join(local_settings_dir, filename)
    expected_s3_filepath = os.path.join(s3_settings_dir, filename)
    fs_mock.put_file = MagicMock()

    # act
    with tmp_settings_context_manager(
        fireatlas.settings,
        LOCAL_PATH=local_settings_dir,
        S3_PATH=s3_settings_dir
    ):
        FireIO.copy_from_local_to_s3(expected_local_filepath, fs_mock)

        # assert
        fs_mock.put_file.assert_called_with(expected_local_filepath, expected_s3_filepath)

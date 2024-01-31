import pandas as pd
import numpy as np
import pytest
import FireClustering


@pytest.mark.parametrize(
    "data, max_thresh_km, expected",
    [
        (
            pd.DataFrame({"x": [0, 1], "y": [0, 1]}),
            1,
            np.array([np.array([1]), np.array([0])], dtype=object),
            # testing pixels both within distances of each other
        ),
        (
            pd.DataFrame({"x": [0, 10000], "y": [0, 11000]}),
            5,
            np.array([np.array([]), np.array([])], dtype=object),
            # testing pixels that are not within distances of each other
        ),
    ],
)
def test_compute_all_spatial_distances(data, max_thresh_km, expected):
    result = FireClustering.compute_all_spatial_distances(data, max_thresh_km)
    assert np.array_equal(result, expected)


@pytest.mark.parametrize(
    "data, max_thresh_km",
    [
        (
            pd.DataFrame(
                {
                    "Lat": [0.9, 2],
                    "Lon": [0.5, 0.7],
                    "FRP": [0.5, 0.7],
                    "Sat": ["SNPP", "SNPP"],
                    "DT": [0.47, 0.44],
                    "DS": [0.46, 0.4],
                    "YYYYMMDD_HHMM": ["2023-08-28 00:03:00", "2023-08-28 00:03:00"],
                    "ampm": ["AM", "AM"],
                }
            ),
            1,
        ),
        (
            pd.DataFrame(
                {
                    "Lat": [0.9, 1000],
                    "Lon": [0.5, 2000.7],
                    "FRP": [0.5, 0.7],
                    "Sat": ["SNPP", "SNPP"],
                    "DT": [0.47, 0.44],
                    "DS": [0.46, 0.4],
                    "YYYYMMDD_HHMM": ["2023-08-28 00:04:00", "2023-08-28 00:04:00"],
                    "ampm": ["AM", "AM"],
                }
            ),
            1,
        ),
    ],
)
def test_do_clustering(data, max_thresh_km, expected):
    clustered_df = FireClustering.do_clustering(input_df, max_thresh_km)
    import pdb; pdb.set_trace()
    assert np.array_equal(clustered_df, expected)

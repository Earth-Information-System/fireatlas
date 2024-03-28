import pandas as pd
import numpy as np
import pytest
from fireatlas import FireClustering


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
    "data, max_thresh_km, should_cluster",
    [
        (
            # when < three pixels always cluster
            # NOTE: this seems wrong!
            pd.DataFrame(
                {
                    "Lon": [1, 20000],
                    "Lat": [1, 20000],
                    "x": [1, 20000],
                    "y": [1, 20000],
                    "FRP": [0.5, 0.7],
                    "Sat": ["SNPP", "SNPP"],
                    "DT": [0.47, 0.44],
                    "DS": [0.46, 0.4],
                    "datetime": ["2023-08-28 00:03:00", "2023-08-28 00:03:00"],
                    "ampm": ["AM", "AM"],
                }
            ),
            1,
            True,
        ),
        (
            # cluster things close together
            pd.DataFrame(
                {
                    "Lon": [1, 2, 300000],
                    "Lat": [1, 2, 300000],
                    "x": [1, 2, 300000],
                    "y": [1, 2, 300000],
                    "FRP": [0.5, 0.7, 0.9],
                    "Sat": ["SNPP", "SNPP", "SNPP"],
                    "DT": [0.47, 0.44, 0.5],
                    "DS": [0.46, 0.4, 0.5],
                    "datetime": [
                        "2023-08-28 00:03:00",
                        "2023-08-28 00:03:00",
                        "2023-08-28 00:03:00",
                    ],
                    "ampm": ["AM", "AM", "AM"],
                }
            ),
            1,
            True,
        ),
        (
            # do not cluster anything so far aparam
            pd.DataFrame(
                {
                    "Lon": [1, 200000, 300000],
                    "Lat": [1, 200000, 30000],
                    "x": [1, 200000, 300000],
                    "y": [1, 200000, 300000],
                    "FRP": [0.5, 0.7, 0.9],
                    "Sat": ["SNPP", "SNPP", "SNPP"],
                    "DT": [0.47, 0.44, 0.5],
                    "DS": [0.46, 0.4, 0.5],
                    "datetime": [
                        "2023-08-28 00:03:00",
                        "2023-08-28 00:03:00",
                        "2023-08-28 00:03:00",
                    ],
                    "ampm": ["AM", "AM", "AM"],
                }
            ),
            1,
            False,
        ),
    ],
)
def test_do_clustering(data, max_thresh_km, should_cluster):
    clustered_df = FireClustering.do_clustering(data, max_thresh_km)
    if should_cluster:
        if isinstance(clustered_df, list):
            # FireCluster.do_clustering will send back indices of all pixels as a list when number of pixels < 3
            assert len(clustered_df) == 2
        else:
            cluster_series = clustered_df.groupby("initial_cid").size()
            assert cluster_series.count() == 2
    else:
        # each pixel is a cluster
        cluster_series = clustered_df.groupby("initial_cid").size()
        assert cluster_series.count() == 3

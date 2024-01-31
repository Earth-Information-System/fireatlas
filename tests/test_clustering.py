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
        ),
        (
            pd.DataFrame({"x": [0, 10000], "y": [0, 11000]}),
            5,
            np.array([np.array([]), np.array([])], dtype=object),
        ),
    ],
)
def test_compute_all_spatial_distances(data, max_thresh_km, expected):
    result = FireClustering.compute_all_spatial_distances(data, max_thresh_km)
    assert np.array_equal(result, expected)

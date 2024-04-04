import pytest

from fireatlas import postprocess


@pytest.mark.parametrize("location", ["s3", "local"])
def test_all_dir(location):
    tst = (2023, 11, 9, "AM")
    region = ["TESTING123", None]
    all_dir = postprocess.all_dir(tst, region, location=location)
    assert "/FEDSoutput-v3/TESTING123/2023" in all_dir

    if location == "s3":
        assert all_dir.startswith("s3://")
    else:
        assert "/data/" in all_dir


@pytest.mark.parametrize("location", ["s3", "local"])
def test_snapshot_folder(location):
    tst = (2023, 11, 9, "AM")
    ted = (2023, 11, 9, "PM")
    region = ["TESTING123", None]
    snapshot_folder = postprocess.snapshot_folder(tst, ted, region, location=location)
    assert "/FEDSoutput-v3/TESTING123/2023/Snapshot" in snapshot_folder


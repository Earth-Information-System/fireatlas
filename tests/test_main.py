import FireConsts
from FireMain import maybe_remove_static_sources
from FireTypes import Region
from shapely.geometry import Polygon
import pytest


@pytest.mark.parametrize(
    "original_region, should_remove_static_sources_bool",
    [
        (["Test1", [0, 0, 1, 1]], True),  # remove static sources
        (
            ["Test2", [0, 0, 1, 1]],
            False,
        ),  # do not remove static sources and compare inputs
    ],
)
def test_maybe_remove_static_sources(
    static_source_dir_fake,
    original_region: Region,
    should_remove_static_sources_bool: bool,
):
    # arrange
    if not should_remove_static_sources_bool:
        FireConsts.remove_static_sources_bool = False
    else:
        FireConsts.remove_static_sources_bool = True

    # act
    result_region = maybe_remove_static_sources(original_region, static_source_dir_fake)

    # assert
    if not should_remove_static_sources_bool:
        assert result_region == original_region
    else:
        assert result_region != original_region
        name, geom = result_region
        assert isinstance(geom, Polygon)
        # check if there are any interior rings (holes)
        assert len(geom.interiors) > 0

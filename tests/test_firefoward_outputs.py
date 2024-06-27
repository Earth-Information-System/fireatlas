import pytest
import fireatlas


def test_fireforward_v2_vs_v3(
    tmp_settings_context_manager,
    test_data_dir,
):
    with tmp_settings_context_manager(
        fireatlas.settings,
        LOCAL_PATH=test_data_dir
    ):
        # now all local reads come from fireatlas_nrt/tests/data instead of fireatlas_nrt/data
        pass

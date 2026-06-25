"""Tests for the Settings field validators in fireatlas.FireConsts."""

import pytest
from pydantic import ValidationError

from fireatlas.FireConsts import Settings


def test_local_path_strips_trailing_slash():
    s = Settings(LOCAL_PATH="/some/data/")
    assert s.LOCAL_PATH == "/some/data"


def test_local_path_without_trailing_slash_unchanged():
    s = Settings(LOCAL_PATH="/some/data")
    assert s.LOCAL_PATH == "/some/data"


def test_s3_path_strips_trailing_slash():
    s = Settings(S3_PATH="s3://bucket/prefix/")
    assert s.S3_PATH == "s3://bucket/prefix"


def test_s3_path_without_s3_prefix_raises():
    with pytest.raises(ValidationError):
        Settings(S3_PATH="bucket/prefix")


def test_s3_path_valid_unchanged():
    s = Settings(S3_PATH="s3://bucket/prefix")
    assert s.S3_PATH == "s3://bucket/prefix"

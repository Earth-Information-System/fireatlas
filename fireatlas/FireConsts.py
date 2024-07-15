""" FireConsts
This is the module containing all constants used in this project as well as the
running controls
"""

from typing import Literal
import os

import fsspec
from pydantic_settings import BaseSettings
from pydantic import Field, validator

from fireatlas.FireTypes import Location

root_dir = os.path.abspath(os.path.dirname(os.path.dirname(__file__)))


class Settings(BaseSettings):
    # read in all env vars prefixed with `FEDS_` they can be in a .env file
    model_config = {
        "env_file": ".env",
        "extra": "ignore",
        "env_prefix": "FEDS_",
    }
    # ------------------------------------------------------------------------------
    # where data is stored
    # ------------------------------------------------------------------------------
    LOCAL_PATH: str = Field(
        os.path.join(root_dir, "data"),
        description="absolute path to where local data are stored",
    )
    S3_PATH: str = Field(
        "s3://maap-ops-workspace/shared/gsfc_landslides",
        description="s3 path where remote data are stored",
    )

    INPUT_DIR: str = Field(
        "FEDSinput", description="directory where input data is stored"
    )
    PREPROCESSED_DIR: str = Field(
        "FEDSpreprocessed", description="directory where preprocessed data is stored"
    )
    OUTPUT_DIR: str = Field(
        "FEDSoutput-v3", description="directory where output data is stored"
    )

    READ_LOCATION: Location = Field(
        "s3",
        description="Final storage place for written files. This is where everything reads from",
    )

    LOG_FILENAME: str = Field("running.log", description="Where to write logs to.")

    # ------------------------------------------------------------------------------
    # spatiotemporal constraints of fire objects
    # ------------------------------------------------------------------------------
    # spatial parameters used for fire pixel clustering
    EARTH_RADIUS_KM: float = Field(6371.0, description="earth radius, km")

    LARGEFIRE_FAREA: int = Field(
        4, description="fire area threshold for determining large fires"
    )
    EPSG_CODE: Literal[9311, 32610, 3571] = Field(
        9311,
        description="epsg projection code ( 3571: North Pole LAEA; 32610: WGS 84 / UTM zone 10N; 9311: US National Atlas Equal Area)",
    )

    # temporal parameters for fire object definition
    maxoffdays: int = Field(
        5,
        description="fire becomes inactive after this number of consecutive days without active fire detection",
    )
    limoffdays: int = Field(
        20,
        description="fire keeps sleeper status even at inactive but with inactive dates smaller than this value",
    )
    CONNECTIVITY_CLUSTER_KM: float = Field(
        0.7, description="the connectivity spatial threshold for initial clustering, km"
    )
    CONNECTIVITY_SLEEPER_KM: float = Field(
        1, description="the connectivity spatial threshold (to previous fire line), km"
    )

    # ------------------------------------------------------------------------------
    # shape parameters
    # ------------------------------------------------------------------------------
    valpha: int = Field(1000, description="alpha parameter, m")

    # VIIRS pixel size
    VIIRSbuf: float = Field(187.5, description="fire perimeter buffer, m")
    fpbuffer: int = Field(
        200, description="buffer use to determine fire line pixels, m"
    )
    flbuffer: int = Field(
        500,
        description="buffer for fire line pixels (radius) to intersect fire perimeter, m",
    )

    extbuffer: int = Field(
        1000, description="buffer to define interior/exterior region, m"
    )
    area_VI: float = Field(0.141, description="area of each 375m VIIRS pixel, km2")

    # MODIS pixel size
    MCD64buf: float = Field(231.7, description="MODIS fire perimeter buffer, m")

    # fire source data
    FIRE_SOURCE: Literal["SNPP", "NOAA20", "VIIRS", "BAMOD"] = Field(
        "NOAA20", description="fire source data"
    )
    FIRE_NRT: bool = Field(True, description="whether to use NRT data")
    FIRE_SENSOR: Literal["viirs", "mcd64"] = Field("viirs", description="fire sensor")

    # ------------------------------------------------------------------------------
    # static fire parameters
    # ------------------------------------------------------------------------------
    remove_static_sources: bool = Field(
        True, description="remove areas with known flaring/gas sources from region"
    )
    remove_static_sources_sourcefile: str = Field(
        "VIIRS_Global_flaring_d.7_slope_0.029353_2017_web_v1.csv",
        description="File where static sources are stored",
    )
    remove_static_sources_buffer: float = Field(
        0.01, description="Buffer around static source points. Units defined by epsg"
    )

    remove_static_small_fires: bool = Field(
        False, description="remove small fires with high pixel density"
    )

    # ------------------------------------------------------------------------------
    # other options
    # ------------------------------------------------------------------------------
    # fire tracking options
    expand_only: bool = Field(
        False,
        description="if set to true, only expand existing fires (no new fire objects created)",
    )
    number_of_multi_proc_workers: int = Field(
        3, description="number of dask process workers to use"
    )
    export_to_veda: bool = Field(
        False, description="whether to export data from MAAP to VEDA s3"
    )
    N_DASK_WORKERS: int = Field(
        6, description="How many dask workers to use for Run."
    )

    # ------------------------------------------------------------------------------
    # fire type related parameters
    # ------------------------------------------------------------------------------
    FTYP_OPT: Literal["preset", "CA", "global"] = Field(
        "CA", description="fire type option"
    )
    CONT_OPT: Literal["preset", "CA", "global"] = Field(
        "CA", description="continuity threshold option"
    )

    @validator("LOCAL_PATH")
    def local_path_must_not_end_with_slash(cls, v: str) -> str:
        if v.endswith("/"):
            v = v[:-1]
        return v

    @validator("S3_PATH")
    def s3_path_must_start_with_s3(cls, v: str) -> str:
        if not v.startswith("s3://"):
            raise ValueError("S3_PATH must start with s3://")
        if v.endswith("/"):
            v = v[:-1]
        return v

    @property
    def fs(self):
        return fsspec.filesystem(self.READ_LOCATION, use_listings_cache=False)

    @property
    def dirextdata(self):
        return f"{self.get_path()}/{self.INPUT_DIR}/"

    @property
    def diroutdata(self):
        return f"{self.get_path()}/{self.OUTPUT_DIR}/"

    def get_path(self, location: Location = None):
        """Path to data - dependent on specified location or READ_LOCATION"""
        if location is None:
            location = self.READ_LOCATION

        if location == "local":
            return self.LOCAL_PATH
        else:
            return self.S3_PATH


FTYP_preset = 2
FTYP = {
    "preset": {2: "Forest"},  # use 'Forest' for all fires
    "CA": {
        0: "Other",
        1: "Urban",
        2: "Forest wild",
        3: "Forest manage",
        4: "Shrub wild",
        5: "Shrub manage",
        6: "Agriculture",
    },  # use algorithm in CAFEDS
    "global": {
        0: "Other",
        1: "Temp Forest",
        2: "Trop Forest",
        3: "Bore Forest",
        4: "Savana",
        5: "Agriculture",
        6: "Deforestation",
    },  #  use algorithm proposed for global study
}

CONT = {
    "preset": {2: 1},  # use 1 for all fires
    "CA": {
        0: 1,
        1: 1,
        2: 2.5,
        3: 5,
        4: 5,
        5: 5,
        6: 1,
    },  # fire type dependent CONNECTIVITY_THRESHOLD_KM
}

FTYPCLR_CA = {
    0: "grey",
    1: "rosybrown",
    2: "darkolivegreen",
    3: "olive",
    4: "saddlebrown",
    5: "sandybrown",
    6: "darkviolet",
}  # colors used for each fire type

""" FireConsts
This is the module containing all constants used in this project as well as the
running controls
"""

from typing import Literal, Optional, Tuple
import os
import warnings

import fsspec
from pydantic_settings import (
    BaseSettings, 
    SettingsConfigDict,
    PydanticBaseSettingsSource, 
    YamlConfigSettingsSource,
)
from pydantic import Field, validator, field_validator


from fireatlas.FireTypes import Location

root_dir = os.path.abspath(os.path.dirname(os.path.dirname(__file__)))

DOTENV_ABS_PATH = os.path.join(os.path.dirname(__file__), ".env")
YAML_FILENAME = "run_config.yaml"
YAML_ABS_PATH = os.path.join(os.path.dirname(__file__), YAML_FILENAME)

class Settings(BaseSettings):

    if os.path.exists(YAML_ABS_PATH):
        model_config = SettingsConfigDict(
            yaml_file=YAML_ABS_PATH, 
            yaml_config_section="settings",
            env_file=DOTENV_ABS_PATH, 
            extra="ignore",
            env_prefix="FEDS_" # note: expects env vars as FEDS_env_var_name
        )
    else:
        model_config = SettingsConfigDict(
            env_file=DOTENV_ABS_PATH, 
            extra="ignore",
            env_prefix="FEDS_"
        )

    # Settings resolution order (in ascending order of priority): 
    # 1. Starts with default field values provided in FireConsts.py
    # 2. Overrides with settings from run_config.yaml if available
    # 3. Overrides with environment variables from .env file if available
    # 4. Overrides with environment variables from the system env if set
    # 5. Overrides with settings passed to the Settings class initializer if passed
    # Result: init > system env vars > .env file > yaml file > default values
    @classmethod
    def settings_customise_sources(
        cls, 
        settings_cls: type[BaseSettings],
        init_settings: PydanticBaseSettingsSource, 
        env_settings: PydanticBaseSettingsSource, 
        dotenv_settings: PydanticBaseSettingsSource, 
        file_secret_settings: PydanticBaseSettingsSource,
    ) -> tuple[PydanticBaseSettingsSource, ...]:
        return (
            init_settings, 
            env_settings, 
            dotenv_settings,
            YamlConfigSettingsSource(settings_cls), 
            )



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

    REGIONS_DIR: str = Field(
        "run_definitions",
        description="directory where region definitions are stored." 
    )

    READ_LOCATION: Location = Field(
        "s3",
        description="Final storage place for written files. This is where everything reads from",
    )

    LOG_FILEPATH: str = Field(
        os.path.join(root_dir, "running.log"),
        description="Absolute path to the log file."
    )

    ENV_META_FILEPATH: str = Field(
        os.path.join(root_dir, "env_metadata.txt"),
        description="Absolute path to the environment metadata file."
    )
    
    # ------------------------------------------------------------------------------
    # spatiotemporal constraints of fire objects
    # ------------------------------------------------------------------------------
    # spatial parameters used for fire pixel clustering
    EARTH_RADIUS_KM: float = Field(6371.0, description="earth radius, km")

    LARGEFIRE_FAREA: int = Field(
        4, description="fire area threshold for determining large fires"
    )

    EPSG_CODE: int = Field(
        9311,
        description="epsg projection code ( 3571: North Pole LAEA; 32610: WGS 84 / UTM zone 10N; 9311: US National Atlas Equal Area)",
    )

    @field_validator("EPSG_CODE")
    @classmethod
    def check_epsg(cls, epsg: int):
        allowed = (3571, 32610, 9311, 6933)
        if epsg not in allowed:
            warnings.warn(
                f"EPSG projection code {epsg} not recognized as one of: {allowed}. (A new code can be registered in FireConsts.py if needed.) The code should only be run with a projected coordinate system."
            )
        return epsg

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
    # OPTIONAL: run parameters 
    # Can be passed from run_config.yaml if using FireRunArchiveCoordinator.py
    # or passed from the command line for all other scripts
    # ------------------------------------------------------------------------------

    TST: Optional[Tuple[int, int, int, Literal["AM", "PM"]]] = Field(
        default=None, 
        description="start time as [year, month, day, 'AM'/'PM']"
    )

    TED: Optional[Tuple[int, int, int, Literal["AM", "PM"]]] = Field(
        default=None, 
        description="end time as [year, month, day, 'AM'/'PM']"
    )

    @field_validator("TST", "TED", mode="after")
    @classmethod
    def _tuple_to_list(cls, v):
        # v is already validated as a tuple of (int, int, int, str)
        if v is None:
            return None
        return list(v)

    RUN_NAME: Optional[str] = Field(
        default=None,
        description="Run name, e.g. 'ArchiveCONUS' or 'ArchiveCONUS_test'."
    )

    REGION_SHAPEFILE: Optional[str] = Field(
        default=None, 
        description="Name of the file that holds a shapefile that defiens this region. " 
        "Assumes that this file is in the FEDSinput/run_definitions/RUN_NAME/ directory."
    )

    REGION_BBOX: Optional[list[float]] = Field(
        default=None,
        description="Bounding box of the region, e.g. [-126,24,-61,49]."
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

    # ------------------------------------------------------------------------------
    # fire data source parameters
    # ------------------------------------------------------------------------------
    
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

    # compute settings

    N_DASK_WORKERS: int = Field(6, description="How many dask workers to use for Run.")

    ARCHIVE_RUN_JOB_SIZE: int = Field(
        10, 
        description="How many days to run in each archive job chunk."
    )

    # NIFC matching options
    DO_NIFC_MATCHING: bool = Field(
        False, 
        description="If True, reads from the NIFC incident database for current year and adds cols with info for matching fires to the combinedLargefire perimeter fgb output."
    )

    NIFC_MATCHING_ACTIVE_ONLY: bool = Field(
        False, 
        description="If True, uses 'WFIGS Current' NIFC database. Else, uses 'WFIGS {current year} to date'."
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

    # ------------------------------------------------------------------------------


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

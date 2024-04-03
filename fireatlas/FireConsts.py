""" FireConsts
This is the module containing all constants used in this project as well as the
running controls
"""

from typing import Literal
import os
from fireatlas.FireEnums import FireSource, EPSG
from functools import partial
from pydantic_settings import BaseSettings
from pydantic import Field


def get_env_var_as_type(name, cast_to_type=int, default=None):
    try:
        if cast_to_type == bool:
            return os.environ[name].lower() in ["true", "1", "t", "y", "yes"]
        return cast_to_type(os.environ[name])
    except KeyError:
        if default is not None:
            return default
        raise ValueError(f"Environment variable '{name}' not found")
    except ValueError:
        raise ValueError(
            f"Environment variable '{name}' could not be converted to an '{cast_to_type}'"
        )


# ------------------------------------------------------------------------------
# project directories
# ------------------------------------------------------------------------------

projnm = "FEDStest"  # project name
dirhome = os.environ.get("HOME")  # get system home directory

dirproject = os.path.abspath(os.path.dirname(os.path.dirname(__file__)))
dirdata_local_path = os.path.join(dirproject, "data")

dirdata = "./"  # project directory -- only used For logging location
dirdata_s3_bucket = "maap-ops-workspace"
dirdata_s3_path = "shared/gsfc_landslides"


def get_dirdata(
    dirname: Literal["FEDSinput", "FEDSpreprocessed", "FEDSoutput-v3"],
    location: Literal["s3", "local"],
):
    if location == "local":
        return f"{dirdata_local_path}/{dirname}/"
    else:
        return f"s3://{dirdata_s3_bucket}/{dirdata_s3_path}/{dirname}/"


get_dirextdata = partial(get_dirdata, dirname="FEDSinput")
get_dirprpdata = partial(get_dirdata, dirname="FEDSpreprocessed")
get_diroutdata = partial(get_dirdata, dirname="FEDSoutput-v3")

# final storage place for written files. This is where everything reads from
READ_LOCATION = "s3"
dirextdata = get_dirextdata(location=READ_LOCATION)
dirprpdata = get_dirprpdata(location=READ_LOCATION)
diroutdata = get_diroutdata(location=READ_LOCATION)

# lakedir = 'D:/fire_atlas/Data/GlobalSurfaceWater/vector/'

class Settings(BaseSettings):
    # ------------------------------------------------------------------------------
    # spatiotemporal constraints of fire objects
    # ------------------------------------------------------------------------------
    # spatial parameters used for fire pixel clustering
    EARTH_RADIUS_KM: float = Field(6371.0, description="earth radius, km")

    # temporal parameters for fire object definition
    maxoffdays: int = Field(5, description=(
        "fire becomes inactive after this number of "
        "consecutive days without active fire detection"
    ))
    limoffdays: int = Field(20, description=(
        "fire keeps sleeper status even at inactive but with "
        "inactive dates smaller than this value"
    ))
    CONNECTIVITY_CLUSTER_KM: float = Field(0.7, description=(
        "the connectivity spatial threshold for initial clustering, km"
    ))
    CONNECTIVITY_SLEEPER_KM: float = Field(1, description=(
        "the connectivity spatial threshold (to previous fire line), km"
    ))

    # spatial parameters for fire object categorization
    LARGEFIRE_FAREA: int = Field(4, description="fire area threshold for determining large fires")

    # fire tracking options
    expand_only: bool = Field(False, description=(
        "if set to true, only expand existing fires (no new fire objects created)"
    ))

    # ------------------------------------------------------------------------------
    # shape parameters
    # ------------------------------------------------------------------------------

    valpha: int = Field(1000, description="alpha parameter, m")

    # VIIRS pixel size
    VIIRSbuf: float = Field(187.5, description="fire perimeter buffer, m")
    fpbuffer: int = Field(200, description="buffer use to determine fire line pixels, m")
    flbuffer: int = Field(500, description="buffer for fire line pixels (radius) to intersect fire perimeter, m")

    extbuffer: int = Field(1000, description="buffer to define interior/exterior region, m")
    area_VI: float = Field(0.141, description="area of each 375m VIIRS pixel, km2")

    # MODIS pixel size
    MCD64buf: float = Field(231.7, description="MODIS fire perimeter buffer, m")
    
    # fire source data
    FIRE_SOURCE: Literal["SNPP", "NOAA20", "VIIRS", "BAMOD"] = Field("VIIRS", description="fire source data")
    FIRE_NRT: bool = Field(True, description="whether to use NRT data")
    FIRE_SENSOR: Literal["viirs", "mcd64"] = Field("viirs", description="fire sensor")

# ------------------------------------------------------------------------------
# fire type related parameters
# ------------------------------------------------------------------------------

# fire type options
FTYP_opt = get_env_var_as_type("FTYP_OPT", cast_to_type=int, default=1)
# 0: preset ftype for all fires;
# 1: use CA type classifications
# 2: proposed global fire types
CONT_opt = get_env_var_as_type("CONT_OPT", cast_to_type=int, default=1)
# 0: preset continuity threshold for all fires;
# 1: use CA type classifications dependent values
# 2: use global fire types and size dependent values

# For FTYP_opt = 0, use 'forest' for all fires
FTYP_preset = [2, "Forest"]  # the default ftype
CONT_preset = 1  # km, default continuity threshold

# For FTYP_opt = 1, use algorithm in CAFEDS
FTYP_CA = {
    0: "Other",
    1: "Urban",
    2: "Forest wild",
    3: "Forest manage",
    4: "Shrub wild",
    5: "Shrub manage",
    6: "Agriculture",
}  # fire type names
FTYPCLR_CA = {
    0: "grey",
    1: "rosybrown",
    2: "darkolivegreen",
    3: "olive",
    4: "saddlebrown",
    5: "sandybrown",
    6: "darkviolet",
}  # colors used for each fire type
CONT_CA = {
    0: 1,
    1: 1,
    2: 2.5,
    3: 5,
    4: 5,
    5: 5,
    6: 1,
}  # preset fire type dependent CONNECTIVITY_THRESHOLD_KM

# For FTYP_opt = 2, use algorithm proposed for global study
FTYP_Glb = {
    0: "Other",
    1: "Temp Forest",
    2: "Trop Forest",
    3: "Bore Forest",
    4: "Savana",
    5: "Agriculture",
    6: "Deforestation",
}  # fire type names

# ------------------------------------------------------------------------------
# other options
# ------------------------------------------------------------------------------
epsg = get_env_var_as_type(
    "EPSG_CODE", cast_to_type=int, default=EPSG.CONUS_EQ_AREA.value
)
# epsg projection code ( 3571: North Pole LAEA; 32610: WGS 84 / UTM zone 10N; 9311: US National Atlas Equal Area)

remove_static_sources_bool = (
    True  # remove areas with known flaring/gas sources from region
)
remove_static_sources_sourcefile = (
    "VIIRS_Global_flaring_d.7_slope_0.029353_2017_web_v1.csv"
)
remove_static_sources_buffer = (
    0.01  # Buffer around static source points. Units defined by epsg.
)

opt_rmstatfire = False  # do the removal of small fires with high pixel density

number_of_multi_proc_workers = 3

export_to_veda = False

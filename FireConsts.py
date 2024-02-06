""" FireConsts
This is the module containing all constants used in this project as well as the
running controls
"""

import os
import FireEnums


def get_env_var_as_type(name, cast_to_type=int, default=None):
    try:
        if cast_to_type == bool:
            return os.environ[name].lower() in ['true', '1', 't', 'y', 'yes']
        return cast_to_type(os.environ[name])
    except KeyError:
        if default is not None:
            return default
        raise ValueError(f"Environment variable '{name}' not found")
    except ValueError:
        raise ValueError(f"Environment variable '{name}' could not be converted to an '{cast_to_type}'")

# ------------------------------------------------------------------------------
# project directories
# ------------------------------------------------------------------------------

projnm = "FEDStest"  # project name
dirhome = os.environ.get("HOME")  # get system home directory

# run at iMac
dirdata = './'  # project directory -- only used For logging location
# run at MacBook
# dirdata = os.path.join(dirhome,'GoogleDrive','My','My.Research','UCI','ProjectData','CAFEDS','v2.0')

dirextdata = '/projects/shared-buckets/gsfc_landslides/FEDSinput/'  # exterior input data directory
dirprpdata = 'tests/data/processed/'
dirtmpdata = "tests/data/tmp/"
diroutdata = "tests/data/out/"

# lakedir = 'D:/fire_atlas/Data/GlobalSurfaceWater/vector/'

# ------------------------------------------------------------------------------
# spatiotemporal constraints of fire objects
# ------------------------------------------------------------------------------

# spatial parameters used for fire pixel clustering
EARTH_RADIUS_KM = 6371.0  # earth radius, km

# temporal parameters for fire object definition
maxoffdays = 5  # fire becomes inactive after this number of consecutive days without active fire detection
limoffdays = 20  # fire keeps sleeper status even at inactive but with inactive dates smaller than this value
CONNECTIVITY_CLUSTER_KM = 0.7  # the connectivity spatial threshold for initial clustering, km
CONNECTIVITY_SLEEPER_KM = 1  # the connectivity spatial threshold (to previous fire line), km

# fire tracking options
expand_only = (
    False  # if set to true, only expand existing fires (no new fire objects created)
)

# ------------------------------------------------------------------------------
# shape parameters
# ------------------------------------------------------------------------------

# alpha shape
valpha = 1000  # alpha parameter, in m (default)
lalpha = 100
stralpha = "1km"  # string of alpha parameter

# VIIRS pixel size
VIIRSbuf = 187.5  # fire perimeter buffer (deg), corresponding to 375m/2 at lat=30
fpbuffer = 200  # buffer use to determine fire line pixels (deg), ~200m
flbuffer = (
    500  # buffer for fire line pixels (radius) to intersect fire perimeter (deg), ~500m
)
extbuffer = 1000  # buffer to define interior/exterior region, 1000 m
area_VI = 0.141  # km2, area of each 375m VIIRS pixel

# MODIS pixel size
MCD64buf = 231.7  # MODIS fire perimeter buffer (deg), corresponding to 463.31271653 m/2

# fire source data
firesrc = get_env_var_as_type('FIRE_SOURCE', cast_to_type=str, default=FireEnums.FireSource.VIIRS.value)  # source - ['SNPP', 'NOAA20', 'VIIRS', 'BAMOD']:
firenrt = get_env_var_as_type('FIRE_NRT', cast_to_type=bool, default=True) # NRT - True, False
firessr = "viirs"  # sensor - 'mcd64'

# ------------------------------------------------------------------------------
# fire type related parameters
# ------------------------------------------------------------------------------

# fire type options
FTYP_opt = get_env_var_as_type('FTYP_OPT', cast_to_type=int, default=1)
# 0: preset ftype for all fires;
# 1: use CA type classifications
# 2: proposed global fire types
CONT_opt = get_env_var_as_type('CONT_OPT', cast_to_type=int, default=1)
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
epsg = get_env_var_as_type('EPSG_CODE', cast_to_type=int, default=FireEnums.EPSG.CONUS_EQ_AREA.value)
# epsg projection code ( 3571: North Pole LAEA; 32610: WGS 84 / UTM zone 10N; 9311: US National Atlas Equal Area)

remove_static_sources_bool = True  # remove areas with known flaring/gas sources from region
remove_static_sources_sourcefile = "VIIRS_Global_flaring_d.7_slope_0.029353_2017_web_v1.csv"
remove_static_sources_buffer = 0.01 # Buffer around static source points. Units defined by epsg. 

opt_rmstatfire = True  # do the removal of small fires with high pixel density

number_of_multi_proc_workers = 3

export_to_veda = False

""" FireConsts
This is the module containing all constants used in this project as well as the
running controls
"""

import os

# ------------------------------------------------------------------------------
# project directories
# ------------------------------------------------------------------------------

projnm = 'FEDStest'
dirhome = os.environ.get('HOME')  # home directory

# run at iMac
dirdata = os.path.join(dirhome, 'Jacaranda', 'ProjectData', projnm)  # project directory
# run at MacBook
# dirdata = os.path.join(dirhome,'GoogleDrive','My','My.Research','UCI','ProjectData','CAFEDS','v2.0')

dirextdata = os.path.join(dirdata, 'ExtData')  # exterior input data directory
dirtmpdata = os.path.join(dirdata, 'Temp')     # temporary data directory
diroutdata = os.path.join(dirdata, 'Output')   # output data directory

# lakedir = 'D:/fire_atlas/Data/GlobalSurfaceWater/vector/'

# ------------------------------------------------------------------------------
# spatiotemporal constraints of fire objects
# ------------------------------------------------------------------------------

# spatial parameters used for fire pixel clustering
EARTH_RADIUS_KM = 6371.0  # earth radius, km
SPATIAL_THRESHOLD_KM = 1  # threshold of spatial distance (between pixels, km)
                          # used to do fire initial clustering
MAX_THRESH_CENT_KM = 50   # threshold of spatial distance (between clusters centroid, km) used to filter nearby clusters

# temporal parameters for fire object definition
maxoffdays = 5            # fire becomes inactive after this number of consecutive days without active fire detection
limoffdays = 20           # this is the additional threshold that involves checking for active fireline

# sleep fire related parameters
sleeperthresh = 1   # km, distance a sleeper is allowes to spread since its last active fire detection

# fire tracking options
expand_only = False  # if set to true, only expand existing fires (no new fire objects created)
# ------------------------------------------------------------------------------
# shape parameters
# ------------------------------------------------------------------------------

# alpha shape
valpha = 1000      # alpha parameter, in m (default)
lalpha = 100
stralpha = '1km'   # string of alpha parameter

# VIIRS pixel size
VIIRSbuf = 187.5   # fire perimeter buffer (deg), corresponding to 375m/2 at lat=30
fpbuffer = 200     # buffer use to determine fire line pixels (deg), ~200m
flbuffer = 500     # buffer for fire line pixels (radius) to intersect fire perimeter (deg), ~500m
extbuffer = 1000   # buffer to define interior/exterior region, 1000 m
area_VI = 0.141    # km2, area of each 375m VIIRS pixel

# MODIS pixel size
MCD64buf = 231.7   # MODIS fire perimeter buffer (deg), corresponding to 463.31271653 m/2

# ------------------------------------------------------------------------------
# fire type related parameters
# ------------------------------------------------------------------------------

# fire type options
FTYP_opt = 1      # 0: preset ftype for all fires;
                  # 1: use CA type classifications
                  # 2: proposed global fire types
CONT_opt = 1      # 0: preset continuity threshold for all fires;
                  # 1: use CA type classifications dependent values
                  # 2: use global fire types and size dependent values

# For FTYP_opt = 0, use 'forest' for all fires
FTYP_preset = [2,'Forest']  # the default ftype
CONT_preset = 1   # km, default continuity threshold

# For FTYP_opt = 1, use algorithm in CAFEDS
FTYP_CA = {
         0:'Other',
         1:'Urban',
         2:'Forest wild',
         3:'Forest manage',
         4:'Shrub wild',
         5:'Shrub manage',
         6:'Agriculture'}           # fire type names
FTYPCLR_CA = {
         0:'grey',
         1:'rosybrown',
         2:'darkolivegreen',
         3:'olive',
         4:'saddlebrown',
         5:'sandybrown',
         6:'darkviolet'}             # colors used for each fire type
CONT_CA = {
         0:1, 1:1, 2:2.5, 3:5, 4:5, 5:5, 6:1
         }  # preset fire type dependent CONNECTIVITY_THRESHOLD_KM

# For FTYP_opt = 2, use algorithm proposed for global study
FTYP_Glb = {
         0:'Other',
         1:'Temp Forest',
         2:'Trop Forest',
         3:'Bore Forest',
         4:'Savana',
         5:'Agriculture',
         6:'Deforestation'}           # fire type names

""" FireConsts
This is the module containing all constants used in this project
"""

# project directories
#import os
dirhome = 'D:/fire_atlas'   # home directory
dirpjcode = dirhome + '/1_code/'     # project code directory
dirpjdata = dirhome + '/2_pipeline/'     # output data directory
dirextdata = dirhome + '/Data/' # exterior input data directory
mcd64dir = dirextdata + 'burned_grid/'
lakedir = 'D:/fire_atlas/Data/GlobalSurfaceWater/vector/'

# extents of test regions
ext_SL = [-119, 60, -110, 64] # test area for slave lake 2014
ext_NT = [-136, 60, -102, 79] # Northwest territories
ext_AKsmall = [-156, 65.5, -152, 66.3] # small test area for AK 2015
ext_AK = [-168, 51, -140, 72] # all of AK
ext_NY = [147, 65, 154, 70] # test area for Northern Yakutia 2020
ext_YAK = [106, 55, 163, 72] # test area for yakutia 2020
ext_ssib = [89, 48, 132.1, 67] # test area southern siberia 2018 reference data
ext_all = {'AK2015':ext_AKsmall, 'AK':ext_AK, 'SL2014':ext_SL, 'NT':ext_NT,
           'NY2020':ext_NY, 'YAK':ext_YAK, 'Sib2018':ext_ssib}

# parameters used for fire pixel clustering
EARTH_RADIUS_KM = 6371.0  # earth radius, km
SPATIAL_THRESHOLD_KM = 1  # threshold of spatial distance (between pixels, km) used to do fire initial clustering
#MAX_THRESH_CENT_KM = 50   # threshold of spatial distance (between clusters centroid, km) used to filter nearby clusters

# fire type and visualization parameters
lc_dict = { # dictionary for simplifying ESA CCI land cover types
           0:0, 20:0, 21:0, 22:0, # no data, water, ice, bare
           1:1, 2:1, 3:1, 4:1, # cropland
           5:2, 6:3, 7:4, 8:5, 9:6, # forest
           10:7, 11:7, 12:7, # shrubs & mosaic landscape
           13:8, 14:9, 15:10, # grassland, lichen, sparse vegetation
           16:11, 17:11, 18:12, # flooded forests and shrubs
           19:13} # urban
# FTYP = {0:'Other', 1:'Crop', 2:'broad-evergreen', 3:'broad-deciduous', 4:'needle-evergreen', 5:'needle-deciduous', 
#         6:'mixed forest', 7:'shrubs/mosaic', 8:'grassland', 9:'lichen/moss', 10:'sparse vegetation',
#         11:'flood-forest', 12:'flood-shrubs',13:'urban'} 
FTYP = {0:'Other', 1:'Crop', 2:'forest-peat', 3:'forest-no peat', 4:'shrubs/mosaic/open-peat', 
        5:'shrubs/mosaic/open-no peat', 6:'grassland:peat', 7:'grassland-no peat', 8:'urban'} 
# FTYP = {0:'Other', 1:'Urban', 2:'Forest wild', 3:'Forest manage', 4:'Shrub wild', 5:'Shrub manage', 6:'Agriculture'}      # fire type names
# FTYPCLR = {0:'grey', 1:'rosybrown', 2:'darkolivegreen', 3:'olive', 4:'saddlebrown', 5:'sandybrown', 6:'darkviolet'}        # colors used for each fire type

# temporal and spatial distances for fire object definition
maxoffdays = {0:4, 1:5, 2:10, 3:10, 4:10, 5:5, 6:10, 7:3, 8:5}   # fire becomes inactive after this number of consecutive days without active fire detection
limoffdays = 20 # this is the additional threshold that involves checking for active fireline
CONNECTIVITY_THRESHOLD_KM = {0:1, 1:1, 2:2.5, 3:2.5, 4:2.5, 5:2.5, 6:1, 7:1, 8:1} # km, corresponding to fire type (forest/shrub: 5km; other: 1km)
sleeperthresh = 1 #km, distance a sleeper is allowes to spread since its last active fire detection

# alpha parameter
valpha = 800      # 1000 m (default)
lalpha = 100
stralpha = '1km'   # string of alpha parameter

# lengths and areas related to VIIRS pixl data
VIIRSbuf = 187.5   # fire perimeter buffer (deg), corresponding to 375m/2 at lat=30
fpbuffer = 200     # buffer use to determine fire line pixels (deg), ~200m
flbuffer = 500     # buffer for fire line pixels (radius) to intersect fire perimeter (deg), ~500m
extbuffer = 1000   # buffer to define interior/exterior region, 1000 m
area_VI = 0.141    # km2, area of each 375m VIIRS pixel

# mcd bffer (not used)
MCD64buf = 231.7   # fire perimeter buffer (deg), corresponding to 463.31271653 m/2

# parameter used for determine large fires
falim = 4    # large fire area threshold (in km2)

# diagnostic data name and types saved in geojson files (in addition to geometries)
dd = {'fireid':'int',                  # id
      'mergid':'int',               # this is the id in the large fire database
      #'clat':'float',               # centroid latitude   -> centroid[0]
      #'clon':'float',               # centroid longitude  -> centroid[1]
      'ftype':'int',                # fire type
      'isactive':'int',             # active status
      't_inactive':'int',           # how long has it been inactive
      'isignition':'int',           # is this a new ignition?
      'invalid':'int',              # invalid status
      'n_pixels':'int',             # number of total pixels
      'n_newpixels':'int',          # number of new pixels
      'farea':'float',              # fire size
      'fperim':'float',             # fire perimeter length
      'flinelen':'float',           # active fire front line length
      'duration':'float',           # fire duration
      'pixden':'float',             # fire pixel density
      'meanFRP':'float',            # mean FRP of the new fire pixels
      'tst_year':'int',             # t_st[0]
      'tst_month':'int',
      'tst_day':'int',
      'tst_ampm':'str',
      'ted_year':'int',             # t_ed[0]
      'ted_month':'int',
      'ted_day':'int',
      'ted_ampm':'str',
      'ted_doy':'int',
      }

# dd = {'fid':'int',                  # id
#       'clat':'float',               # centroid latitude   -> centroid[0]
#       'clon':'float',               # centroid longitude  -> centroid[1]
#       'ftype':'int',                # fire type
#       'isactive':'int',             # active status
#       'invalid':'int',              # invalid status
#       'n_pixels':'int',             # number of total pixels
#       'n_newpixels':'int',          # number of new pixels
#       'farea':'float',              # fire size
#       'fperim':'float',             # fire perimeter length
#       'flinelen':'float',           # active fire front line length
#       'duration':'float',           # fire duration
#       'pixden':'float',             # fire pixel density
#       'meanFRP':'float',            # mean FRP of the new fire pixels
#       'tst_year':'int',             # t_st[0]
#       'tst_month':'int',
#       'tst_day':'int',
#       'tst_ampm':'str',
#       'ted_year':'int',             # t_ed[0]
#       'ted_month':'int',
#       'ted_day':'int',
#       'ted_ampm':'str',
#       }
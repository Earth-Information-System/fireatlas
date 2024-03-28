from enum import Enum


class FireSource(Enum):
    VIIRS = "VIIRS"
    SNPP = "SNPP"
    NOAA20 = "NOAA20"


class EPSG(Enum):
    CONUS_EQ_AREA = 9311
    GLOBAL_EQ_AREA = 6933
    HI_LAT = 3571
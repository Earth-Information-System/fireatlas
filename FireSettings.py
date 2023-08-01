from typing import Optional
from enum import Enum
import pydantic


class FireSource(Enum):
    VIIRS = "VIIRS"
    SNPP = "SNPP"
    NOAA20 = "NOAA20"
    
class EPSG(Enum):
    CONUS_EQ_AREA = 9311
    GLOBAL_EQ_AREA = 6933
    HI_LAT = 3571

class RuntimeSettings(pydantic.BaseSettings):
    
    number_of_multi_proc_workers: int = 4

    epsg: EPSG.CONUS_EQ_AREA

    ftyp_opt: int = 1

    cont_opt: int = 1

    firesrc: FireSource.VIIRS
    
    firenrt: bool = True # NRT - True, False
    
    class Config:
        env_prefix = "FIRE_"
        env_file = ".env"
        

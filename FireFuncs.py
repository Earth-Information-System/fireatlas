from FireLog import logger
from utils import timed

def get_CONNECTIVITY_FIRE(fire):
    """ set the CONNECTIVITY_FIRE_KM value (km) for a given fire. Within this
    buffer of the fire exterior, new pixels will be appended to the fire.

    Constants
    ---------
    CONT_opt : int
        option number of CONNECTIVITY_FIRE_KM calculation
    CONT_preset : float
        the preset CONNECTIVITY_FIRE_KM value for option 0
    CONT_CA : dict
        the preset CONNECTIVITY_FIRE_KM lookup table for option 1
    FTYP_Glb : dict
        the preset CONNECTIVITY_FIRE_KM lookup table for option 2

    Parameters
    ----------
    fire : fire object
        the fire associated with the CONNECTIVITY_FIRE_KM

    Returns
    -------
    v : float
        the CONNECTIVITY_FIRE_KM (km)
    """
    from FireConsts import CONT_opt

    if CONT_opt == 0:  # use preset values
        from FireConsts import CONT_preset

        return CONT_preset
    elif CONT_opt == 1:  # use lookup table from CAFEDS
        from FireConsts import CONT_CA

        return CONT_CA[fire.ftype]
    elif CONT_opt == 2:  # use lookup table for global universal run
        from FireConsts import FTYP_Glb

        ftype = FTYP_Glb[fire.ftype]
        
        fnpix = min(fire.n_pixels, 25)  # set an upper limit 
        
        if ftype == "Temp Forest":
            v = fnpix * 2 / 25 + 0.8  # 0.8 (0) - 2.8 (25)
        elif ftype == "Trop Forest":
            v = fnpix * 0.7 / 25 + 0.7  # 0.7 (0) - 1.4 (25)
        elif ftype == "Bore Forest":
            v = fnpix * 3.2 / 25 + 1.0  # 1.0 (0) - 4.2 (25)
        elif ftype == "Savana":
            v = fnpix * 2.5 / 25 + 0.9  # 0.9 (0) - 3.4 (25)
        else:
            v = 0.7
        return v


def set_ftype(fire):
    """ set fire type and dominant LCT for newly activated fires

    Parameters
    ----------
    fire : fire object
        the fire associated with the CONNECTIVITY_FIRE_KM

    """
    import numpy as np
    from FireConsts import FTYP_opt
    import FireIO
    from datetime import date
    import random

    # 0 - use preset ftype (defined as FTYP_preset in FireConsts) for all fires
    # 1 - use the CA type classification (dtermined using LCTmax)
    # 2 - use the global type classfication (need more work...)

    if FTYP_opt == 0:  # use preset ftype for all fires
        from FireConsts import FTYP_preset

        ftype = FTYP_preset[0]
    elif FTYP_opt == 1:  # use CA type classifications (determined using LCTmax)
        # update or read LCTmax; calculated using all newlocs
        # TODO: not sure why we need the `else` branch here for sampling?!
        # so just increasing this so we don't have to sample 1000 points
        if fire.n_newpixels < 10000000:
            uselocs = fire.newlocs_geo
        else:
            # we can do a random sample of 1000 new pixels (it's likely going to be a forest fire anyways)
            uselocs = random.sample(fire.newlocs_geo, 1000)

        # get all LCT for the fire pixels
        vLCT = FireIO.get_LCT_CONUS(uselocs)
        try:
            # extract the LCT with most pixel counts
            LCTmax = max(set(vLCT), key=vLCT.count)
        except:
            logger.info('No LCT data available, setting ftype to 0...')
            ftype = 0
            return ftype

        # determine the fire type using the land cover type and stFM1000
        if LCTmax in [0, 11, 31]:  #'NoData', 'Water', 'Barren' -> 'Other'
            ftype = 0
        elif LCTmax in [23]:  # 'Urban' -> 'Urban'
            ftype = 1
        elif LCTmax in [82]:  # 'Agriculture' -> 'Agriculture'
            ftype = 6
        elif LCTmax in [42]:  # 'Forest' ->
            if fire.stFM1000 > 12:  # 'Forest manage'
                ftype = 3
            else:  # 'Forest wild'
                ftype = 2
        elif LCTmax in [52, 71]:  # 'Shrub', 'Grassland' ->
            if fire.stFM1000 > 12:  # 'Shrub manage'
                ftype = 5
            else:  # 'Shrub wild'
                ftype = 4
        else:
            logger.info(f"Unknown land cover type {LCTmax}. Setting ftype to 0.")
            ftype = 0
    elif FTYP_opt == 2:  # global type classification
         # update or read LCTmax; calculated using all newlocs
        if fire.n_newpixels < 1000:
            uselocs = fire.newlocs_geo
        else:
            # we can do a random sample of 1000 new pixels (it's likely going to be a forest fire anyways)
            uselocs = random.sample(fire.newlocs_geo, 1000)
        
        vLCT = FireIO.get_LCT_Global(
            uselocs
        )  # call get_LCT to get all LCT for the fire pixels
        
        try:
            LCTmax = max(
            set(vLCT), key=vLCT.count)  # extract the LCT with most pixel counts
        except:
            logger.info('No LCT data available, setting ftype to 0...')
            ftype = 0
            return ftype
        
        #print("LCTMAX:",LCTmax)
        if LCTmax in [0, 50, 60, 70, 80, 90, 100, 200]:
            ftype = 0
        # ^^^ current catch-all for 'Other'. 
        # See: https://developers.google.com/earth-engine/datasets/catalog/COPERNICUS_Landcover_100m_Proba-V-C3_Global
        
        elif LCTmax in [20]: # Shrub --> Savanna
            ftype = 4
        
        elif LCTmax in [40]: # Agriculture class
            ftype = 5
        
        else: # begin forested class
            lat_list = [xy[1] for xy in uselocs]
            lat_mean = abs(np.nanmean(lat_list))
            # Forest Classifications based on: https://ucmp.berkeley.edu/exhibits/biomes/forests.php
            if lat_mean < 23.5:
                ftype = 2 # tropical
            elif lat_mean < 50:
                ftype = 1 # temperate
            else:
                ftype = 3 # boreal
    return ftype


def set_ftypename(fire):
    """ return Fire type name

    Parameters
    ----------
    fire : fire object
        the fire associated with the CONNECTIVITY_FIRE_KM

    Returns
    -------
    ftname : str
        the fire type name of the given fire

    """
    from FireConsts import FTYP_opt

    if FTYP_opt == 0:  # use preset ftype for all fires
        from FireConsts import FTYP_preset

        ftname = FTYP_preset[0]
    elif FTYP_opt == 1:  # use CA type classifications (determined using LCTmax)
        from FireConsts import FTYP_CA

        ftname = FTYP_CA[fire.ftype]
    elif FTYP_opt == 2:  # global type classification
        from FireConsts import FTYP_Glb

        ftname = FTYP_Glb[fire.ftype]
    return ftname

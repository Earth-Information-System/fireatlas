import numpy as np


def get_CONNECTIVITY_FIRE(fire):
    """set the CONNECTIVITY_FIRE_KM value (km) for a given fire. Within this
    buffer of the fire exterior, new pixels will be appended to the fire.

    Constants
    ---------
    CONT_OPT : str
        continuity threshold option
    CONT : float
        the lookup table for CONNECTIVITY_FIRE_KM values
    FTYP : dict
        the lookup table for fire type names

    Parameters
    ----------
    fire : fire object
        the fire associated with the CONNECTIVITY_FIRE_KM

    Returns
    -------
    v : float
        the CONNECTIVITY_FIRE_KM (km)
    """
    from fireatlas import FireConsts, settings

    if settings.CONT_OPT in ["preset", "CA"]:
        return FireConsts.CONT[settings.CONT_OPT][fire.ftype]
    elif settings.CONT_OPT == "global":
        ftype_name = FireConsts.FTYP["global"][fire.ftype]

        fnpix = min(fire.n_pixels, 25)  # set an upper limit

        if ftype_name == "Temp Forest":
            v = fnpix * 2 / 25 + 0.8  # 0.8 (0) - 2.8 (25)
        elif ftype_name == "Trop Forest":
            v = fnpix * 0.7 / 25 + 0.7  # 0.7 (0) - 1.4 (25)
        elif ftype_name == "Bore Forest":
            v = fnpix * 3.2 / 25 + 1.0  # 1.0 (0) - 4.2 (25)
        elif ftype_name == "Savana":
            v = fnpix * 2.5 / 25 + 0.9  # 0.9 (0) - 3.4 (25)
        else:
            v = 0.7
        return v


def set_ftype(fire):
    """set fire type and dominant LCT for newly activated fires

    Parameters
    ----------
    fire : fire object
        the fire associated with the CONNECTIVITY_FIRE_KM

    """
    from fireatlas import FireConsts, FireIO, settings
    from fireatlas.FireLog import logger

    # "preset" - use preset ftype (defined as FTYP in FireConsts) for all fires
    # "CA" - use the CA type classification (determined using LCTmax)
    # "global" - use the global type classification (need more work...)

    if settings.FTYP_OPT == "preset":  # use preset ftype for all fires
        ftype = FireConsts.FTYP_preset
    elif settings.FTYP_OPT == "CA":  
        # update or read LCTmax; calculated using all newlocs
        if fire.n_newpixels < 1000:
            uselocs = fire.newlocs_geo
        else:
            # we can do a random sample of 1000 new pixels (it's likely going to be a forest fire anyways)
            uselocs = fire.newlocs_geo[
                np.random.choice(fire.newlocs_geo.shape[0], 1000, replace=False), :
            ]

        # get all LCT for the fire pixels
        vLCT = FireIO.get_LCT_CONUS(uselocs)
        try:
            # extract the LCT with most pixel counts
            LCTmax = max(set(vLCT), key=vLCT.count)
        except:
            logger.info("No LCT data available, setting ftype to 0...")
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
    elif settings.FTYP_OPT == "global":
        # update or read LCTmax; calculated using all newlocs
        if fire.n_newpixels < 1000:
            uselocs = fire.newlocs_geo
        else:
            # we can do a random sample of 1000 new pixels (it's likely going to be a forest fire anyways)
            uselocs = fire.newlocs_geo[
                np.random.choice(fire.newlocs_geo.shape[0], 1000, replace=False), :
            ]

        vLCT = FireIO.get_LCT_Global(
            uselocs
        )  # call get_LCT to get all LCT for the fire pixels

        try:
            LCTmax = max(
                set(vLCT), key=vLCT.count
            )  # extract the LCT with most pixel counts
        except:
            logger.info("No LCT data available, setting ftype to 0...")
            ftype = 0
            return ftype

        # print("LCTMAX:",LCTmax)
        if LCTmax in [0, 50, 60, 70, 80, 90, 100, 200]:
            ftype = 0
        # ^^^ current catch-all for 'Other'.
        # See: https://developers.google.com/earth-engine/datasets/catalog/COPERNICUS_Landcover_100m_Proba-V-C3_Global

        elif LCTmax in [20]:  # Shrub --> Savanna
            ftype = 4

        elif LCTmax in [40]:  # Agriculture class
            ftype = 5

        else:  # begin forested class
            lat_list = [xy[1] for xy in uselocs]
            lat_mean = abs(np.nanmean(lat_list))
            # Forest Classifications based on: https://ucmp.berkeley.edu/exhibits/biomes/forests.php
            if lat_mean < 23.5:
                ftype = 2  # tropical
            elif lat_mean < 50:
                ftype = 1  # temperate
            else:
                ftype = 3  # boreal
    return ftype


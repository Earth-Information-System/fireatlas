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
        fnpix = fire.n_pixels
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


def get_CONNECTIVITY_CLUSTER():
    """ get the connectivity spatial threshold for initial clustering

    Returns:
    --------
    CONNECTIVITY_CLUSTER_KM : float
        the connectivity spatial threshold for initial clustering, km
    """
    CONNECTIVITY_CLUSTER_KM = 0.7
    return CONNECTIVITY_CLUSTER_KM


def get_CONNECTIVITY_SLEEPER():
    """ get the CONNECTIVITY_SLEEPER_KM value

    Returns:
    --------
    CONNECTIVITY_SLEEPER_KM : float
        the connectivity spatial threshold (to previous fire line), km
    """
    CONNECTIVITY_SLEEPER_KM = 1  # 1km
    return CONNECTIVITY_SLEEPER_KM


def set_ftype(fire):
    """ set fire type and dominant LCT for newly activated fires

    Parameters
    ----------
    fire : fire object
        the fire associated with the CONNECTIVITY_FIRE_KM

    """
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
        if fire.n_newpixels < 1000:
            uselocs = fire.newlocs_geo
        else:
            # we can do a random sample of 1000 new pixels (it's likely going to be a forest fire anyways)
            uselocs = random.sample(fire.newlocs_geo, 1000)

        vLCT = FireIO.get_LCT_CONUS(
            uselocs
        )  # call get_LCT to get all LCT for the fire pixels
        try:
            LCTmax = max(
            set(vLCT), key=vLCT.count)  # extract the LCT with most pixel counts
        except:
            print('No LCT data available, setting ftype to 0...')
            ftype = 0
            return ftype
        # get and record fm1000 value at ignition
        ignct = fire.ignlocsMP.centroid  # ignition centroid
        loc = (ignct.y, ignct.x)  # (lat,lon)
        t = date(*fire.t_st[:-1])
        try: stFM1000 = FireIO.get_FM1000(t, loc)
        except:
            #print('FM1000 data is unavailable at this time.')
            stFM1000 = 0

        # determine the fire type using the land cover type and stFM1000
        if LCTmax in [0, 11, 31]:  #'NoData', 'Water', 'Barren' -> 'Other'
            ftype = 0
        elif LCTmax in [23]:  # 'Urban' -> 'Urban'
            ftype = 1
        elif LCTmax in [82]:  # 'Agriculture' -> 'Agriculture'
            ftype = 6
        elif LCTmax in [42]:  # 'Forest' ->
            if stFM1000 > 12:  # 'Forest manage'
                ftype = 3
            else:  # 'Forest wild'
                ftype = 2
        elif LCTmax in [52, 71]:  # 'Shrub', 'Grassland' ->
            if stFM1000 > 12:  # 'Shrub manage'
                ftype = 5
            else:  # 'Shrub wild'
                ftype = 4
        else:
            print(f"Unknown land cover type {LCTmax}. Setting ftype to 0.")
            ftype = 0
    elif FTYP_opt == 2:  # global type classification
       
        # TODO @ Kat: call to GCT function
        # TODO for runs, modify FTYP_opt logic
        
        # -------
        
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
            print('No LCT data available, setting ftype to 0...')
            ftype = 0
            return ftype
        # get and record fm1000 value at ignition
        ignct = fire.ignlocsMP.centroid  # ignition centroid
        loc = (ignct.y, ignct.x)  # (lat,lon)
        t = date(*fire.t_st[:-1])
        try: stFM1000 = FireIO.get_FM1000(t, loc)
        except:
            #print('FM1000 data is unavailable at this time.')
            stFM1000 = 0

        # determine the fire type using the land cover type and stFM1000
        if LCTmax in [0, 11, 31]:  #'NoData', 'Water', 'Barren' -> 'Other'
            ftype = 0
        elif LCTmax in [23]:  # 'Urban' -> 'Urban'
            ftype = 1
        elif LCTmax in [82]:  # 'Agriculture' -> 'Agriculture'
            ftype = 6
        elif LCTmax in [42]:  # 'Forest' ->
            if stFM1000 > 12:  # 'Forest manage'
                ftype = 3
            else:  # 'Forest wild'
                ftype = 2
        elif LCTmax in [52, 71]:  # 'Shrub', 'Grassland' ->
            if stFM1000 > 12:  # 'Shrub manage'
                ftype = 5
            else:  # 'Shrub wild'
                ftype = 4
        else:
            print(f"Unknown land cover type {LCTmax}. Setting ftype to 0.")
            ftype = 0
        
        # ftype = 1  # need more work here...

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

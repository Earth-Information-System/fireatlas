""" FireGdf_sfs
Module for creating geojson summary for temporal evolution of each single fire
        (large active fires  up to the present time step)
Large fires are fires with area > falim

List of geojson types
---------------------
* fperim(''): fire basic attributes and perimeter geometry
* fline('FL'): active fire line geometry
* NFP('NFP'): new fire pixels

List of functions
-----------------
* make_fire_history
* merge_fires
* save_gdf_1fire
* save_gdf_trng

Modules required
----------------
* FireObj
* FireIO
* FireConsts
"""


def find_largefires(allfires, falim=4):
    """ Given an allfires object, extract large fires to be used for single fire
    recording. Only active fires or sleepers are considered.

    Parameters
    ----------
    allfires : Allfires obj
        the allfires object at a given time step
    falim : float
        the area threshold (km) for large fire

    Returns
    -------
    large_ids : list
        the list of large fire IDs
    """
    large_ids = [
        fid
        for fid in (allfires.fids_active + allfires.fids_sleeper)
        if allfires.fires[fid].farea > falim
    ]
    return large_ids


def find_mergefires_ct(allfires, fid):
    """ Find fires that are merged to fire (fid) at time t and started before t

    Parameters
    ----------
    allfires : Allfires obj
        the allfires object at a given time step
    fid : int
        the target fire ID
    t : tuple (y,m,d,ampm)
        the target merging time

    Returns
    -------
    fids_m : list
        the list of fire IDs that are merged to fire (fid) at time t
    """
    t = allfires.t
    # fids_m = list(sorted([h[0] for h in allfires.heritages if
    #             ((allfires.fires[h[0]].t_ed == t) & (allfires.fires[h[0]].t_st != t)
    #              & (h[1] == fid))]))
    fids_m = []
    for h in allfires.heritages:
        if h[1] == fid:
            if h[0] in allfires.fids:
                f_m = allfires.fires[h[0]]
                if (f_m.t_ed == t) & (f_m.t_st != t):
                    fids_m.append(h[0])
    fids_m = list(sorted(fids_m))
    return fids_m


def merge_fires(gdf_1d, fid, t, sfkeys):
    """ merge fire perimeters of all fires with the same merge id
    Parameters
    ----------b
    gdf : geopandas DataFrame
        the gdf containing all entries that should be merged
    fid: the mergeid of the fire
    Returns
    -------
    gdf_diss : geopandas DataFrame
        the gdf containing a merged geometry and summary stats
    """
    from datetime import date
    import FireObj, FireTime

    # dissolve the dataframe
    gdf_1d["mergeid"] = fid
    gdf_diss = gdf_1d.dissolve(by="mergeid")
    # gdf_sf.loc[t_m,'geometry'] = gdf_sf.loc[t_m,'geometry'].union(gs_t_m['geometry'])

    # replace some values with sums
    if "n_pixels" in sfkeys:
        gdf_diss.loc[fid, "n_pixels"] = sum(gdf_1d.n_pixels)
    if "n_newpixels" in sfkeys:
        gdf_diss.loc[fid, "n_newpixels"] = sum(gdf_1d.n_newpixels)
    if "farea" in sfkeys:
        gdf_diss.loc[fid, "farea"] = sum(gdf_1d.farea)
    if "fperim" in sfkeys:
        gdf_diss.loc[fid, "fperim"] = sum(gdf_1d.fperim)
    if "flinelen" in sfkeys:
        gdf_diss.loc[fid, "flinelen"] = sum(gdf_1d.flinelen)
    if "t" in sfkeys:
        gdf_diss.loc[fid, "t"] = FireTime.t2dt(t)

    # weighted average computed for averages
    if ("pixden" in sfkeys) & ("farea" in sfkeys):
        gdf_1d = gdf_1d.assign(
            pixweigh=gdf_1d["pixden"] * gdf_1d["farea"]
        )  ## IS THIS CORRECT?
        gdf_diss.loc[fid, "pixden"] = sum(gdf_1d.pixweigh) / sum(gdf_1d.farea)
    if ("meanFRP" in sfkeys) & ("n_pixels" in sfkeys):
        gdf_1d = gdf_1d.assign(FRPweigh=gdf_1d["meanFRP"] * gdf_1d["n_pixels"])
        gdf_diss.loc[fid, "meanFRP"] = sum(gdf_1d.FRPweigh) / sum(gdf_1d.n_pixels)

    # here need some column dropping...
    gdf_diss = gdf_diss.reset_index().drop(columns="mergeid")

    return gdf_diss


def make_sf(t, regnm, layer, fids_m, fid, sfkeys):
    """ read snapshot gpkg files, extract records for 1 or more fires, merge them,
    and create gdf ready for making single fire time series

    Parameters
    ----------
    t : tuple, (int,int,int,str)
        the year, month, day and 'AM'|'PM' of time
    regnm : str
        the study region name (used for reading previous day's data)
    fids_m : list
        the list of fire IDs needed to merge (1-element list is ok)
    sfkeys : list
        the attributes (geometry excluded) used for making single fire time series

    Returns
    -------
    gdf_1d : Geopandas DataFrame
        the 1-row gdf for the single fire
    """
    import FireIO, FireTime
    import pandas as pd
    import geopandas as gpd
    from FireConsts import epsg

    if layer != "nfplist":
        # extract rows for fires merged to the target fire
        gdf = FireIO.load_gpkgobj(
            t, regnm, layer=layer
        )  # read daily gdf: active and sleeper only

        # extract rows with fireID in fids_m only
        idx = [i for i in fids_m if i in list(gdf.index)]
        if len(idx) == 0:  # if no rows satistying condition, return None
            return None

        gdf_1d = gdf.loc[idx]  # extract all fires merged to fid
        gdf_1d = gdf_1d[
            sfkeys + ["geometry"]
        ]  # drop some columns not used in single fire gpkg

        gdf_1d = gdf_1d.reset_index().drop(columns="fireID")

        # if there are multiple rows (fires), merge them
        if len(gdf_1d) > 1:
            gdf_1d = merge_fires(gdf_1d, fid, t, sfkeys)

    # # change index to date; maybe no need to drop FireID?
    # gdf_1d.index = [FireTime.t2dt(t)]

    else:  # 'nfplist' layer uses a different approach from fire object directly
        # extract rows for fires merged to the target fire
        allfires = FireIO.load_fobj(t, regnm, activeonly=True)

        # loop over all fids_m to make nfplist at t
        gdf_1d = None
        for fid in fids_m:
            if fid in allfires.fids_valid:
                f = allfires.fires[fid]
                # gdf_1f = make_nfplist_gdf(allfires.fires[fid])
                plist = []
                for p in f.newpixels:
                    p_attrs = []
                    for k in sfkeys:
                        p_attrs.append(getattr(p, k))
                    plist.append(p_attrs)

                if len(plist) > 0:
                    df = pd.DataFrame(plist)
                    df.columns = sfkeys
                    gdf_1f = gpd.GeoDataFrame(
                        df,
                        crs="epsg:" + str(epsg),
                        geometry=gpd.points_from_xy(df.x, df.y),
                    )

                    if gdf_1d is None:
                        gdf_1d = gdf_1f
                    else:
                        gdf_1d = gdf_1d.append(gdf_1f)

    return gdf_1d


def make_sfts_1f(allfires, fid, regnm, dd, layer="perimeter"):
    """ create the single large fire gpkg file at time t.
    In this approximate approach, the snapshots for fid between t_st and t_ed are
    extracted and concategated. In addition, the snapshots for fires merged to
    this fire (with heritage target as fid) are also added (except for the t_ed
    when the fire with fid already contains the contributions from these merged
    fires).

    Parameters
    ----------
    allfires : Allfires obj
        the allfires object at a given time step
    fid : int
        the target fire ID
    regnm : str
        the study region name (used for reading previous day's data)

    Returns
    -------
    gdf_all : Geopandas DataFrame
        the full time series gdf for the single fire
    """

    import FireTime, FireIO

    # single fire attributes list from input dd
    sfkeys = list(dd.keys())

    # extract the fire
    f = allfires.fires[fid]

    # find all fires merged to this fire at any time step
    fids_m = list(sorted([h[0] for h in allfires.heritages if h[1] == fid]))

    # loop over all merged fires
    endloop = False  # flag to control the ending of the loop
    t = list(f.t_st)  # t is the time (year,month,day,ampm) for each step
    while endloop == False:

        gdf_1d = make_sf(t, regnm, layer, fids_m + [fid], fid, sfkeys)

        # append daily row to gdf_all
        if FireTime.t_dif(t, f.t_st) == 0:
            gdf_all = gdf_1d
        else:
            gdf_all = gdf_all.append(gdf_1d, ignore_index=True)

        #  - if t reaches ted, set endloop to True to stop the loop
        if FireTime.t_dif(t, f.t_ed) == 0:
            endloop = True

        #  - update t with the next time stamp
        t = FireTime.t_nb(t, nb="next")

    return gdf_all


def update_sfts_1f(allfires, fid, regnm, layer="perimeter"):
    """ update the single large fire gpkg file at time t

    Parameters
    ----------
    allfires : Allfires obj
        the allfires object at a given time step
    fid : int
        the target fire ID
    regnm : str
        the study region name (used for reading previous day's data)

    Returns
    -------
    gdf_all : Geopandas DataFrame
        the full time series gdf for the single fire
    """
    import FireTime, FireIO
    import geopandas as gpd

    # the target single large fire and target time
    f = allfires.fires[fid]
    t = allfires.t

    # # attributes name and data  types (in addition to geometries)
    if layer == "perimeter":
        dd = {
            "n_pixels": "int",  # number of total pixels
            "n_newpixels": "int",  # number of new pixels
            "farea": "float",  # fire size
            "fperim": "float",  # fire perimeter length
            "flinelen": "float",  # active fire front line length
            "duration": "float",  # fire duration
            "pixden": "float",  # fire pixel density
            "meanFRP": "float",  # mean FRP of the new fire pixels
            "t": "datetime64",
        }
    elif layer == "fireline":
        dd = {
            "t": "datetime64",
        }
    elif layer == "newfirepix":
        dd = {
            "t": "datetime64",
        }
    elif layer == "nfplist":
        dd = {
            "x": "float",
            "y": "float",
            "frp": "float",
            "DS": "float",
            "DT": "float",
            "ampm": "str",
            "datetime": "datetime64",
            "sat": "str",
        }
    sfkeys = list(dd.keys())

    # try to read corresponding sf gpkg file at previous time step (gdf_sf_pt)
    t_pt = FireTime.t_nb(t, nb="previous")
    gdf_sf_pt = FireIO.load_gpkgsfs(t_pt, fid, regnm, layer=layer)

    # if no gdf_sf_pt, create time series using the make_fire_history()
    # if the running is fast, no need to use the code in the 'else' part...
    if gdf_sf_pt is None:
        gdf_all = make_sfts_1f(allfires, fid, regnm, dd, layer=layer)
    # when gdf_sf_pt is present, to save time, read it and add fires just merged at t
    else:
        # use gdf_sf_pt as basis
        gdf_all = gdf_sf_pt

        # find all fires merged to this fire at present time step
        fids_m = find_mergefires_ct(allfires, fid)

        # read allfires object at previous time step
        allfires_pt = FireIO.load_fobj(t_pt, regnm, activeonly=True)

        # loop over all merged fires and add to gdf_all
        for fid_m in fids_m:
            # create time series for merged fire at previous time step using make_fire_history()
            gdf_all_m = make_sfts_1f(allfires_pt, fid_m, regnm, dd, layer=layer)

            # combine gdf_sf_pt and gdf_all_m
            gdf_all = gdf_all.append(gdf_all_m, ignore_index=True)

        # also add current fire record at present time time step (for fid at t only)
        gdf_ct = make_sf(t, regnm, layer, [fid], fid, sfkeys)
        # for k,tp in dd.items():
        #     gdf_ct[k] = gdf_ct[k].astype(tp)
        gdf_all = gdf_all.append(gdf_ct, ignore_index=True)

    # 4. force the correct dtypes
    for k, tp in dd.items():
        gdf_all[k] = gdf_all[k].astype(tp)

    # gdf_all = gdf_all.reset_index()
    return gdf_all


def save_sfts_1f(
    fid,
    regnm,
    allfires=None,
    t=None,
    layers=["perimeter", "fireline", "newfirepix", "nfplist"],
):
    """Wrapper to create and save gpkg files at a time step for a large fire

    Parameters
    ----------
    fid : int
        the fire ID
    regnm : str
        region name
    allfires : Allfires object
        the object containing the fire
    t : tuple, (int,int,int,str)
        the year, month, day and 'AM'|'PM' of time
    """
    import FireIO

    # if allfires is not given, read it using t
    if allfires == None:
        allfires = FireIO.load_fobj(t, regnm, activeonly=True)

    # create and save fperim, fline, nfp time series
    if "perimeter" in layers:
        gdf_fperim = update_sfts_1f(allfires, fid, regnm, layer="perimeter")
        FireIO.save_gpkgsfs(allfires.t, fid, regnm, gdf_fperim=gdf_fperim)

    if "fireline" in layers:
        gdf_fline = update_sfts_1f(allfires, fid, regnm, layer="fireline")
        FireIO.save_gpkgsfs(allfires.t, fid, regnm, gdf_fline=gdf_fline)

    if "newfirepix" in layers:
        gdf_nfp = update_sfts_1f(allfires, fid, regnm, layer="newfirepix")
        FireIO.save_gpkgsfs(allfires.t, fid, regnm, gdf_nfp=gdf_nfp)

    if "nfplist" in layers:
        gdf_nfplist = update_sfts_1f(allfires, fid, regnm, layer="nfplist")
        FireIO.save_gpkgsfs(allfires.t, fid, regnm, gdf_nfplist=gdf_nfplist)


def save_sfts_all(t, regnm, layers=["perimeter", "fireline", "newfirepix", "nfplist"]):
    """Wrapper to create and save gpkg files at a time step for all large fires

    Parameters
    ----------
    t : tuple, (int,int,int,str)
        the year, month, day and 'AM'|'PM' of time
    regnm : str
        region name
    """
    import FireIO

    # read allfires object
    allfires = FireIO.load_fobj(t, regnm, activeonly=True)

    # find all large active fires and sleepers
    large_ids = find_largefires(allfires)

    # loop over all fires, and create/save gpkg files for each single fire
    for fid in large_ids:
        save_sfts_1f(fid, regnm, allfires=allfires, layers=layers)


def save_sfts_trng(
    tst, ted, regnm, layers=["perimeter", "fireline", "newfirepix", "nfplist"]
):
    """Wrapper to create and save all gpkg files for a time period

    Parameters
    ----------
    tst : tuple, (int,int,int,str)
        the year, month, day and 'AM'|'PM' at start time
    ted : tuple, (int,int,int,str)
        the year, month, day and 'AM'|'PM' at end time
    regnm : str
        region name
    """
    import FireTime

    # loop over all days during the period
    endloop = False  # flag to control the ending olf the loop
    t = list(tst)  # t is the time (year,month,day,ampm) for each step
    while endloop == False:
        print("Single fire saving", t)

        # create and save all gpkg files at time t
        save_sfts_all(t, regnm, layers=layers)

        # time flow control
        #  - if t reaches ted, set endloop to True to stop the loop
        if FireTime.t_dif(t, ted) == 0:
            endloop = True

        #  - update t with the next time stamp
        t = FireTime.t_nb(t, nb="next")


if __name__ == "__main__":
    """ The main code to record time series of geojson data for a fire
    """

    import time

    t1 = time.time()

    # set the start and end time
    tst = (2020, 9, 5, "AM")
    # ted=(2020,9,10,'PM')
    ted = (2020, 11, 5, "PM")

    # run creation
    save_sfts_trng(tst, ted, regnm="Creek")
    # save_sfts([2020, 9, 19, 'AM'],'Creek')

    t2 = time.time()
    print(f"{(t2-t1)/60.} minutes used to run code")

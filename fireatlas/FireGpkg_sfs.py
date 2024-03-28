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
import os
import warnings
import multiprocessing
import geopandas as gpd
from AT import AT_FRAP
from glob import glob
import pandas as pd
import time

warnings.simplefilter(action='ignore', category=FutureWarning)


def getdd(layer):
    ''' Get attributes names and formats for different gpkg layers
    '''
    # # attributes name and data  types (in addition to geometries)
    if layer == "all":
        dd = {
            "mergeid": "int",  # this is the id in the large fire database
            "invalid": "bool",
            "ftype": "int",  # fire type
            "n_pixels": "int",  # number of total pixels
            "n_newpixels": "int",  # number of new pixels
            "farea": "float",  # fire size
            "fperim": "float",  # fire perimeter length
            "flinelen": "float",  # active fire front line length
            "duration": "float",  # fire duration
            "pixden": "float",  # fire pixel density
            "meanFRP": "float",  # mean FRP of the new fire pixels
            "t_st": "datetime64[ns]",
            "t_ed": "datetime64[ns]",
            "hull": "geometry",
            "fline": "geometry",
            "nfp": "geometry",
        }
    elif layer == "perimeter":
        dd = {
            "n_pixels": "int",  # number of total pixels
            "n_newpixels": "int",  # number of new pixels
            "farea": "float",  # fire size
            "fperim": "float",  # fire perimeter length
            "flinelen": "float",  # active fire front line length
            "duration": "float",  # fire duration
            "pixden": "float",  # fire pixel density
            "meanFRP": "float",  # mean FRP of the new fire pixels
            "t": "datetime64[ns]",
        }
    elif layer == "fireline":
        dd = {
            "t": "datetime64[ns]",
        }
    elif layer == "newfirepix":
        dd = {
            "t": "datetime64[ns]",
        }
    elif layer == "nfplist":
        dd = {
            "x": "float",
            "y": "float",
            "frp": "float",
            "DS": "float",
            "DT": "float",
            "ampm": "str",
            "datetime": "datetime64[ns]",
            "sat": "str",
        }

    return dd

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
    from .FireTime import t2dt

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
        gdf_diss.loc[fid, "t"] = t2dt(t)

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


def make_sf(t, regnm, layer, fids_m, fid):
    """ At a given time step, create single row DataFrame for a target fire fid
        - read snapshot gpkg files for the time step (t, regnm, layer)
        - extract records (sfkeys) for target fire and fires merged to it (fids_m)
        - merge the extracted records
        - create gdf ready for making single fire time series
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
    from .FireIO import load_gpkgobj

    dd = getdd(layer)
    sfkeys = list(dd.keys())

    # extract rows for fires merged to the target fire
    gdf = load_gpkgobj(
        t, regnm, layer=layer
    )  # read daily gdf: active and sleeper only

    # extract rows with fireID in fids_m only
    if gdf is None:
        return None

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

    # # make sure the fireID is correct (this is needed to merge all fires)
    # gdf_1d['fireID'] = fid

    # # change index to date; maybe no need to drop FireID?
    # gdf_1d.index = [FireTime.t2dt(t)]

    return gdf_1d

def make_sf_nfplist(allfires, t, regnm, fids):
    """ At a given time step, create single row DataFrame for all newly detected
        pixels associated with fires fids.
    Parameters
    ----------
    t : tuple, (int,int,int,str)
        the year, month, day and 'AM'|'PM' of time
    regnm : str
        the study region name (used for reading previous day's data)
    fids : list
        the list of fire IDs needed to merge (1-element list is ok)
    sfkeys : list
        the attributes (geometry excluded) used for making single fire time series
    Returns
    -------
    gdf_1d : Geopandas DataFrame
        the 1-row gdf for the single fire
    """
    from .FireIO import load_fobj
    from .FireTime import t2dt
    from .FireConsts import epsg

    dd = getdd('nfplist')
    sfkeys = list(dd.keys())

    if t != allfires.t:
        allfires = load_fobj(t, regnm, activeonly=False)
    # record all nfps for all fires with fids to make nfplist at t
    gdf_1d = None
    for fid in fids:
        if fid in allfires.fids_valid:   # only record newpixels for valid fires
            f = allfires.fires[fid]
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

                # also add t variable (detection 12-hourly time step)
                gdf_1f['t'] = t2dt(t)

                if gdf_1d is None:
                    gdf_1d = gdf_1f
                else:
                    gdf_1d = gdf_1d.append(gdf_1f)
    del allfires
    return gdf_1d

def make_sfts_1f(allfires,f, fid, fids_m, regnm, layer="perimeter"):
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
    from .FireTime import t_dif, t_nb


    # loop over all merged fires
    endloop = False  # flag to control the ending of the loop
    t = list(f.t_st)  # t is the time (year,month,day,ampm) for each step
    gdf_all = None
    while endloop == False:

        if layer == 'nfplist':
            gdf_1d = make_sf_nfplist(allfires, t, regnm, fids_m+[fid])
        else:
            gdf_1d = make_sf(t, regnm, layer, fids_m+[fid], fid)

        # append daily row to gdf_all
        # if FireTime.t_dif(t, f.t_st) == 0:
        if (gdf_1d is not None):
            if (gdf_all is None):
                gdf_all = gdf_1d
            else:
                gdf_all = gdf_all.append(gdf_1d, ignore_index=True)

        #  - if t reaches ted, set endloop to True to stop the loop
        if t_dif(t, f.t_ed) == 0:
            endloop = True

        #  - update t with the next time stamp
        t = t_nb(t, nb="next")

    return gdf_all



def update_sfts_1f(allfires, allfires_pt, fid, regnm, layer="perimeter"):
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
    from .FireIO import load_gpkgsfs

    # the target single large fire and target time
    f = allfires.fires[fid]
    t = allfires.t

    dd = getdd(layer)
    sfkeys = list(dd.keys())

    # try to read small fire file at previous time step (gdf_sf_pt)
    t_pt = allfires_pt.t
    gdf_sf_pt = load_gpkgsfs(t_pt, fid, regnm, layer=layer)

    # if no gdf_sf_pt, create historical time series using the make_sfts_1f()
    # if the running is fast enough, no need to use the code in the 'else' part...
    if gdf_sf_pt is None:

        # find all fires merged to this fire at any time step
        fids_m = list(sorted([h[0] for h in allfires.heritages if h[1] == fid]))
        gdf_all = make_sfts_1f(allfires,f,fid, fids_m, regnm, layer=layer)

    # when gdf_sf_pt is present, to save time, read it and add fires just merged at t
    else:
        # use gdf_sf_pt as basis, and only need to add changes at t and newly merged fires
        gdf_all = gdf_sf_pt

        # find all fires merged to this fire at present time step
        
        fids_m = find_mergefires_ct(allfires, fid)

        # loop over all merged fires and add to gdf_all
        for fid_m in fids_m:
            # create historical time series of all newly merged fires at t_pt
            # using the same approach as the case of gdf_sf_pt==None
            f_pt = allfires_pt.fires[fid_m]
            fids_mm = list(sorted([h[0] for h in allfires_pt.heritages if h[1] == fid_m]))
            gdf_all_m = make_sfts_1f(allfires_pt,f_pt, fid_m, fids_mm, regnm, layer=layer)

            # combine all existing gdf_all_m
            gdf_all = gdf_all.append(gdf_all_m, ignore_index=True)

        # also add current fire record at present time time step (for fid at t only)
        if layer == 'nfplist':
            gdf_ct = make_sf_nfplist(allfires, t, regnm, [fid])
        else:
            gdf_ct = make_sf(t, regnm, layer, [fid], fid)

        # for k,tp in dd.items():
        #     gdf_ct[k] = gdf_ct[k].astype(tp)
        gdf_all = gdf_all.append(gdf_ct, ignore_index=True)
    
    if gdf_all is None:
        print('Warning gdf_all is None. Returning...')
        return gdf_all
    
    # 4. force the correct dtypes
    for k, tp in dd.items():
        gdf_all[k] = gdf_all[k].astype(tp)
    if layer == 'nfplist':
        gdf_all['t'] = gdf_all['t'].astype('datetime64')
    
#     ampm = t[-1]
#     if ampm == 'AM':
#         time = pd.to_datetime(str(t[0])+'-'+str(t[1])+'-'+str(t[2])+'T00:00:00')
#     else: 
#         time = pd.to_datetime(str(t[0])+'-'+str(t[1])+'-'+str(t[2])+'T12:00:00')
    
    # set region col
    gdf_all['region'] = str(regnm)
    
    # set primary key col
    if layer == 'nfplist':
        gdf_all['primarykey'] = gdf_all['region'] + '|' + str(fid) + '|' + gdf_all['t'].apply(lambda x: x.isoformat()) + '|' + gdf_all.index.map(str)
    else: 
        gdf_all['primarykey'] = gdf_all['region'] + '|' + str(fid) + '|' + gdf_all['t'].apply(lambda x: x.isoformat())

    # gdf_all = gdf_all.reset_index()
    return gdf_all


def save_sfts_all(queue: multiprocessing.Queue, t, regnm, layers=["perimeter", "fireline", "newfirepix", "nfplist"]):
    """Wrapper to create and save gpkg files at a time step for all large fires
    Parameters
    ----------
    t : tuple, (int,int,int,str)
        the year, month, day and 'AM'|'PM' of time
    regnm : str
        region name
    """
    from .FireLog import logger
    from .FireTime import t_nb
    from .FireIO import load_fobj
    from .FireObj import Allfires

    tstart = time.time()
    # read allfires object
    logger.info(f'Load allfires..')
    allfires = load_fobj(t, regnm, activeonly=False)
    
    t_pt = t_nb(t, nb="previous")
    
    try:
    # read allfires object at previous time step
        logger.info(f'Load allfires_pt...')
        allfires_pt = load_fobj(t_pt, regnm, activeonly=False)
    except: 
        allfires_pt = Allfires(t_pt)
    # find all large active fires and sleepers
    logger.info(f'Finding largefires...')
    large_ids = find_largefires(allfires)

    # loop over all fires, and create/save gpkg files for each single fire
    for fid in large_ids:
        logger.info(f'Adding {fid} to queue')
        queue.put([allfires, allfires_pt, fid, regnm, layers])

    tend = time.time()
    logger.info(f'Time sending to queue for timestep {t} with multiproc: {(tend-tstart)/60.} minutes')


def worker_save_sfts_1f(queue: multiprocessing.Queue):
    """pop off the queue and call `save_sfts_1f`

    :param queue:
    :return:
    """
    from .FireLog import logger

    while True:
        payload = queue.get()
        # check for poison pill
        if payload is None:
            break
        # unpack
        allfires, allfires_pt, fid, regnm, layers = payload
        logger.info(f'[ WORKER {os.getpid()} ]: writing out LF={fid}')
        save_sfts_1f(allfires, allfires_pt, fid, regnm, layers)


def save_sfts_1f(allfires, allfires_pt, fid, regnm, layers=["perimeter", "fireline", "newfirepix", "nfplist"]):
    """Wrapper to create and save gpkg files at a time step for all large fires
    Parameters
    ----------
    t : tuple, (int,int,int,str)
        the year, month, day and 'AM'|'PM' of time
    regnm : str
        region name
    """
    from .FireLog import logger
    from .FireIO import save_gpkgsfs

    logger.info(f'Generating LF data for fid {fid}')
    tstart = time.time()
    f = allfires.fires[fid]
    
    if "perimeter" in layers:
        
        tstart = time.time()
        gdf_fperim = update_sfts_1f(allfires, allfires_pt, fid, regnm, layer="perimeter")
        tend = time.time()
        logger.info(f"{(tend-tstart)/60.} minutes used for updating perimeter layer.")

        tstart = time.time()
        save_gpkgsfs(allfires.t, fid, regnm, gdf_fperim=gdf_fperim)
        tend = time.time()
        logger.info(f"{(tend-tstart)/60.} minutes used for I/O data perimeter layer.")

    if "fireline" in layers:
        gdf_fline = update_sfts_1f(allfires,allfires_pt, fid, regnm, layer="fireline")
        save_gpkgsfs(allfires.t, fid, regnm, gdf_fline=gdf_fline)

    if "newfirepix" in layers:
        gdf_nfp = update_sfts_1f(allfires,allfires_pt, fid, regnm, layer="newfirepix")
        save_gpkgsfs(allfires.t, fid, regnm, gdf_nfp=gdf_nfp)

    if "nfplist" in layers:
        gdf_nfplist = update_sfts_1f(allfires,allfires_pt, fid, regnm, layer="nfplist")
        save_gpkgsfs(allfires.t, fid, regnm, gdf_nfplist=gdf_nfplist)

    tend = time.time()

    logger.info(f"{(tend-tstart)/60.} minutes used to save Largefire data for fid {fid}.")
        
        
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
    from .FireConsts import number_of_multi_proc_workers
    from .FireTime import t_dif, t_nb
    from .FireLog import logger

    # loop over all days during the period
    endloop = False  # flag to control the ending olf the loop
    t = list(tst)  # t is the time (year,month,day,ampm) for each step

    queue = multiprocessing.Queue()
    workers = []
    for i in range(number_of_multi_proc_workers):
        p = multiprocessing.Process(target=worker_save_sfts_1f, args=(queue,))
        p.start()
        workers.append(p)

    while endloop == False:
        print("Single fire saving", t)
        logger.info('Single fire saving: '+str(t))
        
        tstart = time.time()
        # create and save all gpkg files at time t
        try:
            save_sfts_all(queue, t, regnm, layers=layers)
        except Exception as exc:
            logger.exception(exc)
        tend = time.time()
        logger.info(f"{(tend-tstart)/60.} minutes used to save Largefire data for this time.")

        if t_dif(t, ted) == 0:
            endloop = True

        #  - update t with the next time stamp
        t = t_nb(t, nb="next")

    # add poison pill to stop work
    for i in range(number_of_multi_proc_workers):
        queue.put(None)

    # wait all writer workers
    for p in workers:
        p.join()

def combine_sfts(regnm,yr,addFRAP=False):
    ''' Combine all large fire time series to a gpkg data, with the option to add fire name from FRAP data
    used for a large region
    Parameters
    ----------
    regnm : str
        the name of the region
    yr : int
        year
    addFRAP : bool
        if set to true, add fire name/id from FRAP database (require FRAP data and AT_FRAP function)
    '''
    from .FireIO import gpd_read_file, get_gpkgsfs_dir
    from .FireTime import dt_dif

    strdir = get_gpkgsfs_dir(yr,regnm)

    if addFRAP:
        vFRAP = AT_FRAP.getFRAPyr(yr)
        if vFRAP is not None:
            vFRAP.geometry = vFRAP.geometry.buffer(0)

    # get a list of all fids
    fids = []
    for f in glob(strdir+'/*.gpkg'):
        fid = int(os.path.basename(f).split('.')[0].split('_')[0][1:])
        fids.append(fid)

    # loop over all fids
    gdf_perim_all, gdf_fline_all, gdf_nfp_all, gdf_nfplist_all = None, None, None, None
    for fid in set(fids):
        print(fid,end=',')
        # for each fire, find the last gpkg file
        fnmlast = sorted(glob(strdir+'/F'+str(fid)+'_*.gpkg'))[-1]

        # read the data (four layers)
        gdf_perim = gpd_read_file(fnmlast,layer='perimeter').to_crs(epsg=4326)
        gdf_fline = gpd_read_file(fnmlast,layer='fireline').to_crs(epsg=4326)
        gdf_nfp = gpd_read_file(fnmlast,layer='newfirepix').to_crs(epsg=4326)
        gdf_nfplist = gpd_read_file(fnmlast,layer='nfplist').to_crs(epsg=4326)

        # add fire ID to the gdf
        gdf_perim['FireID'] = fid
        gdf_fline['FireID'] = fid
        gdf_nfp['FireID'] = fid
        gdf_nfplist['FireID'] = fid

        # add duration to fline, nfp, and nfplist
        t0 = pd.to_datetime(gdf_perim.iloc[0].t)
        gdf_fline['duration'] = pd.to_datetime(gdf_fline['t']).apply(lambda x: dt_dif(x,t0))
        gdf_nfp['duration'] = pd.to_datetime(gdf_nfp['t']).apply(lambda x: dt_dif(x,t0))
        gdf_nfplist['duration'] = pd.to_datetime(gdf_nfplist['t']).apply(lambda x: dt_dif(x,t0))

        # add FRAP fire name and id
        if addFRAP:
            if vFRAP is not None:
                gFEDS =  gdf_perim.iloc[-1].geometry
                FRAPid,FRAPfnm = AT_FRAP.get_FRAPidfnm_from_geometry(vFRAP,gFEDS,fid)
                if FRAPid is not None:
                    gdf_perim['FRAPid'] = FRAPid
                    gdf_perim['FRAPfnm'] = FRAPfnm
                    gdf_fline['FRAPid'] = FRAPid
                    gdf_fline['FRAPfnm'] = FRAPfnm
                    gdf_nfp['FRAPid'] = FRAPid
                    gdf_nfp['FRAPfnm'] = FRAPfnm
                    gdf_nfplist['FRAPid'] = FRAPid
                    gdf_nfplist['FRAPfnm'] = FRAPfnm
                else:
                    gdf_perim['FRAPid'] = -1
                    gdf_perim['FRAPfnm'] = ''
                    gdf_fline['FRAPid'] = -1
                    gdf_fline['FRAPfnm'] = ''
                    gdf_nfp['FRAPid'] = -1
                    gdf_nfp['FRAPfnm'] = ''
                    gdf_nfplist['FRAPid'] = -1
                    gdf_nfplist['FRAPfnm'] = ''

        # combine all fids
        if (gdf_perim_all is None):
            gdf_perim_all = gdf_perim
            gdf_fline_all = gdf_fline
            gdf_nfp_all = gdf_nfp
            gdf_nfplist_all = gdf_nfplist
        else:
            gdf_perim_all = gdf_perim_all.append(gdf_perim,ignore_index=True)
            gdf_fline_all = gdf_fline_all.append(gdf_fline,ignore_index=True)
            gdf_nfp_all = gdf_nfp_all.append(gdf_nfp,ignore_index=True)
            gdf_nfplist_all = gdf_nfplist_all.append(gdf_nfplist,ignore_index=True)

    # drop x and y columns for nfplist
    gdf_nfplist_all = gdf_nfplist_all.drop(columns=['x','y'])

    # save combined gdfs
    fnmout = os.path.join(strdir,'largefires.gpkg')
    gdf_perim_all.to_file(fnmout, driver="GPKG", layer="perimeter")
    gdf_fline_all.to_file(fnmout, driver="GPKG", layer="fireline")
    gdf_nfp_all.to_file(fnmout, driver="GPKG", layer="newfirepix")
    gdf_nfplist_all.to_file(fnmout, driver="GPKG", layer="nfplist")
#     return gdf_perim_all, gdf_fline_all, gdf_nfp_all, gdf_nfplist_all

def convert_sfts(regnm,yr,fids):
    ''' Convert a large fire time series to a gpkg data; usually used to process a single fire run
    Parameters
    ----------
    regnm : str
        the name of the region
    yr : int
        year
    fids : list of int
        fire ids
    '''
    from .FireIO import get_gpkgsfs_dir, gpd_read_file
    from .FireTime import dt_dif

    strdir = get_gpkgsfs_dir(yr,regnm)
    
    # loop over all fids
    gdf_perim_all, gdf_fline_all, gdf_nfp_all, gdf_nfplist_all = None, None, None, None
    for fid in set(fids):
        print(fid,end=',')
        # for each fire, find the last gpkg file
        fnmlast = sorted(glob(strdir+'/F'+str(fid)+'_*.gpkg'))[-1]

        # read the data (four layers)
        gdf_perim = gpd_read_file(fnmlast,layer='perimeter').to_crs(epsg=4326)
        gdf_fline = gpd_read_file(fnmlast,layer='fireline').to_crs(epsg=4326)
        gdf_nfp = gpd_read_file(fnmlast,layer='newfirepix').to_crs(epsg=4326)
        gdf_nfplist = gpd_read_file(fnmlast,layer='nfplist').to_crs(epsg=4326)

        # add fire ID to the gdf
        gdf_perim['FireID'] = fid
        gdf_fline['FireID'] = fid
        gdf_nfp['FireID'] = fid
        gdf_nfplist['FireID'] = fid

        # add duration to fline, nfp, and nfplist
        t0 = pd.to_datetime(gdf_perim.iloc[0].t)
        gdf_fline['duration'] = pd.to_datetime(gdf_fline['t']).apply(lambda x: dt_dif(x,t0))
        gdf_nfp['duration'] = pd.to_datetime(gdf_nfp['t']).apply(lambda x: dt_dif(x,t0))
        gdf_nfplist['duration'] = pd.to_datetime(gdf_nfplist['t']).apply(lambda x: dt_dif(x,t0))

        # combine all fids
        if (gdf_perim_all is None):
            gdf_perim_all = gdf_perim
            gdf_fline_all = gdf_fline
            gdf_nfp_all = gdf_nfp
            gdf_nfplist_all = gdf_nfplist
        else:
            gdf_perim_all = gdf_perim_all.append(gdf_perim,ignore_index=True)
            gdf_fline_all = gdf_fline_all.append(gdf_fline,ignore_index=True)
            gdf_nfp_all = gdf_nfp_all.append(gdf_nfp,ignore_index=True)
            gdf_nfplist_all = gdf_nfplist_all.append(gdf_nfplist,ignore_index=True)

    # drop x and y columns for nfplist
    gdf_nfplist_all = gdf_nfplist_all.drop(columns=['x','y'])

    # save combined gdfs
    fnmout = os.path.join(strdir,'largefires.gpkg')
    gdf_perim_all.to_file(fnmout, driver="GPKG", layer="perimeter")
    gdf_fline_all.to_file(fnmout, driver="GPKG", layer="fireline")
    gdf_nfp_all.to_file(fnmout, driver="GPKG", layer="newfirepix")
    gdf_nfplist_all.to_file(fnmout, driver="GPKG", layer="nfplist")

if __name__ == "__main__":
    """ The main code to record time series of geojson data for a fire
    """
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

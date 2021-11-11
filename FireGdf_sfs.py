""" FireGdf_sfs
Module for creating geojson summary for temporal evolution of each single fire
        (large active fires  up to the present time step)
Large fires are fires with area > falim (defined in FireConsts)

List of geojson types
---------------------
* fperim(''): fire basic attributes and perimeter geometry
* fline('FL'): active fire line geometry
* NFP('NFP'): new fire pixels

List of functions
-----------------
* make_fire_history
* update_fire_history
* save_gdf_1t
* save_gdf_trng
* yrend_clean

Modules required
----------------
* FireObj
* FireIO
* FireConsts
"""

def make_fire_history(fire,op=''):
    ''' derive time series of single fire attributes using the half-daily gdf summary

    Parameters
    ----------
    fire : Fire object
        the fire need to update

    Returns
    -------
    gdf_all : geopandas DataFrame
        the gdf containing half daily fire basic attributes and fire perimeter
    '''
    import pandas as pd
    import FireObj,FireIO

    fid = fire.id
    tst,ted = fire.t_st,fire.t_ed  # start and ending time of the fire

    endloop = False  # flag to control the ending of the loop
    t = list(tst)    # t is the time (year,month,day,ampm) for each step
    while endloop == False:
        # print(t)
        # read daily gdf
        gdf = FireIO.load_gdfobj(t,op=op)
        gdf_1d = gdf[gdf.index == fid]

        if len(gdf_1d) == 0:
            gdf_1d.loc[fid] = None

        # append daily row to gdf_all
        if FireObj.t_dif(t,tst)==0:
            gdf_all = gdf_1d
        else:
            gdf_all = gdf_all.append(gdf_1d)

        #  - if t reaches ted, set endloop to True to stop the loop
        if FireObj.t_dif(t,ted)==0:
            endloop = True

        #  - update t with the next time stamp
        t = FireObj.t_nb(t,nb='next')


    # use correct time column as index
    gdf_all.index = FireObj.ftrange(tst,ted)

    # remove fid column to save space
    # gdf_all = gdf_all.drop(columns='fid')

    return gdf_all

def update_fire_history(fire,op=''):
    ''' update time series of a single fire history

    Parameters
    ----------
    fire : Fire object
        the fire need to update

    Returns
    -------
    gdf_sf_1d : geopandas DataFrame
        the gdf containing half daily fire basic attributes and fire perimeter
    '''
    import pandas as pd
    import FireObj,FireIO
    import os

    fid = fire.id
    t = fire.t

    # get gdf_sfs for previous time step
    pt = FireObj.t_nb(t,nb='previous')
    fnm = FireIO.get_gdfobj_sf_fnm(pt, fid, op=op)

    # if prevous t exist, read current summary gdf and append to previous data (to save time)
    if os.path.exists(fnm):
        # read daily gdf
        gdf = FireIO.load_gdfobj(t,op=op)
        gdf_1d = gdf[gdf.index == fid]

        # gdf_1d = FireIO.load_gdfobj(t,op=op).loc[[fid]]
        gdf_1d.index = [FireObj.t2dt(t)]
        # gdf_1d = gdf_1d.drop(columns='fid')

        # read prevous single fire time series
        gdf_sf_pt = FireIO.load_gdfobj_sf(pt,fid,op=op)

        # append the daily gdf to single fire time series
        gdf_sf_1d = gdf_sf_pt.append(gdf_1d)

    # if no prevous t file, create a new time series gdf
    else:
        # call make_fire_history func
        gdf_sf_1d = make_fire_history(fire,op=op)

    return gdf_sf_1d

def save_gdf_1t(t,op=''):
    ''' derive and save time series of fire perimeter and fine line
        for each large active fire at a time step

    Parameters
    ----------
    t : tuple, (int,int,int,str)
        the year, month, day and 'AM'|'PM' at the time
    '''
    import FireIO
    from FireConsts import falim

    # read allfires
    allfires = FireIO.load_fobj(t)

    # get fids for large active fires
    gdf = FireIO.load_gdfobj(t)
    fids_laf = list(gdf[(gdf.farea>falim) & (gdf.isactive==1)].index)

    # create and save gdfs for each large active fire
    if len(fids_laf) > 0:
        # print(f'{fids_laf} large active fires')
        for fid in fids_laf:
            fire = allfires.fires[fid]
            # create gdf
            gdf_ts = update_fire_history(fire,op=op)

            # save gdf
            FireIO.save_gdfobj_sf(gdf_ts,t,fid,op=op)

def save_gdf_trng(tst,ted,fperim=False,fline=False,NFP=False,fall=False):
    ''' Wrapper to create and save gdf files for a time period

    Parameters
    ----------
    tst : tuple, (int,int,int,str)
        the year, month, day and 'AM'|'PM' at start time
    ted : tuple, (int,int,int,str)
        the year, month, day and 'AM'|'PM' at end time
    fperim : bool
        if set to true, save fire perimeters
    fline : bool
        if set to true, save fire lines
    NFP : bool
        if set to true, save new fire pixels
    fall : bool
        if set to true, save all
    '''
    import FireObj, FireIO

    # if fall is True, set all options True
    if fall:
        fperim,fline,NFP = True,True,True

    # loop over all days during the period
    endloop = False  # flag to control the ending of the loop
    t = list(tst)    # t is the time (year,month,day,ampm) for each step
    while endloop == False:
        # print(t)

        # create and save gdfs according to input options
        if fperim:
            save_gdf_1t(t)
        if fline:
            save_gdf_1t(t,op='FL')
        if NFP:
            save_gdf_1t(t,op='NFP')

        # time flow control
        #  - if t reaches ted, set endloop to True to stop the loop
        if FireObj.t_dif(t,ted)==0:
            endloop = True

        #  - update t with the next time stamp
        t = FireObj.t_nb(t,nb='next')

def yrend_clean(tyed):
    ''' At the end of a year, extract single fire files at t_ed of the fire
    '''
    import FireIO
    from FireConsts import falim, dirpjdata
    from datetime import date
    # import shutil
    import os

    # get all single fire ids using the year end data
    # tyed = (year,12,31,'PM')
    year = tyed[0]
    gdf = FireIO.load_gdfobj(tyed,op='')
    allfires = FireIO.load_fobj(tyed)
    fid_sfs = gdf[(gdf['farea']>falim) & (gdf.invalid==0)].index.to_list()

    # for each fire, only save the file at t_ed of the fire,
    for fid in fid_sfs:
        fire = allfires.fires[fid]
        t = fire.t_ed
        d = date(*t[:-1])
        for op in ['','FL','NFP']:

            # fnmin is the single large fire data recorded at t_ed
            fnmin = FireIO.get_gdfobj_sf_fnm(t,fid,op=op)

            # all single large fires data for fid in the year
            fnms = FireIO.get_gdfobj_sf_fnms_year(year,fid,op=op)

            # only keeps the latest file for a fire
            for fnm in fnms:
                if not os.path.samefile(fnm, fnmin):
                    # print(fnm)
                    os.remove(fnm)

    # also remove individual fire data for invalid fires
    fids_term = gdf[gdf.invalid==1].index.to_list()
    for fid in fids_term:
        for op in ['','FL','NFP']:
            fnms = FireIO.get_gdfobj_sf_fnms_year(year,fid,op=op)
            for fnm in fnms:
                os.remove(fnm)

if __name__ == "__main__":
    ''' The main code to record time series of geojson data for a fire
    '''

    import time
    t1 = time.time()
    # set the start and end time
    tst=(2019,1,1,'AM')
    ted=(2019,12,31,'PM')

    # for each day during the period,
    # create and save geojson files for temporal evolution of large fires
    # save_gdf_trng(tst=tst,ted=ted,fall=True)

    yrend_clean(2019)

    t2 = time.time()
    print(f'{(t2-t1)/60.} minutes used to run code')

    # # only run at year end to clean the files
    # for year in range(2012,2020):
    #     print(year)
    #     yrend_clean(year)

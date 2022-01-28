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
* merge_fires
* save_gdf_1fire
* save_gdf_trng

Modules required
----------------
* FireObj
* FireIO
* FireConsts
"""

def make_fire_history(fid, start_end,op=''):
    ''' derive time series of single fire attributes using the half-daily gdf summary

    Parameters
    ----------
    fire : Fire object
        the fire need to update
    start_end: tuple
        start and end date and AM/PM of the fire object

    Returns
    -------
    gdf_all : geopandas DataFrame
        the gdf containing half daily fire basic attributes and fire perimeter
    '''
    import math
    import FireObj,FireIO

    #fid = fire.id
    tst,ted = start_end  # start and ending time of the fire

    endloop = False  # flag to control the ending of the loop
    t_inact = 0      # counts the days a fire is inactive before it spreads again
    t = list(tst)    # t is the time (year,month,day,ampm) for each step
    while endloop == False:
        #print(t)
        # read daily gdf
        gdf = FireIO.load_gdfobj(t,op=op)
        gdf_fid = gdf[gdf.mergid == fid] # here we need to enter merge id instead of index
        
        # skip date if now new pixels
        gdf_fid = gdf_fid[gdf_fid.n_newpixels > 0]
        if len(gdf_fid) == 0:
            t = FireObj.t_nb(t,nb='next')
            t_inact += 0.5 # half-daily time steps
            continue
        
        # merge if several ids
        if len(gdf_fid) > 1:
            gdf_fid = merge_fires(gdf_fid, fid, tst, t)
        
        # change index to date
        gdf_fid.index = [FireObj.t2dt(t)]
        
        # update t_inactive
        gdf_fid.t_inactive = math.floor(t_inact) # inactivity is counted in full days

        # append daily row to gdf_all
        if FireObj.t_dif(t,tst)==0:
            gdf_all = gdf_fid
        else:
            gdf_all = gdf_all.append(gdf_fid)

        #  - if t reaches ted, set endloop to True to stop the loop
        if FireObj.t_dif(t,ted)==0:
            endloop = True

        #  - update t with the next time stamp
        t_inact = 0 # reset inactivity counter
        t = FireObj.t_nb(t,nb='next')


    # use correct time column as index
    # gdf_all.index = FireObj.ftrange(tst,ted)

    # remove fid column to save space
    # gdf_all = gdf_all.drop(columns='fid')

    return gdf_all


def merge_fires(gdf_fid, fid, tst, t):
    ''' merge fire perimeters of all fires with the same merge id

    Parameters
    ----------
    gdf : geopandas DataFrame
        the gdf containing all entries that should be merged
    fid: the mergeid of the fire

    Returns
    -------
    gdf_diss : geopandas DataFrame
        the gdf containing a merged geometry and summary stats
    '''
    from datetime import date
    import FireObj
    
    # dissolve the dataframe
    gdf_diss = gdf_fid.dissolve(by = 'mergid')
    
    # replace some values with sums
    gdf_diss['mergid'] = fid
    gdf_diss.loc[fid,'n_pixels'] = sum(gdf_fid.n_pixels)
    gdf_diss.loc[fid,'n_newpixels'] = sum(gdf_fid.n_newpixels)
    gdf_diss.loc[fid,'farea'] = sum(gdf_fid.farea)
    gdf_diss.loc[fid,'fperim'] = sum(gdf_fid.fperim)
    gdf_diss.loc[fid,'flinelen'] = sum(gdf_fid.flinelen)
    gdf_diss.loc[fid,'tst_year'] = int(tst[0])
    gdf_diss.loc[fid,'tst_month'] = int(tst[1])
    gdf_diss.loc[fid,'tst_day'] = int(tst[2])
    gdf_diss.loc[fid,'tst_ampm'] = tst[3]
    gdf_diss.loc[fid,'ted_year'] = int(t[0])
    gdf_diss.loc[fid,'ted_month'] = int(t[0])
    gdf_diss.loc[fid,'ted_day'] = int(t[0])
    gdf_diss.loc[fid,'ted_ampm'] = t[3]
    gdf_diss.loc[fid,'ted_doy'] = date(*t[:-1]).timetuple().tm_yday
    gdf_diss.loc[fid,'duration'] = FireObj.t_dif(tst,t) + 0.5
    
    # weighted average computed for averages
    gdf_fid = gdf_fid.assign(pixweigh = gdf_fid['pixden']*gdf_fid['farea']) ## IS THIS CORRECT?
    gdf_diss.loc[fid,'pixden'] = sum(gdf_fid.pixweigh)/sum(gdf_fid.farea)
    gdf_fid = gdf_fid.assign(FRPweigh = gdf_fid['meanFRP']*gdf_fid['n_newpixels'])
    gdf_diss.loc[fid,'meanFRP'] = sum(gdf_fid.FRPweigh)/sum(gdf_fid.n_newpixels)
    
    return(gdf_diss)


def save_gdf_1fire(fid, start_end, op=''):
    ''' derive and save time series of fire perimeter and fine line
        for each large active fire at a time step

    Parameters
    ----------
    fid : int
        the large fire id to be written out
    start_end: dict
        dictionary containing start and end dates for all large fires
    '''
    import FireIO
    start_end_fire = start_end[fid]
    gdf_ts = make_fire_history(fid, start_end_fire, op=op)
    t = start_end_fire[1]
    FireIO.save_gdfobj(gdf_ts,t,param='large',fid=fid,op=op)


def save_gdf_trng(ted,fperim=False,fline=False,NFP=False,fall=False):
    ''' Wrapper to create and save gdf files for a time period

    Parameters
    ----------
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
    import FireIO
    from FireConsts import falim

    # if fall is True, set all options True
    if fall:
        fperim,fline,NFP = True,True,True

    # select all large fires (only merge ids)
    skip,merge_ids = zip(*FireIO.load_fobj(ted).heritages)
    allfires = FireIO.load_fobj(ted)
    large_ids = [allfires.fires[i].id for i in range(allfires.number_of_fires) if allfires.fires[i].farea > falim]
    large_ids = [fire for fire in large_ids if fire not in skip]
    start_end = {i:(allfires.fires[i].t_st, allfires.fires[i].t_ed) for i in large_ids}


    # loop over ids
    for fid in large_ids:
        #print(fid)
        # create and save gdfs according to input options
        if fperim:
            save_gdf_1fire(fid = fid, start_end = start_end)
        if fline:
            save_gdf_1fire(fid = fid, start_end = start_end,op='FL')
        if NFP:
            save_gdf_1fire(fid = fid, start_end = start_end,op='NFP')


if __name__ == "__main__":
    ''' The main code to record time series of geojson data for a fire
    '''

    import time
    t1 = time.time()
    # set the start and end time
    tst=(2020,6,1,'AM')
    ted=(2020,8,31,'PM')

    # for each day during the period,
    # create and save geojson files for temporal evolution of large fires
    save_gdf_trng(ted=ted,fperim=True)

    t2 = time.time()
    print(f'{(t2-t1)/60.} minutes used to run code')

    # # only run at year end to clean the files
    # for year in range(2012,2020):
    #     print(year)
    #     yrend_clean(year)

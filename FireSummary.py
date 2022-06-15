""" FireSummary
Modules for creating a summary of all valid fires up to a time step (usually at year end).

List of summary data
--------------------
* time series of fire statistics
sum_items = {'Number of fires':'int',
             'Number of fire pixels':'int',
             'Fire total area':'float',
             'Fire total perimeter':'float',
             'Fire mean duration':'float',
             'Number of active fires':'int',
             'Number of new fire pixels':'int',
             'Active fire total area':'float',
             'Active fire total perimeter':'float',
             'Active fire line total length':'float',
             'Active fire mean FRP':'float',
             'Active fire mean duration':'float',
             }

* list of large fires
{fid : t_ed}

* list of fire heritage
{fid : [fids_merged]}
"""

sum_items = {'Number of fires':'int',
             'Number of fire pixels':'int',
             'Fire total area':'float',
             'Fire total perimeter':'float',
             'Fire mean duration':'float',
             'Number of active fires':'int',
             'Number of new fire pixels':'int',
             'Active fire total area':'float',
             'Active fire total perimeter':'float',
             'Active fire line total length':'float',
             'Active fire mean FRP':'float',
             'Active fire mean duration':'float',
             }

def cal_sumfires_1t(gdf):
    ''' calculate summary attributes separately for all fires

    Parameters
    ----------
    gdf : gpd DataFrame
        the gdf file with daily fire attributes
    t : datetime date
        date of this day

    Returns
    -------
    summary : pandas Series
        updated summary file
    '''
    
    # create new date columns
    gdf['tst'] = gdf['tst_year'].astype(str) + gdf['tst_month'].astype(str).str.zfill(2) + gdf['tst_day'].astype(str).str.zfill(2) + gdf['tst_ampm']
    gdf['datetime'] = gdf['ted_year'].astype(str) + gdf['ted_month'].astype(str).str.zfill(2) + gdf['ted_day'].astype(str).str.zfill(2) + gdf['ted_ampm']
    gdf.drop(columns = ['isactive','t_inactive','isignition','invalid','n_pixels','tst_year','tst_month',
                        'tst_day','tst_ampm','ted_year','ted_month','ted_day','ted_ampm','ted_doy','lake_border_tot'],
             inplace=True)
    
    return gdf

def cal_summean_1t(gdf,name='All'):
    ''' calculate summary attributes for all or all veg-type

    Parameters
    ----------
    df_summary : pandas DataFrame
        input summary file
    gdf : gpd DataFrame
        the gdf file with daily fire attributes
    t : datetime date
        date of this day
    vtype : str
        vege type, 'All'(default) | other veg types in FTYP

    Returns
    -------
    summary : pandas Series
        updated summary file
    '''
    import pandas as pd

    summary = pd.Series(index = sum_items.keys(),dtype='float')
    gdf_valid = gdf[gdf.invalid==0]        # valid fires
    gdf_active = gdf[gdf.isactive==1]       # active fires

    # non-termninated fires
    summary['Number of fires'] = len(gdf_valid)
    summary['Number of fire pixels'] = gdf_valid.n_pixels.sum()
    summary['Fire total area'] = gdf_valid.farea.sum()
    summary['Fire total perimeter'] = gdf_valid.fperim.sum()
    summary['Fire mean duration'] = gdf_valid.duration.mean()

    # active fires
    summary['Number of active fires'] = len(gdf_active)
    summary['Number of new fire pixels'] = gdf_active.n_newpixels.sum()
    summary['Active fire total area'] = gdf_active.farea.sum()
    summary['Active fire total perimeter'] = gdf_active.fperim.sum()
    summary['Active fire line total length'] = gdf_active.flinelen.sum()
    summary['Active fire mean FRP'] = gdf_active.meanFRP.mean()
    summary['Active fire mean duration'] = gdf_active.duration.mean()

    summary.name = name
    return summary

def makesummary_1d(gdf):
    ''' Calculate daily summary from the geopandas data

    Parameters
    ----------
    gdf : geopandas dataframe
        Saved daily fire geopandas data

    Returns
    -------
    df_summary : pandas dataframe
        Daily summary dataframe for different fire types.
    '''
    import pandas as pd
    from FireConsts import FTYP

    # initialize summary dataframe
    vFTYP = list(FTYP.values())
    iFTYP = list(FTYP.keys())
    rows = ['All'] + vFTYP

    df_summary = pd.DataFrame(cal_summean_1t(gdf))

    # derive summary for variable types of fires
    for i in iFTYP:
        ndx_veg = (gdf['ftype']==i)
        gdf_veg = gdf[ndx_veg]
        ts_summary_veg = cal_summean_1t(gdf_veg,name=vFTYP[i])
        df_summary[vFTYP[i]] = ts_summary_veg

    # rotate and set data types
    df_summary = df_summary.T.astype(sum_items)

    return df_summary

def makesummary_ts(dst,ded,region=''):
    ''' calculate daily summary and update the summary to a given date

    Parameters
    ----------
    dst : tuple (year, month, day, ampm)
        start date
    ded : tuple (year, month, day, ampm)
        end date

    Returns
    -------
    new_summary : xr dataset
        the summary dataset for a time period
    '''
    import FireIO, FireObj
    import xarray as xr

    # loop over all days during the period
    new_summary = None
    gdf_sf = []
    endloop = False  # flag to control the ending of the loop
    t = list(dst)    # t is the time (year,month,day,ampm) for each step
    while endloop == False:
        # print(t)

        # Load daily fire data from geopandas file
        if FireIO.check_gdfobj(t,op='',region=region):
            gdf = FireIO.load_gdfobj(t,region=region,geom=False)

            # call makesummary_1d to get daily summary dataframe
            df_summary_1d = makesummary_1d(gdf)
    
            # convert dataframe to xarray dataset
            xr_sum_1d = df_summary_1d.to_xarray().rename({'index':'Veg'})
    
            # add time coordinate
            dt = FireObj.t2dt(t)
            xr_t = xr_sum_1d.assign_coords(time=dt).expand_dims('time')
    
            # record to new_summary dataset
            if new_summary is None:
                new_summary = xr_t
            else:
                new_summary = xr.concat([new_summary,xr_t],dim='time')
            
            # summary of single fires
            gdf_sf.append(cal_sumfires_1t(gdf)) 

        # time flow control
        #  - if t reaches ted, set endloop to True to stop the loop
        if FireObj.t_dif(t,ded)==0:
            endloop = True

        #  - update t with the next time stamp
        t = FireObj.t_nb(t,nb='next')

    return new_summary, gdf_sf

def updatesummary_ts(pt,ded):
    ''' update the summary by adding summary info to all time steps since dst

    Parameters
    ----------
    pt : tuple (year, month, day, ampm)
        previous time step with summary file
    ded : tuple (year, month, day, ampm)
        end date

    Returns
    -------
    new_summary : xr dataset
        the summary dataset for a time period
    '''
    import FireIO, FireObj
    import xarray as xr

    # first time step after pt
    dst = FireObj.t_nb(pt,nb='next')

    # read summary for previous time step
    new_summary = FireIO.load_summary(pt)

    # loop over all days during the period
    endloop = False  # flag to control the ending of the loop
    t = list(dst)    # t is the time (year,month,day,ampm) for each step
    while endloop == False:
        # print(t)

        # Load daily fire data from geopandas file
        gdf = FireIO.load_gdfobj(t)

        # call makesummary_1d to get daily summary dataframe
        df_summary_1d = makesummary_1d(gdf)

        # convert dataframe to xarray dataset
        xr_sum_1d = df_summary_1d.to_xarray().rename({'index':'Veg'})
        # xr_sum_1d = df_summary_1d.T.to_xarray().rename({'index':'Veg'})

        # add time coordinate
        dt = FireObj.t2dt(t)
        xr_t = xr_sum_1d.assign_coords(time=dt).expand_dims('time')

        # record to new_summary dataset
        new_summary = xr.concat([new_summary,xr_t],dim='time')

        # time flow control
        #  - if t reaches ted, set endloop to True to stop the loop
        if FireObj.t_dif(t,ded)==0:
            endloop = True

        #  - update t with the next time stamp
        t = FireObj.t_nb(t,nb='next')
    return new_summary


def save_sum_1d(tst,ted,region=''):
    ''' Calculate summary file at a time step (starting from the year start)

    Parameters
    ----------
    ted : tuple, (int,int,int,str)
        the year, month, day and 'AM'|'PM' at current time
    '''
    import FireIO, FireObj
    import pandas as pd

    # pt is the previous time step when the summary file is present
    pt = FireIO.get_summary_fnm_lt(ted,region=region)

    # if no summary file is present in this year, create one by recording
    if pt is None:
        # tst = (ted[0],1,1,'AM')
        ds_summary, df_sf = makesummary_ts(tst,ted,region)
        
        df_sf = pd.concat(df_sf)
        FireIO.save_summary_sf(df_sf, ted, region=region)
    # if early summary file is present, update based on that file
    else:
        ds_summary = updatesummary_ts(pt,ted)

    # # add fire heritages up to ted
    # ds_summary = add_heritage(ds_summary,ted)
    #
    # # add list of large fires up to ted
    # ds_summary = add_largefirelist(ds_summary,ted)

    # save the summary file
    FireIO.save_summary(ds_summary,ted,region=region)

def save_sum_trng(tst,ted):
    ''' Wrapper to create and save summary files for a time period

    Parameters
    ----------
    tst : tuple, (int,int,int,str)
        the year, month, day and 'AM'|'PM' at start time
    ted : tuple, (int,int,int,str)
        the year, month, day and 'AM'|'PM' at end time
    '''
    import FireObj, FireIO

    # loop over all days during the period
    endloop = False  # flag to control the ending of the loop
    t = list(tst)    # t is the time (year,month,day,ampm) for each step
    while endloop == False:
        # print(t)

        # create and save file as of t
        save_sum_1d(tst,t)

        # time flow control
        #  - if t reaches ted, set endloop to True to stop the loop
        if FireObj.t_dif(t,ted)==0:
            endloop = True

        #  - update t with the next time stamp
        t = FireObj.t_nb(t,nb='next')


def style_summary(df_summary):
    ''' Format the summary dataframe

    Parameters
    ----------
    df_summary : pandas dataframe
        summary data

    Returns
    -------
    df_summary_fmt : pd dataframe
        formated summary data
    '''
    sum_fm    = {'Number of fires':'{:0.0f}',
                 'Number of fire pixels':'{:0.0f}',
                 'Fire total area':'{:0.1f}',
                 'Fire total perimeter':'{:0.2f}',
                 'Fire mean duration':'{:0.2f}',
                 'Number of active fires':'{:0.0f}',
                 'Number of new fire pixels':'{:0.0f}',
                 'Active fire total area':'{:0.1f}',
                 'Active fire total perimeter':'{:0.2f}',
                 'Active fire line total length':'{:0.3f}',
                 'Active fire mean FRP':'{:0.2f}',
                 'Active fire mean duration':'{:0.2f}',
                 }

    df_summary = df_summary[sum_fm.keys()]
    df_summary_fmt = df_summary.style.format(sum_fm)

    return df_summary_fmt

def pct_summary(df_summary):
    ''' Calculate the percentages of some fire attributes from different fire types

    Parameters
    ----------
    df_summary : pandas dataframe
        summary data

    Returns
    -------
    df_pct : pd dataframe
        percentage summary data
    '''
    from FireConsts import FTYP

    # extract some fire attributes only
    fattrs = ['Number of fires','Number of active fires','Number of fire pixels',
              'Number of new fire pixels','Fire total area','Fire total perimeter',
              'Active fire total area','Active fire total perimeter','Active fire line total length']
    df_pct = df_summary[fattrs].T

    # calculate the percentage (and remove original columns)
    for i in range(len(FTYP)):
        df_pct['Pct_'+FTYP[i]] = df_pct[FTYP[i]]/df_pct['All']*100
        df_pct = df_pct.drop(columns=[FTYP[i]])
        df_pct = df_pct.rename(columns={('Pct_'+FTYP[i]):FTYP[i]})
    df_pct = df_pct.drop(columns='All')

    # style the dataframe
    df_pct = df_pct.T.style.format('{:.1f}')

    return df_pct

def add_heritage(ted):
    ''' add heritage info up to ted to the ds_summary

    Parameters
    ----------
    ded : tuple (year, month, day, ampm)
        end date

    Returns
    -------
    df : pd DataFrame
        the df with heritage info
    '''
    import FireIO
    import pandas as pd

    # read the allfires at ted
    allfires = FireIO.load_fobj(ted)

    # derive the dictionary storing all heritage information
    dic_heri = {}
    for f in allfires.heritages:
        if f[1] in dic_heri.keys():
            dic_heri[f[1]] = list(dic_heri[f[1]] + [f[0]])
        else:
            dic_heri[f[1]] = [f[0]]

    # convert dict to pd DataFrame
    v = [','.join([str(elem) for elem in x]) for x in list(dic_heri.values())]
    df = pd.DataFrame(v,index=list(dic_heri.keys()),columns=['Fid_merge'])
    df.index.rename('Fid',inplace=True)

    # save
    FireIO.save_summarycsv(df,ted[0],op='heritage')

    # return df

def add_largefirelist(ted):
    ''' add large fire info up to ted to the ds_summary

    Parameters
    ----------
    ded : tuple (year, month, day, ampm)
        end date

    Returns
    -------
    df : pd DataFrame
        the df with large fire info
    '''
    import FireIO
    from FireConsts import falim
    import pandas as pd

    # read the allfires,gdf at ted
    allfires = FireIO.load_fobj(ted)
    gdf = FireIO.load_gdfobj(ted)

    # get all large fire ids (up to ted)
    fids_lf = list(gdf[(gdf.farea>falim) & (gdf.invalid==0)].index)

    # record the dictionary containing large fire ids and t_ed
    dic_lf = {}
    for fid in fids_lf:
        dic_lf[fid] =  allfires.fires[fid].t_ed

    # return dic_lf
    v = [str(x) for x in list(dic_lf.values())]
    df_lf = pd.DataFrame(v,index=list(dic_lf.keys()),columns=['t_ed'])
    df_lf.index.rename('Fid',inplace=True)

    # save
    FireIO.save_summarycsv(df_lf,ted[0],op='large')

    # return df_lf

if __name__ == "__main__":
    ''' The main code to record time series of summary data in a netcdf file
    The summary data are all data starting from (year,1,1,'AM')
    '''

    from datetime import date
    import time

    t1 = time.time()

    # set the end time
    tst=(2019,1,1,'AM')
    ted=(2019,12,31,'PM')

    # create summary file for one time step
    save_sum_1d(tst,ted)

    # for each day during the period, create summary file
    # save_sum_trng(tst,ted)

    # save fire lists of heritage and large fires
    add_heritage(ted)
    add_largefirelist(ted)

    t2 = time.time()
    print(f'{(t2-t1)/60.} minutes used to run code')

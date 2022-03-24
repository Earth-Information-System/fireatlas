""" FireIO
This module include functions used to read and save data
"""
# ------------------------------------------------------------------------------
#%% Read active fire data
# ------------------------------------------------------------------------------



def save_af(data,d):
    ''' Save a daily allfires object to a pickle file

    Parameters
    ----------
    data : pd with all fires of a month
    t : tuple, (year,month,day,str)
        the day and 'AM'|'PM' during the intialization
    '''

    import pickle
    import os
    from FireConsts import dirpjdata

    # get output file name
    dir_temp = os.path.join(dirpjdata,'temp') + '/' + str(d.year) + '/'
    fnm = os.path.join(dir_temp,'VNP14IMGML.'+d.strftime('%Y%m')+'.C1.05.pkl')

    # check folder
    check_filefolder(fnm)

    # save
    with open(fnm,'wb') as f:
        pickle.dump(data, f)

def read_VNP14ML04_clip(t,region,ext=False):
    ''' Read monthly S-NPP VIIRS fire location (C1.04)

    Parameters
    ----------
    t : tuple, (int,int,int,str)
        the year, month, day and 'AM'|'PM' during the intialization
    region: str,
        test region, currently possible values are:
            'AK': Alaska,
            'AK2015': small test area in Alaska 2015,
            'SL2014': Slave Lake 2014,
            'NY2020': Northern Yakutia 2020, 
            'YAK2020': Yakutia 2020
        (shapefiles not implemented yet)

    Returns
    -------
    vlist : list
        (lat,lon) tuple of all daily active fires
    '''
    from FireConsts import dirextdata, dirpjdata, ext_all

    from datetime import date
    import os
    import pandas as pd
    from shapely.geometry import Point
    import geopandas as gpd
    import pickle
    
    # extent form shapefile
    # shp_name = 'AK.shp'
    # shp = get_any_shp(shp_name) # shapefile of test region
    # ext = list(shp.bounds)
    # if shp_name == 'AK.shp': ext[2] = 140*-1  # Alaska extends over Antimeridian, therefore bounds are screwed up
    
    
    d = date(*t[:-1])
    
    # for the first day of the month we read in the whole file, filter it
    # according to region and data quality and save it to pickle
    
    #check if pickle for that month already exists, then use it
    dir_temp = os.path.join(dirpjdata,'temp') + '/' + str(d.year) + '/'
    fnm_pckl = os.path.join(dir_temp,'VNP14IMGML.'+d.strftime('%Y%m')+'.C1.05.pkl')
    
    if os.path.exists(fnm_pckl):
        with open(fnm_pckl,'rb') as f:
            gdf = pickle.load(f)
    else:
        # set monthly file name
        dirFC = os.path.join(dirextdata,'VNP14IMGML05') + '/'
        fnmFC = os.path.join(dirFC,'VNP14IMGML.'+d.strftime('%Y%m')+'.C1.05.txt')

        # read and extract
        usecols = ['YYYYMMDD','HHMM','Lat','Lon','Line','Sample','FRP','Confidence','Type','DNFlag']
        if os.path.exists(fnmFC):
            # read dataframe
            df = pd.read_csv(fnmFC,parse_dates=['YYYYMMDD'],usecols=usecols,skipinitialspace=True)
    
            # in Collection 04, some FRP values are '*******', causing problems
            if df.dtypes.FRP.name == 'object':
                df.FRP = df.FRP.replace('*******',0).astype('float')
            
            # set extent from region
            if region in ext_all:
                ext = ext_all[region]
            else: # use circumpolar boreal-arctic (>50N)
                ext = [-168, 50, 180, 80]

            newfirepixels = df.loc[(df['Lat'] > ext[1]) & (df['Lat'] < ext[3]) &
                                  (df['Lon'] > ext[0]) & (df['Lon'] < ext[2]) &
                                  ((df['Type'] == 0) | (df['Type'] == 3)) &
                                  ((df['Confidence'] == 'nominal') | (df['Confidence'] == 'high'))]   # use type==0 (vf) nd 3 (offshore)
            
            # there is spurious data on June 24, 2020, 11:06
            if t[0] == 2020:
                newfirepixels = newfirepixels.loc[~((newfirepixels['YYYYMMDD'] == pd.Timestamp(date(2020,6,24))) & 
                                                    (newfirepixels['HHMM'] == 1106))]
            # spatial filter
            point_data = [Point(xy) for xy in zip(newfirepixels['Lon'], newfirepixels['Lat'])]
            gdf = gpd.GeoDataFrame(newfirepixels, geometry=point_data)
            
                
            
            # save to pickle
            save_af(gdf,d)
        else:
            return None
            
    ## for all options: apply temporal filters
    # temporal (daily) filter
    gdf = gdf.loc[(gdf['YYYYMMDD']==d.strftime('%Y-%m-%d'))]

    # overpass time filter (AM or PM)
    # DNFlag - 0: night; 1: day
    # vlh = (df['HHMM']*0.01+df['lon']/15.)  # local time
    # tpm = (vlh > 7) & (vlh < 19)           # if local time in [7,19], set as PM overpass
    tpm = (gdf['DNFlag'] == 'night')
    if t[-1] == 'AM':
        gdf = gdf.loc[~tpm]
    elif t[-1] == 'PM':
        gdf = gdf.loc[tpm]
        
    # convert to list of (lat,lon)
    vlist = gdf[['Lat','Lon','Line','Sample','FRP']].values.tolist()
        
    return vlist

def read_AFPviirs(t,region):
    ''' Read and extract daily active fire pixels from monthly or daily VIIRS fire location data
    Parameters
    ----------
    t : tuple, (year,month,day,str)
        the day and 'AM'|'PM' during the intialization
        
    Returns
    -------
    vlist : list
        (lat,lon,FRP) tuple of all daily active fires
    '''
    
    # read monthly data if it exist
    # vlist = read_VNP14ML_CA(t)
    # vlist = read_VNP14ML04_CA(t)
    vlist = read_VNP14ML04_clip(t,region)
    
    # if monthly data are unavailable, use NRT daily data
    # if vlist is None: vlist = read_VNP14TDL_CA(t)
    
    return vlist

def read_AFP(t,src='viirs', region=None):
    ''' Read and extract daily active fire pixels from csv file or fire locaiton file
    
    Parameters
    ----------
    t : tuple, (int,int,int,str)
        the year, month, day and 'AM'|'PM' during the intialization
    src : str
        'viirs' (location data) | 'csv' (pre-filtered csv data)
        
    Returns
    -------
    vlist : list
        (lat,lon,FRP) tuple of active fires during the time step
    '''
    
    if src == 'viirs':
        vlist = read_AFPviirs(t,region)
        return vlist
    else:
        print('Please set src to viirs')
        return None

def read_mcd64_pixels(year,ext=[]):
    '''Reads all Modis burned area pixels of a year from text file
    
    Parameters
    ----------
    year: int
    ext: list or tuple, (minx, miny, maxx, maxy) extent for filtering pixels
        
    Returns
    -------
    df: pd.dataframe, dataframe containing lat, lon and day of burning 
        from modis burned area '''
        
    from FireConsts import mcd64dir#, dirpjdata
    # import glob
    import numpy as np
    import pandas as pd
    
    # filelist_viirs = glob.glob(dirpjdata+str(year)+'/Snapshot/*_NFP.txt')
    # df = pd.concat([pd.read_csv(i, dtype = {'lat':np.float64, 'lon':np.float64}, 
    #                             usecols = ['lat', 'lon']) 
    #                 for i in filelist_viirs], ignore_index=True)
    # df['sensor'] = 'viirs'
    
    filelist_mcd64 = mcd64dir + 'ba_centroids_'+str(year)+'.csv'
    df = pd.read_csv(filelist_mcd64, dtype = {'lat':np.float64, 'lon':np.float64},
                     usecols = ['lat', 'lon','doy'])
    # filter by bounding box
    if len(ext) == 4:
        minx, miny, maxx, maxy = ext
        df = df.loc[(df['lat']>=miny) & (df['lat']<= maxy) & (df['lon']>=minx) & (df['lon']<=maxx)]
    
    # df = pd.concat([df, df2])
    
    return df


# ------------------------------------------------------------------------------
#%% Read other datasets
# ------------------------------------------------------------------------------
def load_mcd64(year,xoff=0,yoff=0,xsize=None,ysize=None):
    '''get annual circumpolar modis burned/unburned tif
    optional: clip using UL pixel and pixel dimensions
    
    Parameters
    ----------
    year: int
    xoff, yoff: int, UL pixel coordinates for clipping
    xsize,ysize: int, pixels in x and y direction for clipping
        
    Returns
    -------
    arr: np.array, clipped image'''
    
    from FireConsts import mcd64dir
    import gdal
    
    fnm = mcd64dir + 'mcd64_' + str(year) + '.tif'
    ds = gdal.Open(fnm)
    #arr = ds.ReadAsArray(xsize=xsize, ysize=ysize)
    arr = ds.ReadAsArray(xoff=xoff, yoff=yoff, xsize=xsize, ysize=ysize)
    
    return arr

def get_any_shp(filename):
    ''' get shapefile of any region
    '''
    from FireConsts import dirextdata
    import geopandas as gpd
    import os
    
    # find the california shapefile
    dirshape = os.path.join(dirextdata,'shapefiles')
    statefnm = os.path.join(dirshape,filename)
    
    # read the geometry
    shp = gpd.read_file(statefnm).iloc[0].geometry
    
    return shp

def get_LCT(locs):
    ''' Get land cover type for active fires
    
    Parameters
    ----------
    locs : list of lists (nx2)
        lat and lon values for each active fire detection
        
    Returns
    -------
    vLCT : list of ints
        land cover types for all input active fires
    '''
    from FireConsts import dirextdata
    
    from osgeo import gdal
    import pyproj
    import os
    
    # read NLCD 500m data
    fnmLCT = os.path.join(dirextdata,'nlcd_510m.tif')
    dataset = gdal.Open(fnmLCT, gdal.GA_ReadOnly)
    band = dataset.GetRasterBand(1)
    data = band.ReadAsArray()
    gt = dataset.GetGeoTransform()
    srs = dataset.GetSpatialRef()
    pfunc = pyproj.Proj(srs.ExportToProj4())
    
    # extract LCTs for each points
    vLCT = []
    for lat,lon in locs:
        x,y = pfunc(lon,lat)
        i,j = int((x-gt[0])/gt[1]),int((y-gt[3])/gt[5])
        vLCT.append(data[j,i])
        
    return vLCT

def get_FM1000(t,loc):
    ''' Get fm1000 for a point at t
    
    Parameters
    ----------
    t : datetime date
        date
    loc : tuple
        (lat,lon) value
        
    Returns
    -------
    FM1000_loc : list of floats
        fm1000 value for all input active fires
    '''
    from FireConsts import dirextdata
    
    import xarray as xr
    import os
    
    import warnings
    warnings.simplefilter("ignore")
    
    # read annual fm1000 data
    dirGridMET = os.path.join(dirextdata,'GridMET') + '/'
    fnm = dirGridMET + 'fm1000_'+t.strftime('%Y')+'.nc'
    ds = xr.open_dataset(fnm)
    FM1000_all = ds['dead_fuel_moisture_1000hr']
    
    # extract daily data at t
    try:
        FM1000_day = FM1000_all.sel(day=t.strftime('%Y-%m-%d'))
    except:  # if data are not available, use the last available date
        FM1000_day = FM1000_all.isel(day=-1)
        
    # extract data near the given location
    FM1000_loc = FM1000_day.sel(lon=loc[1],lat=loc[0],method='nearest').item()
    
    return FM1000_loc

def get_stFM1000(fhull,locs,t):
    ''' get FM1000 for a fire at a time
    
    Parameters
    ----------
    fhull : geometry
        the hull of the fire
    locs : list of lists (nx2)
        lat and lon values for each active fire detection
    t : tuple, (int,int,int,str)
        the year, month, day and 'AM'|'PM'
    
    Returns
    -------
    FM1000 : list of floats
        fm1000 value for all input active fires
    '''
    import FireClustering
    
    from datetime import date
    
    # centroid (lat,lon) is defined as that of hull (faster) or locs
    if fhull is not None:
        cent = (fhull.centroid.y,fhull.centroid.x)
    else:
        cent = FireClustering.cal_centroid(locs)
        
    # datetime date of the time tuple
    d_st = date(*t[:-1])
    
    # call get_FM1000 function to extract fm1000 for the centroid at a time
    FM1000 = get_FM1000(d_st, cent)
    
    return FM1000

# ------------------------------------------------------------------------------
#%% read and load object, gdf and summary related files
# ------------------------------------------------------------------------------
def check_filefolder(fnm):
    ''' if the folder containing a file does not exist, make it

    Parameters
    ----------
    fnm : str
        file name
    '''
    import os

    # folder name
    dirnm = os.path.dirname(fnm)

    # create the folder if needed
    if not os.path.exists(dirnm):
        os.makedirs(dirnm)

def check_folder(dirnm):
    ''' if the folder does not exist, make it

    Parameters
    ----------
    dirfnm : str
        folder name
    '''
    import os
    
    # create the folder if needed
    if not os.path.exists(dirnm):
        os.makedirs(dirnm)

def get_path(year,region):
    ''' returns the base output directory for the year and region '''
    from FireConsts import dirpjdata
    if region:
        path = dirpjdata+region+'/'+year
    else:
        path = dirpjdata+year
    return path

def correct_dtype(gdf,op=''):
    ''' correct the datatype for gdfs loaded from geojson files

    Parameters
    ----------
    gdf : gpd DataFrame
        the gdf directly read from the geojson file
    op : str
        the option for the geojson data -
        '': fire perimeter and attributes
        'FL': active fire line
        'NFP': new fire pixels

    Returns
    -------
    gdf : gpd DataFrame
        the gdf with correct datatype as given in FireConsts
    '''
    from FireConsts import dd

    # explicitly set the attributes data types
    if op == '':
        for v,tp in dd.items():
            gdf[v] = gdf[v].astype(tp)
    else:
        gdf['fid'] = gdf['fid'].astype('int')

    return gdf

def get_fobj_fnm(t,region=''):
    ''' Return the fire object pickle file name at a time step
    Parameters
    ----------
    t : tuple, (year,month,day,str)
        the day and 'AM'|'PM' during the intialization
    
    Returns
    ----------
    fnm : str
        pickle file name
    '''
    from datetime import date
    
    d = date(*t[:-1])
    path = get_path(d.strftime('%Y'),region)
    fnm = path+'/Serialization/'+d.strftime('%Y%m%d')+t[-1]+'.pkl'
    return fnm

def check_fobj(t,region=''):
    ''' Check if the pickle file storing a daily allfires object exists

    Parameters
    ----------
    t : tuple, (year,month,day,str)
        the day and 'AM'|'PM' during the intialization
    '''
    import os

    fnm = get_fobj_fnm(t,region)
    return os.path.exists(fnm)

def save_fobj(data,t,region=''):
    ''' Save a daily allfires object to a pickle file

    Parameters
    ----------
    data : obj of Allfires class
        daily allfires
    t : tuple, (year,month,day,str)
        the day and 'AM'|'PM' during the intialization
    '''

    import pickle
    # import os

    # get output file name
    fnm = get_fobj_fnm(t,region)

    # check folder
    check_filefolder(fnm)

    # save
    with open(fnm,'wb') as f:
        pickle.dump(data, f)

def load_fobj(t,region=''):
    ''' Load a daily allfires object from a pickle file

    Parameters
    ----------
    t : tuple, (year,month,day,str)
        the day and 'AM'|'PM' during the intialization

    Returns
    ----------
    data : obj of Allfires class
        daily allfires
    '''
    import pickle

    # get file name
    fnm = get_fobj_fnm(t,region)

    # load data
    with open(fnm,'rb') as f:
        data = pickle.load(f)
    return data

def get_gdfobj_fnm(t,op='',region=''):
    ''' Return the fire object pickle file name at a time step
    Parameters
    ----------
    t : tuple, (year,month,day,str)
        the day and 'AM'|'PM' during the intialization
    op : str
        the option for the geojson data -
        '': fire perimeter and attributes
        'FL': active fire line
        'NFP': new fire pixels
    Returns
    ----------
    fnm : str
        gdf file name
    '''
    from datetime import date

    d = date(*t[:-1])
    path = get_path(d.strftime('%Y'),region)
    if op == '':
        fnm = path+'/Snapshot/'+d.strftime('%Y%m%d')+t[-1]+'.gpkg'
    else:
        fnm = path+'/Snapshot/'+d.strftime('%Y%m%d')+t[-1]+'_'+op+'.gpkg'

    return fnm

def check_gdfobj(t,op='',region=''):
    ''' Check if the gpkg file storing a daily allfires attributes exists

    Parameters
    ----------
    t : tuple, (year,month,day,str)
        the day and 'AM'|'PM' during the intialization
    op : str
        the option for the geojson data -
        '': fire perimeter and attributes
        'FL': active fire line
        'NFP': new fire pixels
    '''
    import os
    
    fnm = get_gdfobj_fnm(t,op=op,region=region)

    return os.path.exists(fnm)

def save_gdfobj(gdf,t,param='',fid='',op='',region=''):
    ''' Save geopandas to a gpgk file

    Parameters
    ----------
    gdf : geopandas DataFrame
        daily allfires diagnostic parameters
    t : tuple, (int,int,int,str)
        the year, month, day and 'AM'|'PM'
    param : string
        if empty: save daily allfires diagnostic dataframe
        if 'large': save daily single fire 
        if 'ign': save ignition layer
        if 'final': save final perimeters layer
    fid : int,
        fire id (needed only with param = 'large')
    '''

    # get file name
    if param == '':
        fnm = get_gdfobj_fnm(t,op=op,region=region)
    elif param == 'large':
        fnm = get_gdfobj_sf_fnm(t,fid,op=op,region=region)
    else:
        from datetime import date
        # get path
        d = date(*t[:-1])
        path = get_path(d.strftime('%Y'),region)
        # get file name
        fnm = path+'/Summary/'+param+d.strftime('%Y')+'.gpkg'
        
    # check folder
    check_filefolder(fnm)
    
    # reset the index column so a new fid column is created
    gdf = gdf.reset_index()
    if op == 'FL':
        gdf['fireid'] = gdf['fireid'].astype(int) # data types are screwed up in fline

    # save file
    gdf.to_file(fnm, driver='GPKG')

def load_gdfobj(t='',op='',region=''):
    ''' Load daily allfires diagnostic dataframe as geopandas gdf

    Parameters
    ----------
    t : tuple, (int,int,int,str)
        the year, month, day and 'AM'|'PM'
    Returns
    ----------
    gdf : geopandas DataFrame
        daily allfires diagnostic parameters
    '''
    import geopandas as gpd

    # get file name
    fnm = get_gdfobj_fnm(t,op=op,region=region)

    # read data as gpd DataFrame
    gdf = gpd.read_file(fnm)

    # correct the datatype
    gdf = correct_dtype(gdf,op=op)

    # set fireid as the index
    gdf = gdf.set_index('fireid')

    return gdf

def save_FP_txt(df,t,region=''):
    
    # get filename of new fire pixels product
    fnm = get_gdfobj_fnm(t,op='NFP',region=region)
    fnm = fnm[:-4]+'txt' # change ending to txt
    
    # check folder
    check_filefolder(fnm)
    
    # write out
    if len(df) > 0:
        df.to_csv(fnm)

def load_lake_geoms(t,fid,region=''):
    ''' Load final perimeters as geopandas gdf

    Parameters
    ----------
    t: time tuple,
        needed for the year
    fid : int, 
        the fire fid
    Returns
    ----------
    geoms : tuple (geom, geom)
        geometry of all lake perimeters within the fire
    '''
    import geopandas as gpd
    from datetime import date
    
    d = date(*t[:-1])
    path = get_path(d.strftime('%Y'),region)
    
    # get file name
    fnm_lakes = path+'/Summary/lakes'+d.strftime('%Y')+'.gpkg'
    
    # read data and extract target geometry
    gdf = gpd.read_file(fnm_lakes)
    gdf = gdf.set_index('mergid')
    if fid in gdf.index:
        gdf = gdf.loc[fid]
        geom_lakes = gdf['geometry']
    else:
        geom_lakes = None
        
    return geom_lakes


def get_gdfobj_sf_fnm(t,fid,op='',region=''):
    ''' Return the single fire fire object pickle file name at a time step
    Parameters
    ----------
    t : tuple, (year,month,day,str)
        the day and 'AM'|'PM' during the intialization
    fid : int
        fire id
    Returns
    ----------
    fnm : str
        gdf file name
    '''
    from FireConsts import dirpjdata
    from datetime import date
    d = date(*t[:-1])
    path = get_path(d.strftime('%Y'),region)
    if op == '':
        fnm = path+'/Largefire/F'+str(int(fid))+'_'+d.strftime('%Y%m%d')+t[-1]+'.gpkg'
    else:
        fnm = path+'/Largefire/F'+str(int(fid))+'_'+d.strftime('%Y%m%d')+t[-1]+'_'+op+'.gpkg'

    return fnm

def get_gdfobj_sf_fnms_year(year,fid,op='',region=''):
    ''' Return the single fire fire object pickle file name at a time step
    Parameters
    ----------
    year : int
        year
    fid : int
        fire id
    Returns
    ----------
    fnms : list of str
        geojson file names
    '''
    from FireConsts import dirpjdata
    from datetime import date
    from glob import glob
    
    path = get_path(year,region)
    if op == '':
        fnms = glob(path+'/Largefire/F'+str(int(fid))+'_??????????.gpkg')
    else:
        fnms = glob(path+'/Largefire/F'+str(int(fid))+'_??????????_'+op+'.gpkg')

    return fnms


def load_gdfobj_sf(t,fid,op='',region=''):
    ''' Load single fire daily allfires diagnostic dataframe to a geojson file

    Parameters
    ----------
    t : tuple, (int,int,int,str)
        the year, month, day and 'AM'|'PM'
    fid : int
        fire id
    Returns
    ----------
    gdf : geopandas DataFrame
        time series of daily single fire diagnostic parameters
    '''
    import geopandas as gpd
    import pandas as pd

    # get file name
    fnm = get_gdfobj_sf_fnm(t,fid,op=op,region=region)

    # read data as gpd DataFrame
    gdf = gpd.read_file(fnm)

    # parse time step to Datetime and set as index
    gdf['index'] = pd.to_datetime(gdf['index'])
    gdf = gdf.set_index('index')

    return gdf

def get_summary_fnm(t,region=''):
    ''' Return the fire summary file name at year end
    Parameters
    ----------
    t : tuple, (year,month,day,str)
        the day and 'AM'|'PM' during the intialization
    Returns
    ----------
    fnm : str
        summary netcdf file name
    '''
    from datetime import date
    import os

    d = date(*t[:-1])
    path = get_path(t[0],region)
    fnm = os.path.join(path,'Summary','fsummary_'+d.strftime('%Y%m%d')+t[-1]+'.nc')

    # check folder
    check_filefolder(fnm)

    return fnm

def get_summary_fnm_lt(t,region=''):
    ''' Return the latest time step before current time when summary file exists
    Parameters
    ----------
    t : tuple, (year,month,day,str)
        the day and 'AM'|'PM' during the intialization

    Returns
    ----------
    pt : tuple, (year, month, day, pmpm)
        the lastest time step with summary file
    '''
    import os
    from glob import glob
    import FireObj
    
    path = get_path(t[0],region)
    
    # if there's no summary file for this year, return the first time step of the year
    fnms = glob(os.path.join(path,'Summary','fsummary_*.nc'))
    if len(fnms) == 0:
        return None

    # if there is, find the nearest one
    endloop = False
    pt = FireObj.t_nb(t,nb='previous')
    while endloop == False:
        if os.path.exists(get_summary_fnm(pt,region=region)):
            return pt
        else:
            pt = FireObj.t_nb(pt,nb='previous')
            if pt[0] != t[0]:
                return None

def check_summary(t,region=''):
    ''' Check if the summary file storing a daily allfires object exists

    Parameters
    ----------
    t : tuple, (year,month,day,str)
        the day and 'AM'|'PM' during the intialization
    '''
    import os

    # get file name
    fnm = get_summary_fnm(t,region)

    # return if it's present
    return os.path.exists(fnm)

def save_summary(ds,t,region=''):
    ''' Save summary info as of t a netcdf file

    Parameters
    ----------
    ds : xarray dataset
        year end summary dataset
    t : tuple, (int,int,int,str)
        the year, month, day and 'AM'|'PM'
    '''
    # get file name
    fnm = get_summary_fnm(t,region)

    # check folder
    check_filefolder(fnm)

    # save netcdf file
    ds.to_netcdf(fnm)

def load_summary(t,region=''):
    ''' Load summary info from a netcdf file at t

    Parameters
    ----------
    t : tuple, (int,int,int,str)
        the year, month, day and 'AM'|'PM'

    Returns
    -------
    ds : xarray dataset
        summary dataset
    '''
    import xarray as xr

    # get file name
    fnm = get_summary_fnm(t,region)

    # read data as xarray Dataset
    ds = xr.open_dataset(fnm)

    return ds

def save_summarycsv(df,year,op='heritage',region=''):
    ''' save summary csv files

    Parameters
    ----------
    df : pandas DataFrame
        the data
    year : int
        year
    op : str
        option, 'heritage'|'large'
    '''
    import os
    
    path = get_path(year,region)
    
    fnm = os.path.join(path,'Summary','Flist_'+op+'_'+str(year)+'.csv')
    check_filefolder(fnm)

    df.to_csv(fnm)

def read_summarycsv(year,op='heritage',region=''):
    ''' read summary csv files

    Parameters
    ----------
    year : int
        year
    op : str
        option, 'heritage'|'large'

    Returns
    -------
    df : pandas DataFrame
        the data
    '''
    import pandas as pd
    import os
    
    path = get_path(year,region)
    
    fnm = os.path.join(path,'Summary','Flist_'+op+'_'+str(year)+'.csv')
    df = pd.read_csv(fnm,index_col=0)
    return df

def get_lts_VNP14IMGTDL(year=None):
    from FireConsts import dirextdata
    from datetime import date
    from glob import glob
    import os

    if year == None:
        year = date.today().year

    dirFC = os.path.join(dirextdata,'VNP14IMGTDL') + '/'
    fnms = glob(os.path.join(dirFC,'SUOMI_VIIRS_C2_Global_VNP14IMGTDL_NRT_'+str(year)+'*.txt'))
    fnms.sort()
    DOY_lts = int(os.path.splitext(os.path.basename(fnms[-1]))[0][-3:])

    return DOY_lts

def get_lts_serialization(year=None,region=''):
    ''' get the time with lastest pkl data
    '''
    from FireConsts import dirpjdata
    from datetime import date
    from glob import glob
    import os

    if year == None:
        year = date.today().year
    
    path = dirpjdata+d.strftime('%Y')
    if region:
        path = path+'_'+region

    fnms = glob(path+'/Serialization/*.pkl')

    if len(fnms) > 0:
        fnms.sort()
        fnm_lts = os.path.basename(fnms[-1])

        lts = [int(fnm_lts[0:4]),int(fnm_lts[4:6]),int(fnm_lts[6:8]),fnm_lts[8:10]]
    else:
        lts = None

    return lts

#%% other functions related to read/write

def save2gtif(arr, outfile, cols, rows, geotrans, proj):
    '''write out a geotiff'''
    
    import gdal
    
    driver = gdal.GetDriverByName('GTiff')
    ds = driver.Create(outfile, cols, rows, 1, gdal.GDT_Byte)
    ds.SetGeoTransform(geotrans)
    ds.SetProjection(proj)
    band = ds.GetRasterBand(1)
    band.WriteArray(arr)

def geo_to_polar(lon_arr, lat_arr):
    '''transform lists of geographic lat lon coordinates to polar LAEA grid (projected)'''
    
    import numpy as np
    import pyproj
    
    proj4str = ("epsg:3571")
    p_modis_grid = pyproj.Proj(proj4str)
        
    x_arr, y_arr = p_modis_grid(lon_arr, lat_arr)
    x_arr = np.array(x_arr)
    y_arr = np.array(y_arr)
    
    return x_arr, y_arr

def polar_to_geo(x_arr, y_arr):
    '''transform lists of geographic lat lon coordinates to polar LAEA grid (projected)'''
    
    import numpy as np
    import pyproj
    
    proj4str = ("epsg:3571")
    p_modis_grid = pyproj.Proj(proj4str)
        
    lon_arr, lat_arr = p_modis_grid(x_arr, y_arr,inverse=True)
    lon_arr = np.array(lon_arr)
    lat_arr = np.array(lat_arr)
    
    return lon_arr, lat_arr

def world2Pixel(gt, Xgeo, Ygeo):
    ''' Uses a geomatrix (gdal.GetGeoTransform()) to calculate the pixel 
    location of a geospatial coordinate'''
    
    import numpy as np
    
    gt = list(gt)
    Xpx = np.rint((Xgeo - gt[0]) / gt[1]).astype(int)
    Ypx = np.rint((Ygeo - gt[3]) / gt[5]).astype(int)
    
    return (Xpx, Ypx)

def pixel2World(gt, Xpixel, Ypixel):
    '''Uses a gdal geomatrix (gdal.GetGeoTransform()) to calculate the
    geospatial coordinate of a pixel location'''
    
    Xgeo = gt[0] + Xpixel*gt[1] + Ypixel*gt[2]
    Ygeo = gt[3] + Xpixel*gt[4] + Ypixel*gt[5]
    
    return (Xgeo, Ygeo)
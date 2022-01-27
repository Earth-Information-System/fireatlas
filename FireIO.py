""" FireIO
This module include functions used to read and save data
"""
# ------------------------------------------------------------------------------
# Read active fire data
# ------------------------------------------------------------------------------


def read_VNP14ML04_clip(t, ext = False):
    ''' Read monthly S-NPP VIIRS fire location (C1.04)

    Parameters
    ----------
    t : tuple, (int,int,int,str)
        the year, month, day and 'AM'|'PM' during the intialization

    Returns
    -------
    vlist : list
        (lat,lon) tuple of all daily active fires
    '''
    from FireConsts import dirextdata, dirpjdata

    from datetime import date
    import os
    import pandas as pd
    from shapely.geometry import Point
    import geopandas as gpd
    import pickle


    d = date(*t[:-1])
    
    # for the first day of the month we read in the whole file, filter it
    # according to region and data quality and save it to pickle
    
    #check if pickle for that month already exists, then use it
    dir_temp = os.path.join(dirpjdata,'temp') + '/'
    fnm_pckl = os.path.join(dir_temp,'VNP14IMGML.'+d.strftime('%Y%m')+'.C1.05.pkl')
    
    if os.path.exists(fnm_pckl):
        with open(fnm_pckl,'rb') as f:
            gdf = pickle.load(f)
    else:
        # set monthly file name
        dirFC = os.path.join(dirextdata,'VNP14IMGML05') + '/'
        fnmFC = os.path.join(dirFC,'VNP14IMGML.'+d.strftime('%Y%m')+'.C1.05.txt')

        # read and extract
        usecols = ['YYYYMMDD','HHMM','Lat','Lon','Sample','FRP','Confidence','Type','DNFlag']
        if os.path.exists(fnmFC):
            # read dataframe
            df = pd.read_csv(fnmFC,parse_dates=['YYYYMMDD'],usecols=usecols,skipinitialspace=True)
    
            # in Collection 04, some FRP values are '*******', causing problems
            if df.dtypes.FRP.name == 'object':
                df.FRP = df.FRP.replace('*******',0).astype('float')

            # preliminary spatial filter from shapefile extent and quality filter
            # ext = [-156, 65.5, -152, 66.3] # test area for AK 2015
            #ext = [-119, 60, -110, 64] # test area for slave lake 2014
            ext = [106, 55, 163, 72] # test area for yakutia 2020
            #ext = [147, 65, 154, 70] # test area for large fire in yakutia only
            # shp_name = 'AK.shp'
            # shp = get_any_shp(shp_name) # shapefile of test region
            # ext = list(shp.bounds)
            # if shp_name == 'AK.shp': ext[2] = 140*-1  # Alaska extends over Antimeridian, therefore bounds are screwed up
            newfirepixels = df.loc[(df['Lat'] > ext[1]) & (df['Lat'] < ext[3]) &
                                  (df['Lon'] > ext[0]) & (df['Lon'] < ext[2]) &
                                  ((df['Type'] == 0) | (df['Type'] == 3)) &
                                  ((df['Confidence'] == 'nominal') | (df['Confidence'] == 'high'))]   # use type==0 (vf) nd 3 (offshore)
    
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
    vlist = gdf[['Lat','Lon','FRP']].values.tolist()
        
    return vlist
        

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
    dir_temp = os.path.join(dirpjdata,'temp') + '/'
    fnm = os.path.join(dir_temp,'VNP14IMGML.'+d.strftime('%Y%m')+'.C1.05.pkl')

    # check folder
    check_filefolder(fnm)

    # save
    with open(fnm,'wb') as f:
        pickle.dump(data, f)

def read_VNP14ML04_CA(t, calext=[-125, -114, 32, 42.5]):
    ''' Read monthly S-NPP VIIRS fire location (C1.04)

    Parameters
    ----------
    t : tuple, (int,int,int,str)
        the year, month, day and 'AM'|'PM' during the intialization

    Returns
    -------
    vlist : list
        (lat,lon) tuple of all daily active fires
    '''
    from FireConsts import dirextdata

    from datetime import date
    import os
    import pandas as pd
    from shapely.geometry import Point
    import geopandas as gpd


    d = date(*t[:-1])
    # set monthly file name
    dirFC = os.path.join(dirextdata,'VNP14IMGML04') + '/'
    fnmFC = os.path.join(dirFC,'VNP14IMGML.'+d.strftime('%Y%m')+'.C1.04.txt')

    # read and extract
    # usecols = ['YYYYMMDD','HHMM','lat','lon','sample','FRP','conf','type']
    usecols = ['YYYYMMDD','HHMM','Lat','Lon','Sample','FRP','Confidence','Type','DNFlag']
    if os.path.exists(fnmFC):
        # read dataframe
        df = pd.read_csv(fnmFC,parse_dates=['YYYYMMDD'],usecols=usecols,skipinitialspace=True)

        # in Collection 04, some FRP values are '*******', causing problems
        if df.dtypes.FRP.name == 'object':
            df.FRP = df.FRP.replace('*******',0).astype('float')

        # temporal (daily) filter
        df = df.loc[(df['YYYYMMDD']==d.strftime('%Y-%m-%d'))]

        # overpass time filter (AM or PM)
        # DNFlag - 0: night; 1: day
        # vlh = (df['HHMM']*0.01+df['lon']/15.)  # local time
        # tpm = (vlh > 7) & (vlh < 19)           # if local time in [7,19], set as PM overpass
        tpm = (df['DNFlag'] == 1)
        if t[-1] == 'AM':
            df = df.loc[~tpm]
        elif t[-1] == 'PM':
            df = df.loc[tpm]

        # preliminary spatial filter and quality filter
        newfirepixels = df.loc[(df['Lat'] > calext[2]) & (df['Lat'] < calext[3]) &
                              (df['Lon'] > calext[0]) & (df['Lon'] < calext[1]) &
                              (df['Type'] == 0)]   # use type==0 only

        # spatial filter
        shp_Cal = get_Cal_shp() # CA shapefile
        point_data = [Point(xy) for xy in zip(newfirepixels['Lon'], newfirepixels['Lat'])]
        vmon_gpd = gpd.GeoDataFrame(newfirepixels, geometry=point_data)
        vmon_gpd = vmon_gpd[vmon_gpd['geometry'].within(shp_Cal)]

        # convert to list of (lat,lon)
        vlist = vmon_gpd[['Lat','Lon','FRP']].values.tolist()

        return vlist
    else:
        return None

def read_VNP14TDL_CA(t, calext=[-125, -114, 32, 42.5]):
    ''' Read daily NRT S-NPP VIIRS fire location

    Parameters
    ----------
    t : tuple, (int,int,int,str)
        the year, month, day and 'AM'|'PM' during the intialization
    calext : list
        the box[lon0, lon1, lat0, lat1] of california boundary

    Returns
    -------
    vlist : list
        (lat,lon) tuple of all daily active fires
    '''
    from FireConsts import dirextdata

    import pandas as pd
    import os
    from shapely.geometry import Point
    import geopandas as gpd
    from datetime import date

    d = date(*t[:-1])

    # derive monthly file name with path
    dirFC = os.path.join(dirextdata,'VNP14IMGTDL') + '/'
    fnmFC = os.path.join(dirFC,'SUOMI_VIIRS_C2_Global_VNP14IMGTDL_NRT_'+d.strftime('%Y%j')+'.txt')

    # read and extract
    usecols = ['latitude','longitude','daynight','frp','confidence']
    if os.path.exists(fnmFC):
        # read data
        df = pd.read_csv(fnmFC,usecols=usecols,skipinitialspace=True)

        # overpass time filter (AM or PM), use day/night flag in the NRT data (N: AM; D: PM)
        if t[-1] == 'AM':
            tpm = df['daynight']=='N'
            df = df.loc[tpm]
        elif t[-1] == 'PM':
            tpm = df['daynight']=='D'
            df = df.loc[tpm]

        # preliminary spatial filter (and quality filter)
        newfirepixels = df.loc[(df['latitude'] > calext[2]) & (df['latitude'] < calext[3]) &
                              (df['longitude'] > calext[0]) & (df['longitude'] < calext[1]) &
                              ((df['confidence'] == 'nominal') | (df['confidence'] == 'high'))]

        # spatial filter
        shp_Cal = get_Cal_shp()
        point_data = [Point(xy) for xy in zip(newfirepixels['longitude'], newfirepixels['latitude'])]
        vmon_gpd = gpd.GeoDataFrame(newfirepixels, geometry=point_data)
        vmon_gpd = vmon_gpd[vmon_gpd['geometry'].within(shp_Cal)]

        # convert to list
        vlist = vmon_gpd[['latitude','longitude','frp']].values.tolist()
        return vlist
    else:
        return None

def read_AFPviirs(t):
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
    vlist = read_VNP14ML04_clip(t)
    
    # if monthly data are unavailable, use NRT daily data
    # if vlist is None: vlist = read_VNP14TDL_CA(t)
    
    return vlist

def read_AFP(t,src='viirs'):
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
        vlist = read_AFPviirs(t)
        return vlist
    else:
        print('Please set src to viirs')
        return None

def read_FRP_VNP14TDL_box(t,box):
    ''' Read daily NRT S-NPP VIIRS fire location and FRP for a box in California
    '''
    return read_VNP14TDL_CA(t, calext=[box[0],box[2],box[1],box[3]])

def read_FRP_VNP14ML04_box(t,box):
    ''' Read monthly NRT S-NPP VIIRS fire location and FRP for a box in California
    '''
    return read_VNP14ML04_CA(t, calext=[box[0],box[2],box[1],box[3]])

def read_FRPviirs_box(t,box):
    ''' Read and extract daily FRP from monthly or daily VIIRS fire location data
    
    Parameters
    ----------
    t : tuple, (year,month,day,str)
        the day and 'AM'|'PM' during the intialization
    box : list, [lonmin,latmin,lonmax,latmax]
        bounds used for data input
        
    Returns
    -------
    vlist : list
        (lat,lon,FRP) tuple of all daily FRP
    '''
    
    # read monthly data if it exist
    vlist = read_FRP_VNP14ML04_box(t,box)
    
    # if monthly data are unavailable, use NRT daily data
    if vlist is None: vlist = read_FRP_VNP14TDL_box(t,box)
    
    return vlist

# ------------------------------------------------------------------------------
# Read other datasets
# ------------------------------------------------------------------------------
def get_any_shp(filename):
    ''' get shapefile of California
    '''
    from FireConsts import dirextdata
    import geopandas as gpd
    import os
    
    # find the california shapefile
    dircalshape = os.path.join(dirextdata,'shapefiles')
    statefnm = os.path.join(dircalshape,filename)
    
    # read the geometry
    shp = gpd.read_file(statefnm).iloc[0].geometry
    
    return shp

def get_Cal_shp():
    ''' get shapefile of California
    '''
    from FireConsts import dirextdata
    import geopandas as gpd
    import os
    
    # find the california shapefile
    dircalshape = os.path.join(dirextdata,'Calshape')
    statefnm = os.path.join(dircalshape,'California.shp')
    
    # read the geometry
    shp_Cal = gpd.read_file(statefnm).iloc[0].geometry
    
    return shp_Cal

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
# read and load object, gdf and summary related files
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

def get_fobj_fnm(t):
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
    from FireConsts import dirpjdata
    from datetime import date
    d = date(*t[:-1])
    fnm = dirpjdata+d.strftime('%Y')+'/Serialization/'+d.strftime('%Y%m%d')+t[-1]+'.pkl'
    return fnm

def check_fobj(t):
    ''' Check if the pickle file storing a daily allfires object exists

    Parameters
    ----------
    t : tuple, (year,month,day,str)
        the day and 'AM'|'PM' during the intialization
    '''
    import os

    fnm = get_fobj_fnm(t)
    return os.path.exists(fnm)

def save_fobj(data,t):
    ''' Save a daily allfires object to a pickle file

    Parameters
    ----------
    data : obj of Allfires class
        daily allfires
    t : tuple, (year,month,day,str)
        the day and 'AM'|'PM' during the intialization
    '''

    import pickle
    import os

    # get output file name
    fnm = get_fobj_fnm(t)

    # check folder
    check_filefolder(fnm)

    # save
    with open(fnm,'wb') as f:
        pickle.dump(data, f)

def load_fobj(t):
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
    fnm = get_fobj_fnm(t)

    # load data
    with open(fnm,'rb') as f:
        data = pickle.load(f)
    return data

def get_gdfobj_fnm(t,op=''):
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
    from FireConsts import dirpjdata
    from datetime import date

    d = date(*t[:-1])
    if op == '':
        fnm = dirpjdata+d.strftime('%Y')+'/Snapshot/'+d.strftime('%Y%m%d')+t[-1]+'.geojson'
    else:
        fnm = dirpjdata+d.strftime('%Y')+'/Snapshot/'+d.strftime('%Y%m%d')+t[-1]+'_'+op+'.geojson'

    return fnm

def check_gdfobj(t,op=''):
    ''' Check if the geojson file storing a daily allfires attributes exists

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
    from FireConsts import dirpjdata
    from datetime import date
    import os

    d = date(*t[:-1])
    fnm = get_gdfobj_fnm(t)

    return os.path.exists(fnm)

def save_gdfobj(gdf,t,op=''):
    ''' Save daily allfires diagnostic dataframe to a geojson file

    Parameters
    ----------
    gdf : geopandas DataFrame
        daily allfires diagnostic parameters
    t : tuple, (int,int,int,str)
        the year, month, day and 'AM'|'PM'
    '''

    # get file name
    fnm = get_gdfobj_fnm(t,op=op)

    # check folder
    check_filefolder(fnm)

    # save file
    gdf.to_file(fnm, driver='GeoJSON')

def load_gdfobj(t,op=''):
    ''' Load daily allfires diagnostic dataframe to a geojson file

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
    fnm = get_gdfobj_fnm(t,op=op)

    # read data as gpd DataFrame
    gdf = gpd.read_file(fnm)

    # correct the datatype
    gdf = correct_dtype(gdf,op=op)

    # set fid as the index
    gdf = gdf.set_index('fid')

    return gdf


def get_gdfobj_sf_fnm(t,fid,op=''):
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

    if op == '':
        fnm = dirpjdata+d.strftime('%Y')+'/Largefire/F'+str(int(fid))+'_'+d.strftime('%Y%m%d')+t[-1]+'.geojson'
    else:
        fnm = dirpjdata+d.strftime('%Y')+'/Largefire/F'+str(int(fid))+'_'+d.strftime('%Y%m%d')+t[-1]+'_'+op+'.geojson'

    return fnm

def get_gdfobj_sf_fnms_year(year,fid,op=''):
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

    if op == '':
        fnms = glob(dirpjdata+str(year)+'/Largefire/F'+str(int(fid))+'_??????????.geojson')
    else:
        fnms = glob(dirpjdata+str(year)+'/Largefire/F'+str(int(fid))+'_??????????_'+op+'.geojson')

    return fnms

def save_gdfobj_sf(gdf,t,fid,op=''):
    ''' Save daily single fire allfires diagnostic dataframe to a geojson file

    Parameters
    ----------
    gdf : geopandas DataFrame
        daily allfires diagnostic parameters
    t : tuple, (int,int,int,str)
        the year, month, day and 'AM'|'PM'
    fid : int
        fire id
    '''
    # get file name
    fnm = get_gdfobj_sf_fnm(t,fid,op=op)

    # check folder
    check_filefolder(fnm)

    # save data to file
    gdf.to_file(fnm, driver='GeoJSON')


def load_gdfobj_sf(t,fid,op=''):
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
    fnm = get_gdfobj_sf_fnm(t, fid, op=op)

    # read data as gpd DataFrame
    gdf = gpd.read_file(fnm)

    # parse time step to Datetime and set as index
    gdf['index'] = pd.to_datetime(gdf['index'])
    gdf = gdf.set_index('index')

    return gdf

def load_gdfobj_sf_final(year,fid,op=''):
    ''' Load single fire data at the final day of fire

    Parameters
    ----------
    fid : int
        fire id
    Returns
    ----------
    gdf : geopandas DataFrame
        time series of daily single fire diagnostic parameters
    '''
    import geopandas as gpd
    import pandas as pd

    # get the final day of fire
    t = load_fobj((year,12,31,'PM')).fires[fid].t_ed

    # read data for small fire
    gdf_sf = load_gdfobj_sf(t, fid, op='')

    return gdf_sf

def save_gdfobj_ign(gdf,t):
    ''' Save daily single fire allfires diagnostic dataframe to a geojson file

    Parameters
    ----------
    gdf : geopandas DataFrame
        daily allfires diagnostic parameters
    t : tuple, (int,int,int,str)
        the year, month, day and 'AM'|'PM'
    fid : int
        fire id
    '''
    from datetime import date
    from FireConsts import dirpjdata
    d = date(*t[:-1])
    
    # get file name
    fnm = dirpjdata+d.strftime('%Y')+'/Summary/ignitions'+d.strftime('%Y')+'.geojson'

    # check folder
    check_filefolder(fnm)

    # save data to file
    gdf.to_file(fnm, driver='GeoJSON')



def get_summary_fnm(t):
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
    from FireConsts import dirpjdata
    from datetime import date
    import os

    d = date(*t[:-1])
    fnm = os.path.join(dirpjdata, d.strftime('%Y'),'Summary','fsummary_'+d.strftime('%Y%m%d')+t[-1]+'.nc')

    # check folder
    check_filefolder(fnm)

    return fnm

def get_summary_fnm_lt(t):
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
    from FireConsts import dirpjdata
    from datetime import date
    import os
    from glob import glob
    import FireObj

    # if there's no summary file for this year, return the first time step of the year
    fnms = glob(os.path.join(dirpjdata, str(t[0]),'Summary','fsummary_*.nc'))
    if len(fnms) == 0:
        return None

    # if there is, find the nearest one
    endloop = False
    pt = FireObj.t_nb(t,nb='previous')
    while endloop == False:
        if os.path.exists(get_summary_fnm(pt)):
            return pt
        else:
            pt = FireObj.t_nb(pt,nb='previous')
            if pt[0] != t[0]:
                return None

def check_summary(t):
    ''' Check if the summary file storing a daily allfires object exists

    Parameters
    ----------
    t : tuple, (year,month,day,str)
        the day and 'AM'|'PM' during the intialization
    '''
    import os

    # get file name
    fnm = get_summary_fnm(t)

    # return if it's present
    return os.path.exists(fnm)

def save_summary(ds,t):
    ''' Save summary info as of t a netcdf file

    Parameters
    ----------
    ds : xarray dataset
        year end summary dataset
    t : tuple, (int,int,int,str)
        the year, month, day and 'AM'|'PM'
    '''
    # get file name
    fnm = get_summary_fnm(t)

    # check folder
    check_filefolder(fnm)

    # save netcdf file
    ds.to_netcdf(fnm)

def load_summary(t):
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
    fnm = get_summary_fnm(t)

    # read data as xarray Dataset
    ds = xr.open_dataset(fnm)

    return ds

def save_summarycsv(df,year,op='heritage'):
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
    from FireConsts import dirpjdata
    import os

    fnm = os.path.join(dirpjdata, str(year),'Summary','Flist_'+op+'_'+str(year)+'.csv')
    check_filefolder(fnm)

    df.to_csv(fnm)

def read_summarycsv(year,op='heritage'):
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
    from FireConsts import dirpjdata
    import pandas as pd
    import os

    fnm = os.path.join(dirpjdata, str(year),'Summary','Flist_'+op+'_'+str(year)+'.csv')
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

def get_lts_serialization(year=None):
    ''' get the time with lastest pkl data
    '''
    from FireConsts import dirpjdata
    from datetime import date
    from glob import glob
    import os

    if year == None:
        year = date.today().year

    fnms = glob(dirpjdata+str(year)+'/Serialization/*.pkl')

    if len(fnms) > 0:
        fnms.sort()
        fnm_lts = os.path.basename(fnms[-1])

        lts = [int(fnm_lts[0:4]),int(fnm_lts[4:6]),int(fnm_lts[6:8]),fnm_lts[8:10]]
    else:
        lts = None

    return lts


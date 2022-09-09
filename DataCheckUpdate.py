""" DataUpdate
This module include functions used to check and update needed data files
"""
from FireConsts import dirextdata
import subprocess
from datetime import date, timedelta
import pandas as pd

from FireIO import os_path_exists

def glob2(pathname, **kwargs):
    '''
    Custom implementation of Glob that also works for S3 directories.
    '''
    if pathname.startswith("s3://"):
        import s3fs
        import fnmatch
        import os
        s3 = s3fs.S3FileSystem(anon=False)
        s3_dir = os.path.dirname(pathname)
        fnms = [f"s3://{f}" for f in s3.ls(s3_dir) if fnmatch.fnmatch(f"s3://{f}", pathname)]
        return sorted(fnms)
    else:
        import glob
        return sorted(list(glob.glob(pathname, **kwargs)))

# ------------------------------------------------------------------------------
# Check external dataset
# ------------------------------------------------------------------------------
def check_VNP14IMGML_avail(year, month, ver="C1.05"):
    """ Check if the monthly VIIRS data exist (C1.05 or C1.04)

    Parameters
    ----------
    year : int
        the year
    month : int
        the month
    """

    from FireConsts import dirextdata
    import os
    from datetime import date

    # derive monthly file name with path
    t = date(year, month, 1)
    dirFC = os.path.join(dirextdata, "VNP14IMGML") + "/"
    fnmFC = os.path.join(dirFC, "VNP14IMGML." + t.strftime("%Y%m") + "." + ver + ".txt")

    # check if the file exist
    if os_path_exists(fnmFC):
        print("[Yes] Monthly VNP14IMGML exists at ", dirFC)
        fnms = glob2(dirFC + "VNP14IMGML." + t.strftime("%Y") + "*.txt")
        mons = max([int(os.path.basename(fnm)[15:17]) for fnm in fnms])
        print(f"The latest VNP14IMGML available month for {year} is {mons}")
    else:
        print("[No] Monthly VNP14IMGML is not available at ", dirFC)
        # if the file does not exist, also find the last available month
        fnms = glob2(dirFC + "VNP14IMGML." + t.strftime("%Y") + "*.txt")
        if len(fnms) == 0:
            print("No any monthly VNP14IMGML is not available at ", dirFC)
        else:
            mons = max([int(os.path.basename(fnm)[15:17]) for fnm in fnms])
            print(f"The latest VNP14IMGML available month for {year} is {mons}")


def check_VNP14IMGTDL_avail(year, month, day):
    """ Check if the daily NRT VIIRS data exist

    Parameters
    ----------
    year : int
        the year
    month : int
        the month
    day : int
        the day
    """

    from FireConsts import dirextdata
    import os
    from datetime import date, timedelta

    # derive monthly file name with path
    t = date(year, month, day)
    dirFC = os.path.join(dirextdata, "VNP14IMGTDL") + "/"
    fnmFC = os.path.join(
        dirFC, "SUOMI_VIIRS_C2_Global_VNP14IMGTDL_NRT_" + t.strftime("%Y%j") + ".txt"
    )

    # check if the file exist
    if os_path_exists(fnmFC):
        print("[Yes] Daily VNP14IMGTDL exists!")
    else:
        print("[No] Daily VNP14IMGTDL is not available!")

        # if the file does not exist, also find the last available day
        fnms = glob2(
            dirFC
            + "SUOMI_VIIRS_C2_Global_VNP14IMGTDL_NRT_"
            + t.strftime("%Y")
            + "*.txt"
        )
        if len(fnms) > 0:
            days = max([int(os.path.basename(fnm)[42:45]) for fnm in fnms])
            tmax = date(t.year, 1, 1) + timedelta(days=days - 1)
            print(
                "The latest VNP14IMGTDL available date is "
                + tmax.strftime("%Y%m%d")
                + "(doy = "
                + tmax.strftime("%j")
                + ")"
            )


def check_GridMET_avail(year, month, day):
    """ Check if the GridMET fm1000 data exist

    Parameters
    ----------
    year : int
        the year
    month : int
        the month
    day : int
        the day
    """

    from FireConsts import dirextdata
    import os
    from datetime import date
    import xarray as xr
    import pandas as pd

    import warnings

    warnings.simplefilter("ignore")

    t = date(year, month, day)
    dirGridMET = os.path.join(dirextdata, "GridMET") + "/"
    fnmFM1000 = dirGridMET + "fm1000_" + t.strftime("%Y") + ".nc"

    if os_path_exists(fnmFM1000):
        print("[Yes] Annual GridMET FM1000 exists!")

        # if the file exist, check whether the date is beyond the time range in the data
        ds = xr.open_dataarray(fnmFM1000)
        tmax = pd.to_datetime(ds["day"][-1].values).date()
        print("The last time in GridMET data is " + tmax.strftime("%Y-%m-%d"))
    else:
        print("[No] Annual GridMET FM1000 is not available!")


def check_data_avail(year, month, day):
    """ A tool used to check the availability of NRT fire and FM1000 data:
    1. monthly suomi VIIRS active fire (VNP14IMGML)
    2. daily suomi VIIRS NRT active fire (VNP14IMGTDL)
    3. annual GridMET FM1000 data

    Parameters
    ----------
    year : int
        the year
    month : int
        the month
    day : int
        the day
    """

    # VIIRS VNP14IMGML data
    check_VNP14IMGML_avail(year, month)

    # VIIRS VNP14IMGTDL data
    check_VNP14IMGTDL_avail(year, month, day)

    # GridMET fm1000 data
    check_GridMET_avail(year, month, day)



# ------------------------------------------------------------------------------
# update external dataset
# ------------------------------------------------------------------------------
def wget(url,locdir=None,header=None,nc=False,recursive=False,depth=None,no_remove_listing=False,
         np=False,acclist=None,rejlist=None,nH=False,cut_dirs=None,robots_off=False,mirror=False,
         timestamping=False, no_wget=False):
    ''' Simulate the shell command wget
    Usage: wget(url,locdir=None,nc=False,recursive=False,depth=None,no_remove_listing=False,
         np=False,rejlist=None,nH=False,cut_dirs=None,header=None,robots_off=False,mirror=False)

    Parameters
    ----------
    url : str
        The remote url to be downloaded
    locdir : str (optional)
        The local directory used to save the downloaded files, default = '.'
    header : str (optional)
        The header for query. If set to 'NASA', the API for NASA earthdata will be used to authorization
    nc : bool (optional)
        If True, the previously downloaded local file will be preserved
    recursive : bool (optional)
        If Ture, turn on recursive retrieving.
    depth : str (optional)
        If not None, specify recursion maximum depth level depth ('inf' for infinite recursion)
    no_remove_listing : bool (optional)
        If True, do not remove listings downloaded by Wget.
    np : bool (optional)
        If True, do not ever ascend to the parent directory when retrieving recursively
    acclist : str (optional)
        Specify comma-separated lists of file name suffixes or patterns to accept
    rejlist : str (optional)
        Specify comma-separated lists of file name suffixes or patterns to reject
    nH : bool (optional)
        If True, disable host-prefixed file names
    cut_dirs : str of number (optional)
        The number of directory layers (starting from root) to be cut.
    robots_off : bool (optional)
        If True, turn off the robot exclusion
    mirror : bool (optional)
        If True, turn on options suitable for mirroring.
        It is currently equivalent to ‘-r -N -l inf --no-remove-listing’.
    timestamping : bool (optional)
        If True, turn on timestamping option (-N)

    Returns
    -------
    strcmd : str
        The full wget script with options
    '''

    strcmd = 'wget'

    # -------
    # options

    # '-nc'
    if nc:
        strcmd = ' '.join([strcmd, '-nc'])

    # '-r'
    if recursive:
        strcmd = ' '.join([strcmd, '-r'])

    # '-l depth'
    if depth is not None:
        strcmd = ' '.join([strcmd, '-l', depth])

    # ‘-np’
    if np:
        strcmd = ' '.join([strcmd, '-np'])

    # '--no-remove-listing'
    if no_remove_listing:
        strcmd = ' '.join([strcmd, '--no-remove-listing'])

    # '-A rejlist'
    if acclist is not None:
        strcmd = ' '.join([strcmd, '-A', acclist])

    # '-R rejlist'
    if rejlist is not None:
        strcmd = ' '.join([strcmd, '-R', rejlist])

    # ‘-nH’
    if nH:
        strcmd = ' '.join([strcmd, '-nH'])

    if cut_dirs is not None:
        strcmd = ' '.join([strcmd, '--cut-dirs='+cut_dirs])

    # '-e robots=off'
    if robots_off:
        strcmd = ' '.join([strcmd, '-e', 'robots=off'])

    # ‘-m’
    if mirror:
        strcmd = ' '.join([strcmd, '-m'])

    # '-N'
    if timestamping:
        strcmd = ' '.join([strcmd, '-N'])

    # ---
    # url
    strurl = '"'+url+'"'
    strcmd = ' '.join([strcmd, strurl])

    # -----------
    # '-P locdir'
    if locdir is not None:
        strlocdir = '-P "'+locdir+'"'
        strcmd = ' '.join([strcmd, strlocdir])

    # -------
    # '--header header'
    if header is not None:
        if header=='NASA':
            strheader = '--header "Authorization: Bearer dHJpbmtldDplV05vWlc0eE4wQjFZMmt1WldSMToxNjIyODMxODExOjg1NDMyYTZiZTFjZDFkNzZkZWIxMjc3ODdlYzY2NGUxMmI1NzYyMTU"'
        else:
            strheader = '--header "'+header+'"'
        strcmd = ' '.join([strcmd, strheader])

    # run
    if no_wget is False:
        subprocess.run(strcmd,shell=True)

    return strcmd


def update_VNP14IMGTDL(local_dir=None):

    ''' Batch read and extract update_VJ114IMGTDL data
    Usage : update_VJ114IMGTDL(local_dir)

    Parameters
    ----------
    local_dir : str
        the directory containing the downloaded data
    '''

    # The directory to save VNP14IMGTDL data
    if local_dir == None:
        local_dir=dirextdata+'VIIRS/VNP14IMGTDL/'

    # Derive the date periods needed to download
    today = date.today()
    fnms = glob2(local_dir+'SUOMI_VIIRS_C2_Global_VNP14IMGTDL_NRT_'+str(today.year)+'*.txt')
    if len(fnms)==0:
        ndays = 0
    else:
        doys = [int(d[-7:-4]) for d in fnms]
        ndays = max(doys)
    dstart = date(today.year,1,1) + timedelta(days=ndays)
    dstart = dstart - timedelta(days=1)              # downloaded the last file again to avoid incomplete data

    # Do the download process
    urldir = "https://nrt3.modaps.eosdis.nasa.gov/api/v2/content/archives/FIRMS/suomi-npp-viirs-c2/Global/"
    for d in pd.date_range(dstart,today):
        urlfnm = urldir + "SUOMI_VIIRS_C2_Global_VNP14IMGTDL_NRT_"+d.strftime('%Y%j')+".txt"
        strcmd = wget(url=urlfnm,locdir=local_dir,robots_off=True,no_wget=False,timestamping=True,header='NASA')

def update_VJ114IMGTDL(local_dir=None):

    ''' Batch read and extract update_VJ114IMGTDL data
    Usage : update_VJ114IMGTDL(local_dir)

    Parameters
    ----------
    local_dir : str
        the directory containing the downloaded data
    '''

    # The directory to save VNP14IMGTDL data
    if local_dir == None:
        local_dir=dirextdata+'VIIRS/VJ114IMGTDL/'
    # Derive the date periods needed to download
    today = date.today()
    fnms = glob2(local_dir+'J1_VIIRS_C2_Global_VJ114IMGTDL_NRT_'+str(today.year)+'*.txt')
    if len(fnms)==0:
        ndays = 0
    else:
        doys = [int(d[-7:-4]) for d in fnms]
        ndays = max(doys)
    dstart = date(today.year,1,1) + timedelta(days=ndays)
    dstart = dstart - timedelta(days=1)              # downloaded the last file again to avoid incomplete data

    # Do the download process
    urldir = "https://nrt3.modaps.eosdis.nasa.gov/api/v2/content/archives/FIRMS/noaa-20-viirs-c2/Global/"
    for d in pd.date_range(dstart,today):
        urlfnm = urldir + "J1_VIIRS_C2_Global_VJ114IMGTDL_NRT_"+d.strftime('%Y%j')+".txt"
        strcmd = wget(url=urlfnm,locdir=local_dir,robots_off=True,no_wget=False,timestamping=True,header='NASA')
        
        
def update_GridMET_fm1000(local_dir=None):
    ''' Get updated GridMET data (including fm1000)
    '''
    
    # The directory to save GridMET data
    if local_dir == None:
        local_dir = dirextdata+'GridMET/'

    today = date.today()

    # Do the download process
    urldir = "http://www.northwestknowledge.net/metdata/data/"
    # strvars = ['vpd','pr','tmmn','tmmx','vs','fm100','fm1000','bi','pdsi']
    strvars = ['fm1000']
    for strvar in strvars:
        urlfnm = urldir + strvar + '_' + str(today.year) + '.nc'
        strget = ' '.join(['wget', '-N', '-c', '-nd',urlfnm, '-P', local_dir])
        subprocess.run(strget,shell=True)


if __name__ == "__main__":
    """ Do check or update the data
    """
    # check_data_avail(year,month,day)
    check_VNP14IMGML_avail(2021, 1, ver="C1.05")

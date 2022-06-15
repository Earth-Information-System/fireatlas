""" FireObj
This is the module containing the object definitions used in the system

This module include:
1. MODULES AND UTILITY FUNCTIONS
    a. Other project modules used in this module
    b. Object related utility functions (related to time step)
2. FOUR LAYERS OF OBJECTS
    a. Allfires:  the class of all fire events
    b. Fire:      the class of a fire event
    c. Cluster:   the class of active fire pixel cluster (only for supporting)
    d. FirePixel: the class of an active fire pixel
"""

# ------------------------------------------------------------------------------
# 1. MODULES AND UTILITY FUNCTIONS
# ------------------------------------------------------------------------------

# a. Project modules used
import FireIO
import FireVector
import FireClustering
from FireConsts import maxoffdays,limoffdays,area_VI,fpbuffer,flbuffer

# b. Object-realted utility functions
# The time step in the module is now defined as a list (year, month, day, ampm).
#   The following functions are used to convert times between different formats.
#   t : time steps, tuple (year,month,day,ampm)
#   d : date, datetime.date()
#   ampm : ampm, str()
#   dt : time steps, datetime.datetime()
def t_nb(t,nb='next'):
    ''' Calculate the next or previous time step (year, month, day, ampm)
    Parameters
    ----------
    t : tuple, (int,int,int,str)
        the year, month, day and 'AM'|'PM' for present time
    nb : str, 'next'|'previous'
        option to extract next or previous time step

    Returns
    -------
    t_out : tuple, (int,int,int,str)
        the year, month, day and 'AM'|'PM' for next/previous time
    '''
    from datetime import date, timedelta

    # the next time step
    if nb == 'next':
        # if current time is 'AM', set next time as the current day and 'PM'
        if t[-1] == 'AM':
            t_out = list(t[:-1])
            t_out.append('PM')
        # if current time is 'PM', set next time as the following day and 'AM'
        else:
            d = date(*t[:-1])
            d_out = d + timedelta(days=1)
            t_out = [d_out.year,d_out.month,d_out.day,'AM']

    # the previous time step
    elif nb == 'previous':
        # if current time is 'PM', set previous time as the current day and 'AM'
        if t[-1] == 'PM':
            t_out = list(t[:-1])
            t_out.append('AM')
        # if current time is 'AM', set previous time as the previous day and 'PM'
        else:
            d = date(*t[:-1])
            d_out = d + timedelta(days=-1)
            t_out = [d_out.year,d_out.month,d_out.day,'PM']
    return t_out

def t_dif(t1,t2):
    ''' calculate the time difference between two time steps
    Parameters
    ----------
    t1 : tuple, (int,int,int,str)
        the year, month, day and 'AM'|'PM' for time 1
    t2 : tuple, (int,int,int,str)
        the year, month, day and 'AM'|'PM' for time 2

    Returns
    -------
    dt : float
        time difference in days (t2-t1), half day as 0.5
    '''
    from datetime import date

    # calculate the day difference
    d1 = date(*t1[:-1])
    d2 = date(*t2[:-1])
    dt = (d2-d1).days

    # adjust according to ampm difference
    if t1[-1] != t2[-1]:
        if t1[-1] == 'PM':
            dt -= 0.5
        else:
            dt += 0.5
    return dt

def t2d(t):
    ''' convert a t tuple to date and ampm
    Parameters
    ----------
    t : tuple, (int,int,int,str)
        the year, month, day and 'AM'|'PM'

    Returns
    -------
    d : datetime date
        date
    ampm : str, 'AM'|'PM'
        ampm indicator
    '''
    from datetime import date

    d = date(*t[:-1])     # current date, datetime date
    ampm = t[-1]          # current ampm, 'AM'|'PM'

    return d, ampm

def t2dt(t):
    ''' convert a t tuple to datetime
    Parameters
    ----------
    t : tuple, (int,int,int,str)
        the year, month, day and 'AM'|'PM'

    Returns
    -------
    dt : datetime datetime
        datetime
    '''
    from datetime import datetime
    dlh = {'AM':0,'PM':12}
    dhl = {0:'AM',12:'PM'}

    dt = datetime(*t[:-1],dlh[t[-1]])

    return dt

def d2t(year,month,day,ampm):
    ''' convert year, month, day, ampm to a t tuple
    Parameters
    ----------
    year : int
        year
    month : int
        month
    day : int
        day
    ampm : str, 'AM'|'PM'
        ampm indicator

    Returns
    -------
    t : list, (int,int,int,str)
        the year, month, day and 'AM'|'PM'
    '''
    t = [year,month,day,ampm]
    return t

def dt2t(dt):
    ''' convert datetime to a t tuple
    Parameters
    ----------
    dt : datetime datetime
        datetime
    Returns
    -------
    t : tuple, (int,int,int,str)
        the year, month, day and 'AM'|'PM'
    '''
    dlh = {'AM':0,'PM':12}
    dhl = {0:'AM',12:'PM'}

    t = [dt.year,dt.month,dt.day,dhl[dt.hour]]
    return t

def ftrange(firstday,lastday):
    ''' get datetime range for given first and last t tuples (both ends included)

    Parameters
    ----------
    firstday : tuple, (int,int,int,str)
        the year, month, day and 'AM'|'PM'
    lastday : tuple, (int,int,int,str)
        the year, month, day and 'AM'|'PM'

    Returns
    -------
    trange : pandas date range
        date range defined by firstday and lastday
    '''
    import pandas as pd

    trange = pd.date_range(t2dt(firstday),t2dt(lastday),freq='12h')
    return trange

# ------------------------------------------------------------------------------
# 2. FOUR LAYERS OF OBJECTS
# ------------------------------------------------------------------------------

# a. Object - Allfires
class Allfires:
    """ Class of allfire events at a particular time step
    """

    # initilization
    def __init__(self,t):
        ''' Initiate the object with current time
        Parameters
        ----------
        t : tuple, (int,int,int,str)
            the year, month, day and 'AM'|'PM'
        '''
        # time
        self.t = t_nb(t,nb='previous') # itialize the object with previous time step

        # a list of Fire objects
        self.fires = []

        # list of fire ids which has changes at the current time step
        self.fids_expanded = []   # a list of ids for fires with expansion at current time step
        self.fids_new = []        # a list of ids for all new fires formed at current time step
        self.fids_merged = []     # a list of ids for fires with merging at current time step
        self.fids_invalid = []    # a list of ids for fires invalidated at current time step

        # cumulative recordings
        self.heritages = []       # a list of fire heritage relationships (source, target)

    # properties
    @property
    def cday(self):
        ''' Datetime date of current time step
        '''
        from datetime import date
        return date(*self.t[:-1])

    @property
    def ampm(self):
        ''' Ampm indicator of current time step
        Parameters
        ----------
        ampm : str, 'AM'|'PM'
           the ampm option calculated from t
        '''
        return self.t[-1]

    @property
    def number_of_fires(self):
        ''' Total number of fires (active and inactive) at this time step
        '''
        return len(self.fires)

    @property
    def fids_active(self):
        ''' List of active fire ids
        '''
        return [f.fireID for f in self.fires if f.isactive is True]

    @property
    def fids_sleeper(self):
        ''' List of fire ids that may reactivate
        '''
        return [f.fireID for f in self.fires if f.mayreactivate is True]

    @property
    def number_of_activefires(self):
        ''' Total number of active fires at this time step
        '''
        return len(self.fids_active)

    @property
    def activefires(self):
        ''' List of active fires
        '''
        return [self.fires[fid] for fid in self.fids_active]

    @property
    def fids_valid(self):
        ''' List of valid (non-invalid) fire ids
        '''
        return [f.fireID for f in self.fires if f.invalid is False]

    @property
    def number_of_validfires(self):
        ''' Total number of valid fires at this time step
        '''
        return len(self.fids_valid)

    @property
    def validfires(self):
        ''' List of valid fires
        '''
        return [self.fires[fid] for fid in self.fids_valid]

    @property
    def fids_updated(self):
        ''' List of fire id which is updated at this time step
            (expanded, new, merged, invalid)
        '''
        fids_updated = list(set(self.fids_expanded+self.fids_new+
                                self.fids_merged+self.fids_invalid))
        return fids_updated

    @property
    def fids_ne(self):
        ''' List of fire id which is newly formed or expanded
               at this time step
        '''
        fids_ne = sorted(set(self.fids_expanded+self.fids_new))
        return fids_ne

    # functions to be run before tracking VIIRS active fire pixels at each time step
    def update_t(self,t):
        ''' Update the time (cday and ampm) for the Allfire object.
        Parameters
        ----------
        t : tuple, (int,int,int,str)
            the year, month, day and 'AM'|'PM'
        '''
        self.t = list(t)  # current date and ampm

    def update_t_allfires(self,t):
        ''' Update the time (t) for each Fire object in the Allfire object.
        Parameters
        ----------
        t : tuple, (int,int,int,str)
            the year, month, day and 'AM'|'PM'
        '''
        for f in self.fires:
            f.t = list(t)

    def reset_newpixels(self):
        ''' Reset newpixels to [] for each fire.
        '''
        for f in self.fires:
            f.newpixels = []

    def reset_fids_updated(self):
        ''' Reset the fids used to record changes.
        '''
        self.fids_expanded = []   # a list of ids for fires with expansion at current time step
        self.fids_new = []        # a list of ids for all new fires formed at current time step
        self.fids_merged = []     # a list of ids for fires with merging at current time step
        self.fids_invalid = []    # a list of ids for fires invalidated at current time step

    def newyear_reset(self):
        ''' reset fire ids at the start of a new year
        '''
        # re-id all active fires
        lastyearfires = []
        fidmapping = []
        nfid = 0
        for f in self.activefires:
            ofid = f.fireID
            f.fireID = nfid
            lastyearfires.append(f)
            fidmapping.append((ofid,nfid))
            nfid += 1
        self.fires = lastyearfires

        # clean heritages
        self.heritages = []

        # save the mapping table
        FireIO.save_newyearfidmapping(fidmapping,self.t[0])

    # functions to be run after tracking VIIRS active fire pixels at each time step
    def record_fids_change(self,
         fids_expanded=None, fids_new=None, fids_merged=None, fids_invalid=None):
        ''' Update the list of ids for fires with changes at current time step.
        Parameters
        ----------
        fids_expanded : list
            ids of expanded fires
        fids_new : list
            ids of new formed fires
        fids_merged : list
            ids of fires with other fires merging to them
        fids_invalid : list
            ids of fires invalidated (due to merging with other fires)
        '''
        if fids_expanded:
            self.fids_expanded = fids_expanded     # expanded fires
        if fids_new:
            self.fids_new = fids_new               # new formed fires
        if fids_merged:
            self.fids_merged = fids_merged         # fires with merging with other fires
        if fids_invalid:
            self.fids_invalid = fids_invalid       # fires invalidated due to merging with other fires

    def invalidate_statfires(self):
        ''' If pixel density of an active fire is too large, assume it's static
                fires and invalidate it.
        '''
        for f in self.activefires:
            if (f.pixden > 20) & (f.farea < 20):
                # invalidate the fire
                f.invalid = True

                # add the fire id into the fids_invalid list
                self.fids_invalid.append(f.fireID)

    # def updateLCTmax(self):
    #     ''' update Land cover type (dominant LCT for all fire pixels) for active
    #             fires that is small in size (<1000 pixels)
    #     '''
    #     import FireIO
    #     for f in self.activefires:
    #         if (f.n_pixels < 1000):
    #             # get all LCT for the fire pixels
    #             vLCT = FireIO.get_LCT(f.locs)
    #             # extract the LCT with most pixel counts
    #             LCTmax = max(set(vLCT), key = vLCT.count)
    #             f.LCTmax = LCTmax

    def updateftypes(self):
        ''' update fire type and dominant LCT for active fires
        '''
        for f in self.activefires:
            f.set_ftype()

# b. Object - Fire
class Fire:
    """ Class of a single fire event at a particular time step
    """

    # initilization
    def __init__(self, id, t, pixels, sensor='viirs'):
        ''' Initialize Fire class with active fire pixels locations and t. This
                is only called when fire clusters forming a new Fire object.
        Parameters
        ----------
        id : int
            fire id
        t : tuple, (int,int,int,str)
            the year, month, day and 'AM'|'PM'
        pixels : list (nx5)
            latitude, longitude, line, sample, and FRP values of active fire pixels
        '''
        # initialize fire id
        self.fireID = id
        self.mergeid = id
        self.sensor = sensor

        # initialize current time, fire start time, and fire final time
        tlist = list(t)
        self.t = tlist  # current time
        self.t_st = tlist
        self.t_ed = tlist

        # initialize pixels
        # fpixels = [FirePixel((p[0],p[1]),(p[2],p[3],p[4]),tlist,id) for p in pixels]
        fpixels = pixels
        self.pixels = fpixels      # all pixels
        self.newpixels = fpixels   # new detected pixels
        self.actpixels = fpixels   # new detected pixels of last active fire detection
        self.ignpixels = fpixels   # pixels at ignition

        # initialize hull using the pixels
        locs = [p.loc for p in fpixels]  # list of [lat,lon]
        hull = FireVector.cal_hull(locs, sensor) # the hull from all locs
        self.hull = hull   # record hull

        # initialize the exterior pixels (pixels within the inward extbuffer of
        #    the hull; used for saving time for hull calculation of large fires)
        self.extpixels = FireVector.cal_extpixels(fpixels,hull)

        # always set validate at initialization
        self.invalid = False

        # get and record fm1000 value at ignition (at the time of initilization)
        # self.stFM1000 = FireIO.get_stFM1000(hull,locs,tlist)

        # LCTmax
        # vLCT = FireIO.get_LCT(locs)
        # self.LCTmax = max(set(vLCT), key = vLCT.count)

    # properties
    @property
    def cday(self):
        ''' Current day (datetime date)
        '''
        from datetime import date
        return date(*self.t[:-1])

    @property
    def cdoy(self):
        ''' Current day (datetime date)
        '''
        return self.cday.timetuple().tm_yday

    @property
    def ampm(self):
        ''' Current ampm flag, 'AM'|'PM'
        '''
        return self.t[-1]

    @property
    def duration(self):
        ''' Time difference between first and last active fire detection
        '''
        duration = t_dif(self.t_st,self.t_ed) + 0.5
        return duration

    @property
    def t_inactive(self):
        ''' Time difference between current time and the last active fire detection
        '''
        t_inactive = t_dif(self.t_ed,self.t)
        return t_inactive

    @property
    def isactive(self):
        ''' Fire active status
        '''
        # invalidated fires are always inactive
        if self.invalid:
            return False
        # otherwise, set to True if no new pixels detected for 5 consecutive days
        return (self.t_inactive <= maxoffdays)

    @property
    def mayreactivate(self):
        ''' Fire active status
        '''
        # invalidated fires are always inactive
        if self.invalid:
            return False
        # otherwise, set to True if no new pixels detected for 5 consecutive days
        return (maxoffdays < self.t_inactive <= limoffdays)

    @property
    def isignition(self):
        ''' Is the current timestep the ignition?
        '''
        ign = (t_dif(self.t_st,self.t_ed) == 0)*1
        return ign

    @property
    def locs(self):
        ''' List of fire pixel locations (lat,lon)
        '''
        return [p.loc for p in self.pixels]

    @property
    def n_pixels(self):
        ''' Total number of fire pixels'''
        return len(self.pixels)

    @property
    def newlocs(self):
        ''' List of new fire pixels locations (lat,lon)
        '''
        return [p.loc for p in self.newpixels]

    @property
    def newlocsMP(self):
        ''' MultiPoint shape of newlocs
        '''
        from shapely.geometry import Point
        import geopandas as gpd
        mp = [Point(p[1],p[0]) for p in self.newlocs]
        return gpd.GeoSeries(mp)

    @property
    def newpixelatts(self):
        ''' List of new fire pixels locations (lat,lon)
        '''
        return [(p.frp, p.origin) for p in self.newpixels]

    @property
    def n_newpixels(self):
        ''' Total number of new fire pixels
        '''
        return len(self.newpixels)

    @property
    def extlocs(self):
        ''' List of exterior fire pixel locations (lat,lon)
        '''
        return [p.loc for p in self.extpixels]

    @property
    def n_extpixels(self):
        ''' Total number of exterior fire pixels
        '''
        return len(self.extpixels)

    @property
    def ignlocs(self):
        ''' List of fire pixel locations (lat,lon) at ignition time step
        '''
        return [p.loc for p in self.ignpixels]

    @property
    def n_ignpixels(self):
        ''' Total number of ignition fire pixels
        '''
        return len(self.ignpixels)

    @property
    def farea(self):
        ''' Fire spatail size of the fire event (km2)
        '''
        # get hull
        fhull = self.hull

        # If no hull, return area calculated from number of pixels
        if fhull is None:
            return self.n_pixels * area_VI
        # otherwise, use calConcHarea to calculate area,
        #   but no smaller than area_VI (sometimes calculated hull area is very mall)
        else:
            return max(FireVector.calConcHarea(fhull),area_VI)

    @property
    def pixden(self):
        ''' Fire pixel density (number of pixels per km2 fire area)
        '''
        farea = self.farea
        if farea > 0:
            return self.n_pixels/farea
        else:
            return 0

    @property
    def meanFRP(self):
        ''' Mean FRP of the new fire pixels
        '''
        frps = [p.frp for p in self.newpixels]
        if len(frps) > 0:
            m = sum(frps)/len(frps)
        else:
            m = 0
        return m

    @property
    def centroid(self):
        ''' Centroid of fire object (lat,lon)
        '''
        # get hull
        fhull = self.hull

        if fhull is not None: # centroid of the hull
            cent = (fhull.centroid.y,fhull.centroid.x)
        else: # when no hull, use the centroid of all pixels
            cent = FireClustering.cal_centroid(self.locs)
        return cent

    # @property
    # def ftype(self):
    #     ''' Fire type (as defined in FireConsts) derived using LCTmax and stFM1000
    #     '''
    #     return 2 # set to wild forest type, we'll deal with fire type later
        # get the dominant land cover type
        # LCTmax = self.LCTmax

        # determine the fire type using the land cover type and stFM1000
        # if LCTmax in [11,31]:   # 'Water', 'Barren' -> 'Other'
        #     return 0
        # elif LCTmax in [23]:    # 'Urban' -> 'Urban'
        #     return 1
        # elif LCTmax in [82]:    # 'Agriculture' -> 'Agriculture'
        #     return 6
        # elif LCTmax in [42]:    # 'Forest' ->
        #     stFM1000 = self.stFM1000
        #     if stFM1000 > 12:        # 'Forest manage'
        #         return 3
        #     else:                  # 'Forest wild'
        #         return 2
        # elif LCTmax in [52,71]:    # 'Shurb', 'Grassland' ->
        #     stFM1000 = self.stFM1000
        #     if stFM1000 > 12:        # 'Shrub manage'
        #         return 5
        #     else:                  # 'Shrub wild'
        #         return 4

    def set_ftype(self):
        ''' set fire type and dominant LCT for newly activated fires
        '''
        from FireConsts import FTYP_opt
        # 0 - use preset ftype (defined as FTYP_preset in FireConsts) for all fires
        # 1 - use the CA type classification (dtermined using LCTmax)
        # 2 - use the global type classfication (need more work...)

        if FTYP_opt == 0: # use preset ftype for all fires
            from FireConsts import FTYP_preset
            self.ftype = FTYP_preset
        elif FTYP_opt == 1: # use CA type classifications (determined using LCTmax)
            # update or read LCTmax
            # call get_LCT to get all LCT for the fire pixels
            vLCT = FireIO.get_LCT(self.locs)
            # extract the LCT with most pixel counts
            LCTmax = max(set(vLCT), key = vLCT.count)
            self.LCTmax = LCTmax

            # get and record fm1000 value at ignition (at the time of initilization)
            self.stFM1000 = FireIO.get_stFM1000(self.hull,self.locs,self.t_st)

            # determine the fire type using the land cover type and stFM1000
            if LCTmax in [11,31]:   # 'Water', 'Barren' -> 'Other'
                self.ftype = 0
            elif LCTmax in [23]:    # 'Urban' -> 'Urban'
                self.ftype = 1
            elif LCTmax in [82]:    # 'Agriculture' -> 'Agriculture'
                self.ftype = 6
            elif LCTmax in [42]:    # 'Forest' ->
                stFM1000 = self.stFM1000
                if stFM1000 > 12:        # 'Forest manage'
                    self.ftype = 3
                else:                  # 'Forest wild'
                    self.ftype = 2
            elif LCTmax in [52,71]:    # 'Shurb', 'Grassland' ->
                stFM1000 = self.stFM1000
                if stFM1000 > 12:        # 'Shrub manage'
                    self.ftype = 5
                else:                  # 'Shrub wild'
                    self.ftype = 4
        elif FTYP_opt == 2:  # global type classification
            self.ftype = 1  # need more work here...


    @property
    def ftypename(self):
        ''' Fire type name
        '''
        from FireConsts import FTYP
        return FTYP[self.ftype]

    @property
    def fperim(self):
        ''' Perimeter length of fire hull
        '''
        # get hull
        fhull = self.hull

        if fhull is None:  # if no hull, return zero
            perim = 0
        else:  # otherwise, use the hull length
            perim = fhull.length
        return perim

    @property
    def flinepixels(self):
        ''' List of all fire pixels near the fire perimeter (fine line pixels)
        '''
        from shapely.geometry import Point, MultiLineString

        # get pixels of last active fire detection
        nps = self.actpixels

        # get hull
        fhull = self.hull

        if fhull is None: # if no hull, return empty list
            return []
        else:  # otherwise, extract the pixels nearl the hull
            # if hull is a polygon, return new pixels near the hull
            if fhull.type == 'Polygon':
                # lr = fhull.exterior.buffer(fpbuffer)
                lr = FireVector.addbuffer(fhull.exterior, fpbuffer)
                return [p for p in nps if lr.contains(Point(p.loc[1],p.loc[0]))]

            # if hull is a multipolygon, return new pixels near the hull
            elif fhull.type == 'MultiPolygon':
                # mlr = MultiLineString([x.exterior for x in fhull]).buffer(fpbuffer)
                mlr = MultiLineString([x.exterior for x in fhull])
                mlr = FireVector.addbuffer(mlr,fpbuffer)
                return [p for p in nps if mlr.contains(Point(p.loc[1],p.loc[0]))]

    @property
    def fline(self):
        ''' Active fire line MultiLineString shape (segment of fire perimeter with active fires nearby)
        '''
        from shapely.geometry import MultiLineString#,Polygon,Point,MultiPoint
        from FireConsts import valpha

        if len(self.flinepixels)==0: # this happens is last active pixels are within the fire scar
            return None

        # get fireline pixel locations
        flinelocs = [p.loc for p in self.flinepixels]
        flinelocsMP = FireVector.doMultP(flinelocs,valpha)

        # get the hull
        fhull = self.hull

        # calculate the fire line
        if fhull is None: # if no hull, return None
            return None
        else:  # otherwise, create shape of the active fire line
            if fhull.type == 'MultiPolygon':
                # extract exterior of fire perimeter
                mls = MultiLineString([plg.exterior for plg in fhull])
                # return the part which intersects with  bufferred flinelocsMP
                # return mls.intersection(flinelocsMP.buffer(flbuffer))
                flinelocsMP_buf = FireVector.addbuffer(flinelocsMP,flbuffer)
                return mls.intersection(flinelocsMP_buf)

            elif fhull.type == 'Polygon':
                mls = fhull.exterior
                # return mls.intersection(flinelocsMP.buffer(flbuffer))
                flinelocsMP_buf = FireVector.addbuffer(flinelocsMP,flbuffer)
                return mls.intersection(flinelocsMP_buf)
            else:  # if fhull type is not 'MultiPolygon' or 'Polygon', return flinelocsMP
                return flinelocsMP

    @property
    def flinelen(self):
        ''' The length of active fire line
        '''
        try:
            flinelen = self.fline.length
        except:
            flinelen = 0

        return flinelen

# c. Object - Cluster
class Cluster:
    """ class of active fire pixel cluster at a particular time
    """

    # initilization
    def __init__(self, id, pixels, t, sensor='viirs'):
        ''' initilization

        Parameters
        ----------
        id : int
            cluster id number
        pixels : 3-element list
            (lat, lon, FRP) of AF pixels
        t : tuple, (int,int,int,str)
            the year, month, day and 'AM'|'PM'
        '''
        from datetime import date
        self.cday = date(*t[:-1])  # current date
        self.ampm = t[-1]          # current ampm
        self.id = id
        self.sensor = sensor
        self.pixels = pixels       # (lat,lon, FRP)

    # properties
    @property
    def locs(self):
        ''' List of pixel locations (lat,lon)
        '''
        return [(p.lat,p.lon) for p in self.pixels]

    @property
    def centroid(self):
        ''' Centroid of the cluster (lat, lon)
        '''
        return FireClustering.cal_centroid(self.locs)

    @property
    def n_pixels(self):
        ''' Number of total pixels
        '''
        return len(self.pixels)

    @property
    def hull(self):
        ''' Fire concave hull (alpha shape)
        '''
        hull = FireVector.cal_hull(self.locs, self.sensor)
        return hull

    @property
    def b_box(self):
        ''' Bounding box of concave hull
        '''
        b_box = self.hull.bounds
        return b_box

# d. Object - FirePixel
class FirePixel:
    """ class of an acitve fire pixel, which includes
        loc : location (lat, lon)
        atts : line & sample or viirs pixel, fire radiative power
        t : time (y,m,d,ampm) of record
        origin : the fire id originally recorded (before merging)
    """
    def __init__(self,lat,lon,frp,t,origin):
        self.lat = lat
        self.lon = lon
        self.frp = frp        # frp
        self.t = list(t)      # (year,month,day,ampm)
        self.origin = origin  # originate  fire id

    @property
    def loc(self):
        return (self.lat,self.lon)

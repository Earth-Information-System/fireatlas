""" FireGpkg
Module for creating regional geopackage summary at each time step
List of geojson types
---------------------
* fperim(''): fire basic attributes and perimeter geometry
* fline('FL'): active fire line geometry
* NFP('NFP'): new fire pixels
List of functions
-----------------
* make_gdf_fperim
* make_gdf_fline
* make_gdf_NFP
* save_gdf_1t
* save_gdf_trng
Modules required
----------------
* FireObj
* FireIO
* FireConsts
"""

def getdd(layer):
    ''' Get attributes names and formats for different gpkg layers
    '''
    if layer == "perimeter":
        dd = {
            "mergeid": "int",  # this is the id in the large fire database
            "ftype": "int",  # fire type
            "n_pixels": "int",  # number of total pixels
            "n_newpixels": "int",  # number of new pixels
            "farea": "float",  # fire size
            "fperim": "float",  # fire perimeter length
            "flinelen": "float",  # active fire front line length
            "duration": "float",  # fire duration
            "pixden": "float",  # fire pixel density
            "meanFRP": "float",  # mean FRP of the new fire pixels
            "t": "datetime64[ns]",
            "t_st": "datetime64",
            "t_ed": "datetime64",
        }
    elif layer in ["fireline", "newfirepix"]:
        dd = {
            "mergeid": "int",
            "t": "datetime64[ns]",
            "t_st": "datetime64[ns]",
            "t_ed": "datetime64[ns]",
        }
    return dd



def make_gdf_snapshot(allfires, regnm, layer="perimeter"):
    """ Create gpd DataFrame for fire basic attributes and fire perimeter.
    Method:
    -------
    Extract attributes rom the allfires obj input, and create a gdf DataFrame.
    In order to save the running time, if the gpkg file for previous time step
    is available, read the data and only modify properties of previously active
    fires.
    Parameters
    ----------
    allfires : Allfires object
        the Allfires object to be used to modify gdf
    Returns
    -------
    gdf : geopandas DataFrame
        the gdf containing half daily fire basic attributes and fire perimeter
    """
    import numpy as np
    import geopandas as gpd
    import pandas as pd
    # from FireConsts import dd
    import FireIO, FireObj, FireTime
    from FireConsts import epsg

    # diagnostic data name and types saved in geojson files (in addition to geometries)
    if layer == "perimeter":
        # dd1: group of attributes that can change with time (even without expansion or merging)
        dd1 = {
            "isactive": "int",  # active status
            "isdead": "int",  # dead status
            "t_inactive": "int",  # how long has it been inactive
            "isignition": "int",  # is this a new ignition?
            "mayreactivate": "int",  # sleeper status
            "t": "datetime64",
        }

        # dd2: group of attributes that can only change with fire expansion or merging
        dd2 = {
            "mergeid": "int",  # this is the id in the large fire database
            "ftype": "int",  # fire type
            # 'invalid':'int',              # invalid status
            "n_pixels": "int",  # number of total pixels
            "n_newpixels": "int",  # number of new pixels
            "farea": "float",  # fire size
            "fperim": "float",  # fire perimeter length
            "flinelen": "float",  # active fire front line length
            "duration": "float",  # fire duration
            "pixden": "float",  # fire pixel density
            "meanFRP": "float",  # mean FRP of the new fire pixels
            "t_st": "datetime64",
            "t_ed": "datetime64",
        }
    elif layer == "fireline":
        dd1 = {
            "t": "datetime64",
        }
        dd2 = {
            "mergeid": "int",
            "t_st": "datetime64",
            "t_ed": "datetime64",# this is the id in the large fire database
        }
    elif layer == "newfirepix":
        dd1 = {
            "t": "datetime64",
        }
        dd2 = {
            "mergeid": "int",
            "t_st": "datetime64",
            "t_ed": "datetime64",# this is the id in the large fire database
        }

    dd = {**dd1, **dd2}  # or (dd1 | dd2) for Python 3.9.0 or higher
    # read or initialize the gdf
    t_pt = FireTime.t_nb(allfires.t, nb="previous")
    
    gdf = FireIO.load_gpkgobj(t_pt, regnm, layer=layer)  # try to read data at t_pt

    if gdf is None:  # when no previous time step data, initialize the GeoDataFrame
        gdf = gpd.GeoDataFrame(
                columns=(list(dd.keys()) + ["fireID"]), crs="epsg:" + str(epsg), geometry=[]
            )

        gdf = gdf.set_index("fireID") # set fid column as the index column
        

    # 1. drop rows for fires that are newly invalidated or dead

    # - for newly invalidated fire objects, drop the row
    # if 'invalid' in dd.keys():
    #     for fid in allfires.fids_invalid:  # set newly invalidated object
    #         gdf.loc[fid,'invalid'] = 1
    #     gdf.drop(gdf.index[gdf['invalid'] == 1], inplace = True)  # drop the row
    gdf.drop(allfires.fids_invalid, inplace=True)

    # - for newly dead fires, drop the row
    # if 'isdead' in dd.keys():
    #     for fid,f in allfires.deadfires.items():  # set dead fire object
    #         gdf.loc[fid,'isdead'] = 1
    #     gdf.drop(gdf.index[gdf['isdead'] == 1], inplace = True)  # drop the row
    gdf.drop(allfires.fids_dead, inplace=True)

    # - drop fires with fid not in allfires object anymore (mostly fires turning from mayactive to dead)
    fids_remove = list(set(gdf.index) - set(allfires.fids))
    gdf.drop(fids_remove, inplace=True)

    # 2. modify dd2 for fires with possible modification

    # - for mayactive fires (active+sleeper), copy attributes from fire object to gdf
    for fid, f in allfires.mayactivefires.items():  # loop over active fires
        for k, tp in dd2.items():
            if tp == "datetime64":
                gdf.loc[fid, k] = FireTime.t2dt(getattr(f, k))
            else:
                gdf.loc[fid, k] = getattr(f, k)

    # update the hull of each active fire as the geometry column
    if layer == "perimeter":
        for fid, f in allfires.mayactivefires.items():
            fhull = f.hull
            if fhull.geom_type == "MultiPolygon":
                gdf.loc[fid, "geometry"] = gpd.GeoDataFrame(
                    geometry=[fhull]
                ).geometry.values
                
            else:
                gdf.loc[fid, "geometry"] = fhull

    elif layer == "fireline":
        for fid, f in allfires.activefires.items():  # loop over active fires
            fline = f.fline
            if fline is not None:
                if fline.geom_type == "MultiLineString":
                    gdf.loc[fid, "geometry"] = gpd.GeoDataFrame(
                        geometry=[fline]
                    ).geometry.values
                else:
                    gdf.loc[fid, "geometry"] = fline

    elif layer == "newfirepix":
        for fid, f in allfires.activefires.items():  # loop over active fires
            if f.n_newpixels > 0:
                gdf.loc[fid, "geometry"] = gpd.GeoDataFrame(
                    geometry=[f.newlocsMP]
                ).geometry.values

    # 3. modify dd1 for all fires
    for fid, f in allfires.fires.items():  # loop over all fires
        for k, tp in dd1.items():
            if tp == "datetime64":
                gdf.loc[fid, k] = FireTime.t2dt(getattr(f, k))
            else:
                gdf.loc[fid, k] = getattr(f, k)

    # 4. force the correct dtypes
    for k, tp in dd.items():
        gdf[k] = gdf[k].astype(tp)

    # # 5. drop the columns with no use
    # if 'invalid' in dd.keys(): gdf = gdf.drop(columns='invalid')
    
    # Set region col
    gdf['region'] = str(regnm)
    
    # Set primary key col
    t = allfires.t
    ampm = t[-1]
    if ampm == 'AM':
        time = pd.to_datetime(str(t[0])+'-'+str(t[1])+'-'+str(t[2])+'T00:00:00')
    else: 
        time = pd.to_datetime(str(t[0])+'-'+str(t[1])+'-'+str(t[2])+'T12:00:00')
    
    # primary key is: region + fireID + 12hr slice 
    gdf['primarykey'] = gdf['region'] + '|' + gdf.index.map(str) + '|' + time.isoformat()
    
    if layer == 'perimeter': # apply filter flag on the perimeter layer
        gdf['geom_counts'] = gdf["geometry"].explode(index_parts=True).groupby(['fireID']).nunique() # count number of polygons
        gdf['low_confidence_grouping'] = np.where(gdf['geom_counts']>5, 1, 0) # if more than 5 geometries are present, flag it
    
    return gdf


def save_gdf_1t(t, regnm):
    """ Creat gdf using one Allfires object and save it to a geopackage file at 1 time step.
            This can be used  for fperim, fline, and NFP files.
    Parameters
    ----------
    allfires : Allfires object
        the Allfires object to be used to create gdf
    """
    import FireIO

    # read Allfires object from the saved pkl file
    allfires = FireIO.load_fobj(t, regnm, activeonly=True)
#     allfires = FireIO.load_fobj(t, regnm, activeonly=False)

    # create gdf using previous time step gdf values and new allfires object
    gdf_fperim = make_gdf_snapshot(allfires, regnm, layer="perimeter")
    gdf_fline = make_gdf_snapshot(allfires, regnm, layer="fireline")
    gdf_nfp = make_gdf_snapshot(allfires, regnm, layer="newfirepix")

    FireIO.save_gpkgobj(
        allfires.t, regnm, gdf_fperim=gdf_fperim, gdf_fline=gdf_fline, gdf_nfp=gdf_nfp
    )


def save_gdf_trng(tst, ted, regnm):
    """ Wrapper to create and save gpkg files for a time period
    Parameters
    ----------
    tst : tuple, (int,int,int,str)
        the year, month, day and 'AM'|'PM' at start time
    ted : tuple, (int,int,int,str)
        the year, month, day and 'AM'|'PM' at end time
    """
    import FireObj, FireIO, FireTime

    # loop over all days during the period
    endloop = False  # flag to control the ending olf the loop
    t = list(tst)  # t is the time (year,month,day,ampm) for each step
    while endloop == False:
        print("Snapshot making", t)

        # create and save gdfs according to input options
        save_gdf_1t(t, regnm)

        # time flow control
        #  - if t reaches ted, set endloop to True to stop the loop
        if FireTime.t_dif(t, ted) == 0:
            endloop = True

        #  - update t with the next time stamp
        t = FireTime.t_nb(t, nb="next")


def save_gdf_uptonow(t, regnm):
    """ Create and save all valid fires (active, sleeper, and dead) up to now
    """
    import FireIO, FireTime
    import geopandas as gpd
    from FireConsts import epsg

    # read allfires object
    allfires = FireIO.load_fobj(t, regnm, activeonly=False)

    # initialize gdf
    dd = {
        "isactive": "int",  # active status
        "isdead": "int",  # dead status
        "t_inactive": "int",  # how long has it been inactive
        "isignition": "int",  # is this a new ignition?
        "mayreactivate": "int",  # sleeper status
        "t": "datetime64",
        "mergeid": "int",  # this is the id in the large fire database
        "ftype": "int",  # fire type
        "n_pixels": "int",  # number of total pixels
        "farea": "float",  # fire size
        "fperim": "float",  # fire perimeter length
        "duration": "float",  # fire duration
        "pixden": "float",  # fire pixel density
        "t_st": "datetime64",
        "t_ed": "datetime64",
    }
    gdf = gpd.GeoDataFrame(
        columns=(list(dd.keys()) + ["fireID"]), crs="epsg:" + str(epsg), geometry=[]
    )
    gdf = gdf.set_index("fireID")

    # loop over valid fires
    for fid, f in allfires.validfires.items():
        # assign attributes
        for k, tp in dd.items():
            if tp == "datetime64":
                gdf.loc[fid, k] = FireTime.t2dt(getattr(f, k))
            else:
                gdf.loc[fid, k] = getattr(f, k)
        # assign geometry
        fhull = f.hull
        if fhull.geom_type == "MultiPolygon":
            gdf.loc[fid, "geometry"] = gpd.GeoDataFrame(
                geometry=[fhull]
            ).geometry.values
        else:
            gdf.loc[fid, "geometry"] = fhull

    # make sure the data types are correct
    for k, tp in dd.items():
        gdf[k] = gdf[k].astype(tp)

    # save
    FireIO.save_gpkgobj(allfires.t, regnm, gdf_uptonow=gdf)
    # return gdf


if __name__ == "__main__":
    """ The main code to record daily geojson data for a time period
    """

    import time

    t1 = time.time()
    # set the start and end time
    # tst=(2021,7,13,'AM')
    # ted=(2021,9,15,'PM')
    #
    # # for each day during the period, create and save geojson summary file
    # save_gdf_trng(tst=tst,ted=ted,regnm='Dixie')

    # tst=(2020,10,29,'PM')
    tst = (2020, 1, 1, "PM")
    ted = (2020, 12, 31, "PM")
    # save_gdf_trng(tst=tst,ted=ted,regnm='Creek')
    save_gdf_trng(tst=tst, ted=ted, regnm='California32610')
#     save_gdf_uptonow(ted, "Creek")




    t2 = time.time()
    print(f"{(t2-t1)/60.} minutes used to run code")

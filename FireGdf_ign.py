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

def find_all_ign(tst, ted):
    ''' find all ignition points in the half-daily snapshots and
    save them to one gdf

    Parameters
    ----------
    tst : Start date of fire season
    ted: End date of fire season

    Returns
    -------
    gdf_all : geopandas DataFrame
        the gdf containing half daily fire basic attributes and ignition perimeter
    '''
    import FireObj,FireIO

    endloop = False  # flag to control the ending of the loop
    t = list(tst)    # t is the time (year,month,day,ampm) for each step
    while endloop == False:
        #print(t)
        # read daily gdf
        gdf = FireIO.load_gdfobj(t)
        gdf_ign = gdf[gdf.isignition == 1] # here we need to enter merge id instead of index

        # append daily row to gdf_all
        if FireObj.t_dif(t,tst)==0:
            gdf_all = gdf_ign
        else:
            gdf_all = gdf_all.append(gdf_ign)

        #  - if t reaches ted, set endloop to True to stop the loop
        if FireObj.t_dif(t,ted)==0:
            endloop = True

        #  - update t with the next time stamp
        t = FireObj.t_nb(t,nb='next')

    return gdf_all

def save_gdf_trng(tst,ted):
    ''' Wrapper to create and save all ignitions as gdf

    Parameters
    ----------
    tst : tuple, (int,int,int,str)
        the year, month, day and 'AM'|'PM' at start time
    ted : tuple, (int,int,int,str)
        the year, month, day and 'AM'|'PM' at end time
    '''
    import FireIO
    gdf_ign = find_all_ign(tst, ted)
    FireIO.save_gdfobj_ign(gdf_ign,tst)
        


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
    save_gdf_trng(tst=tst,ted=ted)

    t2 = time.time()
    print(f'{(t2-t1)/60.} minutes used to run code')

# -*- coding: utf-8 -*-
"""
Created on Fri Jan 28 10:52:41 2022

@author: rebec
"""

def find_all_end(tst, ted):
    ''' find all final perimeters points in the half-daily snapshots and
    save them to one gdf

    Parameters
    ----------
    tst : Start date of fire season
    ted: End date of fire season

    Returns
    -------
    gdf_all : geopandas DataFrame
        the gdf containing half daily fire basic attributes and final perimeter
    '''
    import FireObj,FireIO
    
    # initialise list of already checked fire ids
    checked_ids = []
    
    endloop = False  # flag to control the ending of the loop
    creategdf = True # needed since ted cannot be used anymore
    t = list(ted)    # t is the time (year,month,day,ampm) for each step
    while endloop == False:
        #print(t)
        # read daily gdf
        if FireIO.check_gdfobj(t,op=''):
            gdf = FireIO.load_gdfobj(t)
            gdf_active = gdf[gdf.isactive == 1] # here we need to enter merge id instead of index
    
            # append daily row to gdf_all
            if creategdf:
                gdf_all = gdf_active
                creategdf = False
            else:
                # excluse ids that have an older active perimeter
                gdf_active = gdf_active[~gdf_active['mergid'].isin(checked_ids)]
                gdf_all = gdf_all.append(gdf_active)
            
            # add ids that have been written out to checked_ids, these will be skipped next time
            checked_ids.extend(gdf_active['mergid'].tolist())
            checked_ids = list(set(checked_ids))
        
        #  - if t reaches tst, set endloop to True to stop the loop
        if FireObj.t_dif(t,tst)==0:
            endloop = True

        #  - update t with the previous time stamp
        t = FireObj.t_nb(t,nb='previous')

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
    gdf = find_all_end(tst, ted)
    FireIO.save_gdfobj(gdf,tst,param='final')
        


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
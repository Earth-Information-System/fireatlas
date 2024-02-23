import FireMain, FireGpkg, FireGpkg_sfs, FireConsts, FireIO
import os
import geopandas as gpd


def CARB_filterbyfire(perimeter_gdf_path, year, tst=None, ted=None):
    """
    Runs fire_forward within the specified time range and constrains viirs pixels to the supplied fire perimeters.
    Accomodates multiple perimeters (will loop through each one) or a single perimeter.
    Time window will be further filtered based on the fire start/end dates.

    Parameters
    ----------
    tst : tuple
        Min start time (year, month, day, "AM" or "PM").
    ted : tuple
        Max end time (year, month, day, "AM" or "PM").
    perimeter_gdf : geopandas.GeoDataFrame of fire perimeters

    Returns
    -------
    Currently spits out the pkl files in the output directory, 
    but will eventually return gdf/gpkg of sub-daily perimeters + pixels for each fire.
    Waiting to figure out the best way to do that based on the developments from Julia. 
    """

    import FireMain, FireConsts, FireIO
    import os
    import time
    from FireLog import logger

    t1 = time.time()

    # set the start and end time
    if tst is None:
        tst = (year, 1, 1, 'AM')
    if ted is None:
        ted = (year, 12, 31, 'PM')

    # read perimeter data + preprocess
    perimeter_gdf = gpd.read_file(perimeter_gdf_path)
    perimeter_gdf = FireIO.preprocess_polygon(perimeter_gdf)

    # if there's more than one perimeter, loop through each one
    # wrap fire_forward in a try/except block to catch any errors and continue processing
    def Fire_Forward_try(tst, ted, polygon):
        try:
            FireMain.Fire_Forward(tst=tst, ted=ted, polygon=polygon)
        except Exception as e:
            logger.error(e)
            pass

    if len(perimeter_gdf) > 1:
        perimeter_gdf.apply(lambda x: Fire_Forward_try(tst=tst, ted=ted, polygon=x), axis=1)
    else:
        FireMain.Fire_Forward(tst=tst, ted=ted, polygon=perimeter_gdf)
    
    # get total time used
    t2 = time.time()
    print(f'{(t2-t1)/60.} minutes used to run code')

# perimeters_fnm = os.path.join(FireConsts.dirextdata, 'perimeters', 'frap_select.geojson')


if __name__ == "__main__":
    ''' The main code to run time forwarding for a time period
    '''
    perimeters_fnm = os.path.join(FireConsts.dirextdata, 'perimeters', 'frap_select.geojson')
    CARB_filterbyfire(perimeters_fnm, 2020, tst=(2020, 9, 1, 'AM'), ted=(2020, 9, 10, 'PM'))



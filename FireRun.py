""" FireRun
Module to control different runs
"""
def Yearbatchrun(year,tst=None,ted=None,restart=False):
    ''' Run the code for each single year
    '''
    import FireMain,FireGdf,FireGdf_sfs,FireSummary

    import time

    t1 = time.time()
    # set the start and end time
    if tst is None: tst = (year,6,1,'AM')
    if ted is None: ted = (year,8,31,'PM')
    # if year == 2012: tst = (year,1,20,'AM')

    # Run the time forward and record daily fire objects .pkl data and fire attributes .GeoJSON data
    FireMain.Fire_Forward(tst=tst,ted=ted,restart=restart)
    t2 = time.time()
    print(f'{(t2-t1)/60.} minutes to run algorithm')

    # Run to save geojson files for each time step
    FireGdf.save_gdf_trng(tst=tst,ted=ted,fperim=True)
    t3 = time.time()
    print(f'{(t3-t2)/60.} minutes to save geojson files')

    # Run to save large fire geosjon files for each time step
    FireGdf_sfs.save_gdf_trng(tst=tst,ted=ted,fperim=True)
    t4 = time.time()
    print(f'{(t4-t3)/60.} minutes to save large fires')

    # Run to save year end summary and file lists
    FireSummary.save_sum_1d(tst,ted)
    FireSummary.add_heritage(ted)
    FireSummary.add_largefirelist(ted)
    t5 = time.time()
    print(f'{(t5-t4)/60.} minutes to generate summary')

    # Run to clean up large fire geojson files not needed
    FireGdf_sfs.yrend_clean(ted)

    t2 = time.time()
    print(f'{(t2-t1)/60.} minutes to run code')

if __name__ == "__main__":
    ''' The main code to run time forwarding for a time period
    '''
    import sys
    sys.path.insert(1, 'D:/fire_atlas/1_code/fireatlas')
    Yearbatchrun(2014)

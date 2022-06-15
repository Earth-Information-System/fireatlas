""" FireRun
Module to control different runs
"""
# def Yearbatchrun(year,tst=None,ted=None,restart=False):
#     ''' Run the code for each single year
#     '''
#     import FireMain,FireSummary,FireGdf_merge,FireGdf_sfs_merge,FireGdf_ign,FireGdf_final
#
#     import time
#
#     t1 = time.time()
#     # set the start and end time
#     if tst is None: tst = (year,6,1,'AM')
#     if ted is None: ted = (year,8,31,'PM')
#     # if year == 2012: tst = (year,1,20,'AM')
#
#     # Run the time forward and record daily fire objects .pkl data and fire attributes .GeoJSON data
#     FireMain.Fire_Forward(tst=tst,ted=ted,restart=restart,region='AK')
#     t2 = time.time()
#     print(f'{(t2-t1)/60.} minutes to run algorithm')
#
#     # Run to save geojson files for each time step
#     FireGdf_merge.save_gdf_trng(tst=tst,ted=ted,fperim=True)
#     #FireGdf_merge.save_gdf_trng(tst=tst,ted=ted,fall=True)
#     FireGdf_merge.save_gdf_trng(tst=tst,ted=ted,NFP_txt=True)
#     t3 = time.time()
#     print(f'{(t3-t2)/60.} minutes to save gpkg files')
#
#     # Run to save ignition point layer for each time step
#     FireGdf_ign.save_gdf_trng(tst,ted)
#     FireGdf_final.save_gdf_trng(tst,ted)
#     t31 = time.time()
#     print(f'{(t31-t3)/60.} minutes to save ognitions and final perimeters')
#
#     # Run to save large fire geosjon files for each time step
#     #FireGdf_sfs.save_gdf_trng(tst=tst,ted=ted,fperim=True)
#     FireGdf_sfs_merge.save_gdf_trng(ted=ted,fperim=True)
#     t4 = time.time()
#     print(f'{(t4-t31)/60.} minutes to save large fires')
#
#
#
#     # Run to save year end summary and file lists
#     FireSummary.save_sum_1d(tst,ted)
#     FireSummary.add_heritage(ted)
#     FireSummary.add_largefirelist(ted)
#     t5 = time.time()
#     print(f'{(t5-t4)/60.} minutes to generate summary')
#
#     # Run to clean up large fire geojson files not needed
#     # FireGdf_sfs.yrend_clean(ted)
#
#     t2 = time.time()
#     print(f'{(t2-t1)/60.} minutes to run code')

def testRun():
    import sys
    dircode = '/Users/yangchen/GoogleDrive/My/My.Research/UCI/ProjectsFull/California.fire/Code/fireatlas/'
    sys.path.insert(1, dircode)

    import FireIO, FireMain, FireGpkg,FireGpkg_sfs
    # set the start and end time
    tst=(2019,10,24,'AM')
    ted=(2019,11,2,'PM')

    # example of shape file
    # shp_Cal = FireIO.get_Cal_shp()
    # region = ['California',shp_Cal]

    # example of box region
    ext = [-122.9,38.5,-122.6,38.9]
    region = ['Kincade',ext]

    # example of a country
    # region = ['Brazil','Brazil']

    FireMain.Fire_Forward(tst,ted,restart=False,region=region)
    FireGpkg.save_gdf_trng(tst=tst,ted=ted)
    FireGpkg_sfs.save_gdf_trng(ted)


    # # import libraries
    # from FireLog import logger
    # import FireMain
    # import FireObj
    # import FireIO
    # from FireConsts import dirpjdata
    # import os
    # import glob
    #
    #
    # import time
    # t1 = time.time()
    #
    # restart = False
    # # initialize allfires object (using previous day data or an empty fire object)
    # if FireObj.t_dif(tst,(tst[0],1,1,'AM'))==0:  # force restart at the start of a year
    #     restart = True
    # allfires = FireMain.Fobj_init(tst,restart=restart)
    #
    # # loop over all days during the period
    # endloop = False  # flag to control the ending of the loop
    # t = list(tst)    # t is the time (year,month,day,ampm) for each step
    # while endloop == False:
    #     logger.info('')
    #     logger.info(t)
    #
    #     # 1. record existing active fire ids (for the previous time step)
    #     fids_ea = allfires.fids_active
    #
    #     # 2. update allfires and fire object changes due to temporal progression
    #     # all fires
    #     allfires.update_t(t)           # update t for allfires
    #     allfires.reset_fids_updated()  # reset fids_expanded, fids_new, fids_merged, fids_invalid
    #     # for each fire
    #     allfires.update_t_allfires(t)  # update t
    #     allfires.reset_newpixels()     # reset newpixels
    #
    #     # 3. read active fire pixels from VIIRS dataset
    #     t_read = time.time()
    #     afp = FireIO.read_AFP(t,src='viirs',region=region)
    #     t_read2 = time.time()
    #     logger.info(f'reading file {(t_read2-t_read)}')
    #     print(len(afp))
    #     if len(afp) > 0:
    #         t_expand = time.time()
    #         # 4. do fire expansion/creation using afp
    #         allfires = FireMain.Fire_expand_rtree(allfires,afp,fids_ea)
    #         t_expand2 = time.time()
    #         logger.info(f'expanding fires {(t_expand2-t_expand)}')
    #         # 5. do fire merging using updated fids_ne and fid_ea
    #         fids_ne = allfires.fids_ne                         # new or expanded fires id
    #         fids_ea = sorted(set(fids_ea+allfires.fids_new))   # existing active fires (new fires included)
    #         fids_sleep = allfires.fids_sleeper
    #         t_merge = time.time()
    #         if len(fids_ne) > 0:
    #             allfires = FireMain.Fire_merge_rtree(allfires,fids_ne,fids_ea,fids_sleep)
    #         t_merge2 = time.time()
    #         logger.info(f'merging fires {(t_merge2-t_merge)}')
    #     # 9. loop control
    #     #  - if t reaches ted, set endloop to True to stop the loop
    #     if FireObj.t_dif(t,ted)==0:
    #         endloop = True
    #         # correct fire heritage of last time step
    #         # allfires.heritages = correct_nested_ids(allfires.heritages)
    #     # 10. update t with the next time stamp
    #     t = FireObj.t_nb(t,nb='next')
    # # FireMain.Fire_Forward(tst=tst,ted=ted,region=region)
    #
    # t2 = time.time()
    # print(f'{(t2-t1)/60.} minutes used to run code')


if __name__ == "__main__":
    ''' The main code to run time forwarding for a time period
    '''

    testRun()
    # import FireIO

    # shp_Cal = FireIO.get_Cal_shp()
    # region = ['California',shp_Cal]
    # region = ['Cbox',[-125, -114, 32, 42.5]]
    # cty = 'Brazil'
    # region = [cty,FireIO.get_Cty_shp(cty)]
    # t=(2018,9,1,'AM')
    # aflist = FireIO.read_AFP(t,region=region)
    # print(len(aflist))

    # import os
    # if 'GDAL_DATA' not in os.environ:
    #     os.environ['GDAL_DATA'] = r'C:/Users/rebec/anaconda3/envs/py3work/Library/share/gdal'
    #     os.environ['PROJ_LIB'] = r'C:/Users/rebec/anaconda3/envs/fireatlas/Library/share/proj'

    # Yearbatchrun(2020)

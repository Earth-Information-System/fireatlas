import pytest
import fireatlas
import geopandas as gpd
import os


# make sure that v2 converted data are in tests/data correctly
def test_fireforward_v2_vs_v3(
    tmp_settings_context_manager,
    test_data_dir,
):
    
    # set any other settings needed
    with tmp_settings_context_manager(
        fireatlas.settings,
        LOCAL_PATH=test_data_dir, 
        EPSG_CODE=9311, 
        FIRE_SOURCE="SNPP", 
        FIRE_NRT=False,
        remove_static_sources=False, 
        FTYP_OPT = "global", 
        CONT_OPT="global"
    ):

        from fireatlas import FireMain, postprocess, FireTime, preprocess, settings
    
        tst = [2020, 9, 5, "AM"]
        ted = [2020, 9, 10, "PM"]
        region = ("v3_test_data_for_Creek_SNPP", [-119.5, 36.8, -118.9, 37.7])

        # do preprocessing ahead of time and save in test/data/FEDSpreprocessed
        list_of_ts = list(FireTime.t_generator(tst, ted))

        preprocess.preprocess_region(region, force=True) 

        for t in list_of_ts:
            preprocess.preprocess_region_t(t, region=region, read_location="local", force=True)

        # Run, reading from local. Will call preprocess.read_preprocessed with 
        # read_location="local"
        # In test, local will be overwritten with the test/data directory. 
        allfires, allpixels, t_saved = FireMain.Fire_Forward(tst=tst, ted=ted, restart=True, region=region, read_location="local")

        # test allfires outputs here 

        allpixels = postprocess.read_allpixels(tst, ted, region, location="local")
        allfires_gdf = postprocess.read_allfires_gdf(tst, ted, region, location="local")


        # Save largefire outputs
        large_fires = postprocess.find_largefires(allfires_gdf)
        postprocess.save_large_fires_nplist(allpixels, region, large_fires, tst)
        postprocess.save_large_fires_layers(allfires_gdf, region, large_fires, tst, ted)
        
        
        fireline = gpd.read_file(os.path.join(
            settings.LOCAL_PATH, settings.OUTPUT_DIR, region[0], "2020",
            "Largefire", "1", "fireline.fgb"),
        engine='pyogrio')
        
        newfirepix = gpd.read_file(os.path.join(
            settings.LOCAL_PATH, settings.OUTPUT_DIR, region[0], "2020",
            "Largefire", "1", "newfirepix.fgb"),
        engine='pyogrio')
        
        nfplist = gpd.read_file(os.path.join(
            settings.LOCAL_PATH, settings.OUTPUT_DIR, region[0], "2020",
            "Largefire", "1", "nfplist.fgb"),
        engine='pyogrio')
        
        perimeter = gpd.read_file(os.path.join(
            settings.LOCAL_PATH, settings.OUTPUT_DIR, region[0], "2020",
            "Largefire", "1", "perimeter.fgb"),
        engine='pyogrio')
        
        assert allfires_gdf.crs.equals(fireline.crs)
        assert allfires_gdf.crs.equals(newfirepix.crs)
        assert allfires_gdf.crs.equals(nfplist.crs)
        assert allfires_gdf.crs.equals(perimeter.crs)

        lf_fireline = gpd.read_file(
                os.path.join(
                    settings.LOCAL_PATH, settings.OUTPUT_DIR,region[0], "2020", 
                    "CombinedLargefire", "20201105PM", "lf_fireline.fgb"
                ),
        engine='pyogrio')
        
        lf_newfirepix = gpd.read_file(
                os.path.join(
                    settings.LOCAL_PATH, settings.OUTPUT_DIR,region[0], "2020", 
                    "CombinedLargefire", "20201105PM", "lf_newfirepix.fgb"
                ),
        engine='pyogrio')
        
        lf_perimeter = gpd.read_file(
                os.path.join(
                    settings.LOCAL_PATH, settings.OUTPUT_DIR,region[0], "2020", 
                    "CombinedLargefire", "20201105PM", "lf_perimeter.fgb"
                ),
        engine='pyogrio')
        
        assert allfires_gdf.crs.equals(lf_fireline.crs)
        assert allfires_gdf.crs.equals(lf_newfirepix.crs)
        assert allfires_gdf.crs.equals(lf_perimeter.crs)
        
        

        # assert allfires_gdf.crs.equals(ss_fireline.crs)
        # assert allfires_gdf.crs.equals(ss_newfirepix.crs)
        # assert allfires_gdf.crs.equals(ss_perimeter.crs)
        
        

import pytest
import fireatlas
import geopandas as gpd
import os


@pytest.fixture
def v3_run(tmp_settings_context_manager,
    test_data_dir):
    # run v3 for Creek fire and store outputs locally
    
    with tmp_settings_context_manager(
        fireatlas.settings,
        LOCAL_PATH=test_data_dir, 
        EPSG_CODE=9311, 
        FIRE_SOURCE="SNPP", 
        FIRE_NRT=False,
        remove_static_sources=False, 
        FTYP_OPT = "CA", 
        CONT_OPT="CA", 
        S3_PATH="s3://eis-fire-public" 
    ):

        from fireatlas import FireMain, postprocess, FireTime, preprocess, settings
        from fireatlas.FireLog import logger
        import logging

        logger.setLevel(logging.WARNING) 
    
        tst = [2020, 9, 5, "AM"]
        ted = [2020, 11, 5, "PM"]
        
        region = ("v3_test_data_for_Creek_SNPP", [-119.5, 36.8, -118.9, 37.7])

        # Expects that we have done preprocessing of monthly data ahead of time 
        # and saved in test/data/FEDSpreprocessed. 
        # To reproduce, use FireIO.convert_v2_pkl_to_csv on data saved 
        # in FEDStemp by v2. 
        
        list_of_ts = list(FireTime.t_generator(tst, ted))

        preprocess.preprocess_region(region, force=True) 

        for t in list_of_ts:
            preprocess.preprocess_region_t(
                t, region=region, read_location="local", force=True
            )

        # Run, reading from local. Will call preprocess.read_preprocessed with 
        # read_location="local"
        # In test, local will be overwritten with the test/data directory. 
        allfires, allpixels, t_saved = FireMain.Fire_Forward(
            tst=tst, ted=ted, restart=True, region=region, read_location="local"
        )

        allfires_gdf = postprocess.read_allfires_gdf(tst, ted, region, location="local")

        # save outputs 
        large_fires = postprocess.find_largefires(allfires_gdf)
        postprocess.save_large_fires_nplist(allpixels, region, large_fires, tst)
        postprocess.save_large_fires_layers(allfires_gdf, region, large_fires, tst, ted)
        postprocess.save_snapshots(allfires_gdf, region, tst, ted)

        # read in outputs for testing
        outpath = os.path.join(
            settings.LOCAL_PATH, settings.OUTPUT_DIR,
            region[0], "2020"
        )
        
        fireline = gpd.read_file(os.path.join(
            outpath, "Largefire", "1", "fireline.fgb"),
        engine='pyogrio')
        
        newfirepix = gpd.read_file(os.path.join(
            outpath, "Largefire", "1", "newfirepix.fgb"),
        engine='pyogrio')
        
        nfplist = gpd.read_file(os.path.join(
            outpath, "Largefire", "1", "nfplist.fgb"),
        engine='pyogrio')
        
        perimeter = gpd.read_file(os.path.join(
            outpath, "Largefire", "1", "perimeter.fgb"),
        engine='pyogrio')
        
        lf_fireline = gpd.read_file(
                os.path.join(
                    outpath, "CombinedLargefire", "20201105PM",
                    "lf_fireline.fgb"
                ), 
            engine='pyogrio')
        
        lf_newfirepix = gpd.read_file(
                os.path.join(
                    outpath, 
                    "CombinedLargefire", "20201105PM", "lf_newfirepix.fgb"
                ),
        engine='pyogrio')
        
        lf_perimeter = gpd.read_file(
                os.path.join(
                    outpath, 
                    "CombinedLargefire", "20201105PM", "lf_perimeter.fgb"
                ),
        engine='pyogrio')

        ss_fireline = gpd.read_file(
            os.path.join(
                outpath, "Snapshot", "20200908PM", "fireline.fgb"),
            engine='pyogrio')

        ss_newfirepix = gpd.read_file(
            os.path.join(outpath, "Snapshot", "20200908PM", "fireline.fgb"), 
            engine='pyogrio')

        ss_perimeter = gpd.read_file(
            os.path.join(outpath, "Snapshot", "20200908PM", "perimeter.fgb"), 
            engine='pyogrio')
        
        return (
            allfires_gdf, fireline, newfirepix, nfplist, perimeter, lf_fireline, 
            lf_newfirepix, lf_perimeter, ss_fireline, ss_newfirepix, ss_perimeter
        )

@pytest.fixture
def v2_load(tmp_settings_context_manager, test_data_dir):
    # load v2 from file and return 

      with tmp_settings_context_manager(
        fireatlas.settings,
        LOCAL_PATH=test_data_dir
      ):
        
        from fireatlas import settings
        
        largefire_perim = gpd.read_file(
            os.path.join(
                settings.LOCAL_PATH, "v2_outputs", "Generate_Test_data_for_Creek_CLEANER_SNPP", 
                "2020", "Largefire", "1", "perimeter.fgb" 
            ), 
            engine='fiona'
        )

        return largefire_perim

# this test is slow, so skipped by default. to run: pytest --runslow
@pytest.mark.slow
def test_intersection_over_union(v2_load, v3_run):
    # Ensure that v2 and current version are generating aproximately the same 
    # perimeters, given the same inputs
    AVG_IOU_THRESHOLD = 0.999
    MIN_IOU_THRESHOLD = 0.95
    
    v2 = v2_load.set_index('t')
    v3 = v3_run[4].set_index('t')

    #  v2 CRS may be sensitive to dependency versions/env 
    v2 = v2.to_crs(v3.crs)
    
    ious = []
    
    for t, row in v2.iterrows():
        assert t in v3.index
        
        intersection = row.geometry.intersection(v3.loc[t].geometry)
        union = row.geometry.union(v3.loc[t].geometry)
        iou = intersection.area / union.area
        ious.append(iou)
        
    assert (sum(ious) / float(len(ious)) > AVG_IOU_THRESHOLD)
    assert min(ious) > MIN_IOU_THRESHOLD


@pytest.mark.slow
def test_crs_set_on_fgb_outputs(v3_run):
    
    allfires_gdf = v3_run[0]
    
    for output in v3_run[1:]:
        assert allfires_gdf.crs.equals(output.crs)
       
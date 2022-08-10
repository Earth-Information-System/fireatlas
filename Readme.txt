============= Fire object tracking system using VIIRS active fires =============

[NAME]: Fire object tracking system using VIIRS active fires

[VERSION]: 2.0

[AUTHORS]: Yang Chen (yang.chen@uci.edu)
           Rebecca Scholten
           Stijn Hantson
           Casey Graff
           Shane Coffield

[LANGUAGE]: Python 3.0 with external packages of
  * Numpy (1.17.5)
  * Pandas (1.0.1)
  * Geopandas (0.7.0)
  * Xarray (0.15.0)
  * Scipy (1.4.1)
  * Shapely (1.7.1)
  * Gdal (3.0.4)
  * Pyproj (2.5.0)
  * sklearn

[LIST OF FILES]:
- Data package
  * DataCheckUpdate.py: check and update the input data
- Main package (Fire tracking)
  * FireObj.py:        module containing the object definitions
  * FireMain.py:       main module for running the fire object tracking along time
  * FireConsts.py:     constants and options used for the project
  * FireFuncs.py:      functions used to handle optional calculations
  * FireIO.py:         functions used to read and save data
  * FireTime.py:       functions used to handle times
  * FireVector.py:     functions used for vector related calculations
  * FireClustering.py: functions used for doing fire clustering
  * FireLog.py:        module containing all logging info

- Diagnostic package (output, expandable)
  * FireGdf.py:        module for creating regional geojson summary at each time
                         step
  * FireGdf_sfs.py:    module for creating geojson summary for temporal evolution
                         of each single fire
  * FireGdf_ign.py:    module for creating geojson summary for ignition points
                        of each single fire
  * FireSummary.py:    module for creating a summary of all valid fires up to a
                         time step (usually at year end)

- Post-processing package (for individual large fires)
  * Postprocess.py: lake removing; unburnable area removing; ...
  * Convert2geojson.py: convert geopackage to geojson for better web usage

- Webpage package (use output to create webpage)
  * AT_NRTweb.py: NRT fire webpage creation

- Run package (run fire tracking, diagnostic ouput, or webpage creation)
  * FireRun.py
  * FireRunNRT.py

[COPYRIGHT]: This code is open-source and can be free used for research purposes

[CONTACTS]: yang.chen@uci.edu

=================================== End ========================================

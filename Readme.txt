============= Fire object tracking system using VIIRS active fires =============

[NAME]: Fire object tracking system using VIIRS active fires

[VERSION]: 1.0

[AUTHORS]: Yang Chen (yang.chen@uci.edu), Casey Graff, and Shane Coffield

[LANGUAGE]: Python 3.0 with external packages of
  * Numpy (1.17.5)
  * Pandas (1.0.1)
  * Geopandas (0.7.0)
  * Xarray (0.15.0)
  * Scipy (1.4.1)
  * Shapely (1.8.0)
  * Gdal (3.0.4)
  * Pyproj (2.5.0)
  * sklearn

[LIST OF FILES]:
* FireConsts.py:     constants used for the project
* FireIO.py:         functions used to read and save data
* FireVector.py:     functions used for vector related calculations
* FireClustering.py: functions used for doing fire clustering
* FireLog.py:        module containing all logging info
* FireObj.py:        module containing the object definitions
* FireMain.py:       main module for running the fire object tracking along time
* FireGdf.py:        module for creating regional geojson summary at each time
                       step
* FireGdf_sfs.py:    module for creating geojson summary for temporal evolution
                       of each single fire
* FireSummary.py:    module for creating a summary of all valid fires up to a
                       time step (usually at year end)
* FireRun.py:        module to control different runs

[COPYRIGHT]: This code is open-source and can be free used for research purposes

[CONTACTS]: yang.chen@uci.edu

=================================== End ========================================

import glob
import pandas as pd
import geopandas as gpd
from datetime import datetime
import time
from FireConsts import diroutdata

def load_lf(lf_id,file_path,layer='nfplist',drop_duplicate_geometries=False):
    '''Load in Largefire data'''
    try: gdf = gpd.read_file(file_path,layer=layer)
    except Exception as e: 
        #print(e)
        return
    gdf['ID'] = lf_id
    
    if (drop_duplicate_geometries==True) and (layer != 'nfplist'):
        gdf.sort_values('t',ascending=True,inplace=True)
        gdf = gdf.loc[gdf['geometry'].drop_duplicates(keep='first').index]
    return gdf

def main(year):
    # load in NRT Largefire data
    import s3fs
    s3 = s3fs.S3FileSystem(anon=False)
    s3path = str(f'{diroutdata}CONUS_NRT_DPS/{year}/Largefire/')
    lf_files = [f for f in s3.ls(s3path)]
    
    #lf_files = glob.glob(f'{diroutdata}/CONUS_NRT_DPS/{year}/Largefire/*')
    
    lf_files.sort()
    lf_ids = list(set([file.split('Largefire/')[1].split('_')[0] for file in lf_files])) # unique lf ids

    # each largefire id has a file for each timestep which has entire evolution up to that point.
    # the latest one has the most up-to-date info for that fire
    largefire_dict = dict.fromkeys(lf_ids)
    for lf_id in lf_ids:
        
        most_recent_file = 's3://'+[file for file in lf_files if lf_id in file][-1] # most recent file is last!
        
        largefire_dict[lf_id] = most_recent_file
    
    all_lf_pixels = pd.concat([load_lf(lf_id,file_path,layer='nfplist') for lf_id, file_path in largefire_dict.items()],ignore_index=True)
    all_lf_firelines = pd.concat([load_lf(lf_id,file_path,layer='fireline') for lf_id, file_path in largefire_dict.items()],ignore_index=True)
    all_lf_perimeters = pd.concat([load_lf(lf_id,file_path,layer='perimeter') for lf_id, file_path in largefire_dict.items()],ignore_index=True)

    all_lf_pixels.to_file(f'{diroutdata}CONUS_NRT_DPS/LargeFire_Outputs/{year}/lf_pixels_{year}.geojson')
    all_lf_firelines.to_file(f'{diroutdata}CONUS_NRT_DPS/LargeFire_Outputs/{year}/lf_firelines_{year}.geojson')
    all_lf_perimeters.to_file(f'{diroutdata}CONUS_NRT_DPS/LargeFire_Outputs/{year}/lf_perimeters_{year}.geojson')
    
    return

if __name__ == "__main__":
    
    start = time.time()
    
    main(datetime.now().year)

    total = (time.time() - start)
    
    print('Total runtime is',str(total/60),'minutes')

    
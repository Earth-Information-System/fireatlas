import pickle
from FireObj import *
import glob

def identify_pickle_objects(pickle_file):
    with open(pickle_file, 'rb') as file:
        try:
            obj = pickle.load(file)
            print(f'Successfully loaded object of type: {type(obj)}')
        except (EOFError, Exception) as e:
            print(f'Failed to load object: {e}')

# shared-buckets/gsfc_landslides/FEDSoutput-s3-conus/CONUS_NRT_DPS/2023/Serialization/20230101AM_full.pkl
for item in glob.glob('shared-buckets/gsfc_landslides/FEDSoutput-s3-conus/CONUS_NRT_DPS/2023/Serialization/*.pkl'):
    identify_pickle_objects(item)

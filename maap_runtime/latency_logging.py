import sys
import requests
import pytz
from datetime import datetime, timezone, timedelta


def fire_api_query(base = "https://openveda.cloud/api/features/collections/", collection = "public.eis_fire_lf_fireline_nrt", region = "&region=CONUS"):
  foo = requests.get(f"{base}{collection}/items?f=geojson{region}&sortby=-t")
  if(foo.status_code == 200):
    return(foo.json())
  else:
    foo.raise_for_status()


def get_time_difference_in_data(api_data, overpass_cadence = 12, baseline_latency = 12, some_time_buffer = 0.20, eastern_timezone_region = "US/Eastern"):
  '''
  api_data (dict): json output from an api call to the FEDS api.
  overpass_cadence (float): expected number hours between data collections
  baseline_latency (float): hours of expected latency from overpass to distribution of data through api. 
  some_time_buffer (float): hours added to overpass and latency as a buffer. 
  eastern_timezone_region (str): The eastern-most timezone of the region in question. This must be compatible with pytz.timezone().  
  '''
  now_utc = datetime.now(timezone.utc)
  api_t = api_data['features'][0]['properties']['t']
  api_raw = api_data['features'][0]['properties']['t']
  api_t = datetime.strptime(api_t, "%Y-%m-%dT%H:%M:%S")
  tz = pytz.timezone(eastern_timezone_region) # Converting from local solar to tz at the eastern most timezone of region, which would have the longest latency. 
  api_t = tz.localize(api_t)
  tz_utc = pytz.timezone('UTC') ## Now, we can convert it to UTC, because we want the UTC time of an Eastern overpass. 
  api_t = api_t.astimezone(tz_utc)
  time_diff = now_utc - api_t
  hour_diff = (time_diff.seconds/(60))/60
  if(hour_diff > (overpass_cadence + baseline_latency + some_time_buffer_thresh)):
    # Alert
    print(f"At {now_utc.strftime("%Y-%m-%d %H:%M:%S")} UTC, the API displayed {api_raw}. There were {round(hour_diff, 2)} hours between check time in UTC and the last API data time, or {round(hour_diff - overpass_cadence, 2)} hours since last {eastern_timezone_region} satellite overpass. ")
  
    # make sure calling process gets an bad exit code so it bubbles as failure
    sys.exit(1)




if __name__ == '__main__':
  most_recent_data = fire_api_query()
  get_time_difference_in_data(api_data = most_recent_data)

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


def get_time_difference_in_data(api_data, overpass_cadence = 12, baseline_latency = 10, some_time_buffer = 0, eastern_timezone_region = "'US/Eastern'"):
  '''
  api_data (dict): json output from an api call to the FEDS api.
  overpass_cadence (float): expected number hours between data collections
  baseline_latency (float): hours of expected latency from overpass to distribution of data through api. 
  some_time_buffer (float): hours added to overpass and latency as a buffer. 
  eastern_timezone_region (str): The eastern-most timezone of the region in question. This must be compatible with pytz.timezone().  
  '''
    now_utc = datetime.now(timezone.utc)


    return sorted(filtered_jobs, key=lambda l: l[0])




if __name__ == '__main__':
    jobs = list_jobs()
    failed_jobs = filter_jobs_last_hour(jobs['jobs'], 'job-failed')
    print(f"[ FOUND ]: {len(failed_jobs)} failed jobs")
    for job in failed_jobs:
        job_id, cmd, tag = job
        print(f"{job_id} running command='{cmd}' with tag='{tag}' failed")
    if failed_jobs:
        # make sure calling process gets an bad exit code so it bubbles as failure
        sys.exit(1)

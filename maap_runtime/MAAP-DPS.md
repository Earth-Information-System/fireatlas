# Instructions for running on MAAP DPS

## Setup

- Make sure you can access `https://repo.ops.maap-project.org`, the MAAP instance of GitLab.

## Configuring the algorithm

All algorithm configuration is specified in the `./maap_runtime/<algorithm-dir-name>/algorithm_config.yaml` file.
Make sure that:

- The `version` matches the current Git branch.
- The `repository_url` is an HTTPS URL to the MAAP repository (GitHub will not work for security reasons).
- The `run_command` starts with the name of the repository (i.e., the last part of the repository URL -- probably, `fireatlas`; note, **not** the name of the current directory, the branch, or anything like that).

## Registering the algorithm

1. Do at least a partial interactive test to make sure that your run script will actually run.
  Debugging the MAAP DPS is tricky and time-consuming, so try to identify and squash as many bugs as you can from the ADE.

2. Commit all changes to your workspace **and push them** to the corresponding branch of the **MAAP GitLab repository**.
  The ADE will not let you register an algorithm if you have uncommitted changes.
  However, it's easy to forget to push the changes to MAAP when debugging interactively.

3. If you have previously registered this algorithm under the same name **and version**, delete it.
  This will prevent you from accidentally running an older version of the algorithm.
  From the MAAP ADE, in the top settings bar, go to "DPS/MAS Operations" --> "Delete Algorithm".
  Find your algorithm (**note the version!**) and select it to delete.
    - Alternatively, using the MAAP Python library, you can use code like the following:

      ```python
      from maap.maap import MAAP
      maap = MAAP(maap_host='api.ops.maap-project.org')
      alg = "eis-fire-feds-v2"  # Algorithm name
      version = "conus-dps-4"   # Current branch
      del_result = maap.deleteAlgorithm(f"{alg}_ubuntu:{version}")
      print(del_result.text)  # Confirm that this is code 200 and says "successfully delted job ..."
      ```

4. In the MAAP ADE file browser, right click on the `./maap_runtime/<algorithm-dir-name>/algorithm_config.yaml` file and click "Register as MAS Algorithm".
  In the popup menu, the only option you should change is the MAAP queue; select the smallest node necessary to run your algorithm.
  Submit your selection.

5. You should see a popup window with some JSON information about your algorithm.
  Copy this into a plain text file for reference.
    - The most useful thing in here is the job build log URL, which will look something like: `https://repo.ops.maap-project.org/root/register-job/-/jobs/817/raw`. This URL is used to track the build progress of your algorithm. You will not be able to run the latest version of your algorithm until this build is completed; this usually takes about 10-20 minutes.
    - Note that depending on permissions, you may not actually have access to the build log, in which case you will just have to wait until the algorithm re-appears in the list of algorithms. You can try Python code like the following to monitor the status:

      ```python
      import time
      from maap.maap import MAAP
      maap = MAAP(maap_host="api.ops.maap-project.org")
      alg = "eis-fire-feds-v2"  # Algorithm name
      version = "conus-dps-4"   # Current branch
      while True:
        desc = maap.describeAlgorithm(f"{alg}_ubuntu:{version}")
        if desc.status_code == 200:
          break
        else:
          print(".", end="")    # Primitive 'progress bar', to reassure you that it's still going
          time.sleep(5)  # Sleep for 5 seconds, then try again
      print("Algorithm ready!")
      ```

## Running the algorithm

Only proceed with this step once you are sure the current version of your algorithm is registered and build (see previous steps).

- From the MAAP ADE interface, go to "DPS/MAS Operations" --> "Execute DPS Job".
  Select your algorithm (note the version!) from the drop-down menu and proceed.
  Assuming your algorithm doesn't have any inputs, you should see a popup saying this; proceed.
  You'll see a final popup with the job ID -- **copy this into a plain text file before closing this window**;
  you will need this job ID to check the status of your job, but once this window closes, there is no other way to retrieve the job ID (this is a MAAP bug).

- Alternatively, you can submit jobs using the MAAP Python interface, using code like the following:

      ```python
      import time
      from maap.maap import MAAP
      maap = MAAP(maap_host="api.ops.maap-project.org")
      alg = "eis-fire-feds-v2"  # Algorithm name
      version = "conus-dps-4"   # Current branch
      job = maap.submitJob(
        identifier=f"job-{alg}_ubuntu:{version}",
        username="myuser",  # NOTE: Replace `myuser` with your username
        algo_id=f"{alg}_ubuntu",
        version=version
      )
      print(job)
      ```

You can check on the status of your job using the "DPS/MAS Operations" --> "Get DPS Job Status" menu item, or `maap.getJobStatus` Python function.

Here some code that (1) submits a MAAP DPS job; (2) continuously reports on the job's status while it's running; and (3) prints out the output in a (more) human-readable format once the job is finished.

```python
from maap.maap import MAAP
maap = MAAP(maap_host='api.ops.maap-project.org')
import xml.etree.ElementTree as ET
import time

def getStatus(job):
    jobid = job["job_id"]
    st = ET.fromstring(maap.getJobStatus(jobid).content)
    result = st.find("{http://www.opengis.net/wps/2.0}Status").text
    return result

def getOutput(job):
    jobid = job["job_id"]
    e = ET.fromstring(maap.getJobResult(jobid).content)
    wps = "{http://www.opengis.net/wps/2.0}"
    result = e.find(f"{wps}Output/{wps}Data").text
    return result

def sleepJob(job, sleeptime=5):
    while not getStatus(job) in ("Failed", "Succeeded"):
        current_status = getStatus(job)
        print(f'Status: {current_status}')
        while getStatus(job) == current_status:
            print('.', end='')
            time.sleep(sleeptime)
        print("\n", end='')
    print(f'Status: {getStatus(job)}')

alg = "eis-fire-feds-v2"  # Algorithm name
version = "conus-dps-4"   # Current branch

submit = maap.submitJob(
    identifier=f"job-{alg}_ubuntu:{version}",
    username="ashiklom",
    algo_id=f"{alg}_ubuntu",
    version=version
)
print(submit)

# Continuously monitor job status until the job succeeds or fails
sleepJob(submit)

# Print the output of the job once it finishes
print(getOutput(job))
```

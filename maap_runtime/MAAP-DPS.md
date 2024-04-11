# Async Task Runner

## Overview:

Currently the EIS team uses the MAAP-DPS async executor (talked about below) to run their tasks. However, for reasons
related to portability, data volume and simplistic idioms we've decided to leverage Dask for parallelism instead of DPS. 
What that means in practice is that you can run this algorithm on any cloud compute instance (EC2, GCP) 
or cluster (k8s, JupyterHub) taking advantage of the number of cores as long as it has Dask installed (which is a requirement in
`pyproject.toml`). In the future we might tweaks things slightly so you can run against remote Dask clusters as well.

---

# MAAP DPS

## DPS Docs:

https://docs.maap-project.org/en/latest/technical_tutorials/dps_tutorial/dps_tutorial_demo.html#Data-Processing-System-(DPS)-Tutorial-A-to-Z

## Registering the Algorithm for DPS Runs:

1. Make any necessary code changes on your feature branch or main branch


2. Decide on a tag name that is unique (it can be a semvar version or just text). For purposes of this documentation we'll call the new tag `testing123`


3. Change the `algorithm_version: testing123` line in `/maap_runtime/<algo-name>/algorithm_config.yaml` files to point to this new tag 


4. Git commit all your changes and push remotely

  ```bash
  git add -u
  git commit -m "updates with new tag"
  git push origin <branch-name>
  ```

5. Then tag and push your tags remotely

  ```bash
  git tag -f testing123
  git push origin -f testing123
  ```

6. Then register your algorithm. Really the output of this step is just building a new DPS image. Either you can manually do this as talked about in the DPS docs above or automate it like in `/maap_runtim/register-all.ipynb`


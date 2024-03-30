# Getting Started

### Install `fireatlas` as a Package

We use `pyproject.toml` with `pip` to "deploy" fireatlas locally inside your environment:

```bash
$ git clone <this-repo> fireatlas

$ cd fireatlas/

# note that '-e' allows you to edit your code and all changes will be available without reinstall
$ pip install -e .

# or if you need testing infrastructure too 
$ pip install -e '.[test]'
```

### Running the Algorithm

There are two different environments where you can run all algorithm steps: JupyterHub and MAAP-DPS


#### JupyterHub

1. If you have accesss to JupyterHub and have installed the conda environment from `env.yml` then you can walk through the notebooks interactively in `/notebooks` to run each step

2. Or if you'd rather just run one script that coordinates all the work for you in the right order with dask parallelism then you can run `fireatlas/FireRunDaskCoordinator.py`



#### DPS

DPS is an async task runner. For more info [look over here](../maap_runtime/MAAP-DPS.md)


1. You can run individual steps asynchronously on DPS via the notebook at `/maap_runtime/dps-coordinator.ipynb`


2. Or you can kick off all steps in a single DPS job via the last step in the notebook at `/maap_runtime/dps-coordinator.ipynb`






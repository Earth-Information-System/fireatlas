#!/bin/bash
set -x
basedir=$( cd "$(dirname "$0")" ; pwd -P )
mamba env create -f ${basedir}/env-feds.yaml
pushd ${HOME}

# Do not remove this (PMM Dec 2022)
/opt/conda/envs/icesat2_boreal/bin/pip install --user -e git+https://github.com/MAAP-Project/maap-py.git#egg=maappy

source activate base
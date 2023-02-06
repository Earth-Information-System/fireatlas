#!/bin/bash
set -euxo pipefail
basedir=$( cd "$(dirname "$0")" ; pwd -P )
mamba env create -f ${basedir}/env-feds.yaml
pushd ${HOME}

# required to get all the MAAP(y) things
/opt/conda/envs/feds/bin/pip install --user -e git+https://github.com/MAAP-Project/maap-py.git#egg=maappy

conda info --envs
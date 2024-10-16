#!/bin/bash
set -eo pipefail
basedir=$( cd "$(dirname "$0")"; pwd -P )
echo "Basedir: $basedir"
echo "Initial working directory: $(pwd -P)"
echo "conda: $(which conda)"
echo "Python: $(which python)"
python --version
# the vanilla image uses conda version 23.10
# where mamba should be default resolver
pushd "$basedir"

conda create -n "fire_env" python=3.11
source activate fire_env

echo "Doing fireatlas install next..."
/opt/conda/envs/fire_env/bin/pip install -e ..
/opt/conda/envs/fire_env/bin/pip install "git+https://github.com/MAAP-Project/maap-py.git@develop"


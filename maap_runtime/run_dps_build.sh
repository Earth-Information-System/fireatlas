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

conda create -n "fire_env"
source activate fire_env
conda install pip 
which -a pip

echo "Doing fireatlas install next..."
/opt/conda/envs/fire_env/bin/pip install -e ..
/opt/conda/envs/fire_env/bin/pip install "git+https://github.com/MAAP-Project/maap-py.git@develop"


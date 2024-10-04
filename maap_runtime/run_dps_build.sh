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
source activate /opt/conda/envs/python

echo "Doing fireatlas install next..."
/opt/conda/envs/python/bin/pip install -e ..
/opt/conda/envs/python/bin/pip install "git+https://github.com/MAAP-Project/maap-py.git@develop"


#!/bin/bash
# set -euxo pipefail
set -eo pipefail
basedir=$( cd "$(dirname "$0")"; pwd -P )
echo "Basedir: $basedir"
echo "Initial working directory: $(pwd -P)"
echo "conda: $(which conda)"
echo "Python: $(which python)"
python --version
# the vanilla image uses conda version 23.10
# where mamba should be default resolver
conda env create -f env-feds-2024.yml

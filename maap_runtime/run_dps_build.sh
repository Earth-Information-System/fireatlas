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
source activate /opt/conda/envs/vanilla
# vanilla:3.1.5 makes all our code just hang for reasons know one knows about
# so we explicitly upgrade just s3fs here out-of-band of pyproject.toml and env.yml
conda update -c conda-forge s3fs -y > update_log.txt 2>&1
cat update_log.txt
conda list | grep s3fs
echo "Doing fireatlas install next..."
/opt/conda/envs/vanilla/bin/pip install -e ..
/opt/conda/envs/vanilla/bin/pip install "git+https://github.com/MAAP-Project/maap-py.git@develop"


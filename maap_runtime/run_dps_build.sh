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

# Build from the pinned lock, NOT the loose env.yml
conda env create -f ../env.lock.yml
source activate fire_env

# The lock pins every dependency but deliberately omits fireatlas itself, so
# install the package from the repo checked out in this image. --no-deps
# because the lock already provides dependencies. -e (editable) so that
# root_dir in FireConsts.py still resolves to this checkout instead of a
# copy under site-packages, which run_dps_cli.sh's log-copy step relies on.
echo "Installing fireatlas..."
pip install --no-deps -e "$basedir/.."

# Fail loudly rather than discovering an inconsistent environment at
# runtime on DPS.
pip check

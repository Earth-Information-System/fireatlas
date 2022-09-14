#!/bin/bash
# set -euxo pipefail
set -eo pipefail

basedir=$( cd "$(dirname "$0")"; pwd -P )
echo "Basedir: $basedir"
echo "Initial working directory: $(pwd -P)"

echo "conda: $(which conda)"

# Trying to resolve conda permissiosn issue
# https://gitter.im/conda/conda?at=5dc427aa2f8a034357513172
export CONDA_PKGS_DIRS="$basedir/.conda"
mkdir -p "$CONDA_PKGS_DIRS"

conda env create -f "$basedir/env-feds.yml" -p "$basedir/env-feds"
source activate "$basedir/env-feds"

echo "Python: $(which python)"
python --version
echo "pip: $(which pip)"
pip --version

echo "Starting algorithm in subshell"
(
cd "$basedir"
echo "Switched to directory: $(pwd -P)"
python -u c "import FireRun; FireRun.CArun()"
)
echo "Done!"

exit

#!/bin/bash
# set -euxo pipefail
set -eo pipefail

basedir=$( cd "$(dirname "$0")"; pwd -P )
echo "Basedir: $basedir"
echo "Initial working directory: $(pwd -P)"

echo "conda: $(which conda)"

# Only create conda environment if it's not currently active. This allows for
# interactive testing of this script.
if [[ $CONDA_PREFIX != "/projects/env-feds" ]]; then
  # Trying to resolve conda permissiosn issue
  # https://gitter.im/conda/conda?at=5dc427aa2f8a034357513172
  export CONDA_PKGS_DIRS="/projects/.conda"
  mkdir -p "$CONDA_PKGS_DIRS"
  conda env create -f "$basedir/env-feds.yml" -p "/projects/env-feds"
  source activate "/projects/env-feds"
fi

echo "Python: $(which python)"
python --version

echo "Starting algorithm in subshell"
(
cd "$basedir"
echo "Running in directory: $(pwd -P)"
# python FireRun.py $1
python -u -c "import FireRun; FireRun.CreekSamplerun()"
# python -u -c "import FireRun; FireRun.CArun()"
)
echo "Done!"

exit

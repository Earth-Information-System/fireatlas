#!/bin/bash
# set -euxo pipefail
set -eo pipefail
mkdir output
basedir=$( cd "$(dirname "$0")"; pwd -P )
echo "Basedir: $basedir"
echo "Initial working directory: $(pwd -P)"
echo "conda: $(which conda)"
echo "Python: $(which python)"
python --version
source activate /opt/conda/envs/env-feds-dask
conda install --solver=libmamba -c conda-forge scalene
echo "Starting algorithm in subshell"
(
pushd "$basedir"
{ # try
  echo "Running in directory: $(pwd -P)"
  scalene --cli --no-browser --reduced-profile --html --column-width 180 --outfile ${PWD}/output/profile.html --- FireRun.py $1
  popd
  echo "Copying log to special output dir"
  cp "$basedir/running.log" ./output
  #cp "$basedir/dask-report.html" ./output

} || { # catch
  popd
  echo "Copying log to special output dir"
  cp "$basedir/running.log" ./output
  #cp "$basedir/dask-report.html" ./output
}
)
echo "Done!"

exit

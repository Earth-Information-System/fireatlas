#!/bin/bash
# set -euxo pipefail
set -xeo pipefail
output_dir=${PWD}/output
mkdir "${output_dir}"
basedir=$( cd "$(dirname "$0")"; pwd -P )
echo "Basedir: $basedir"
echo "Initial working directory: $(pwd -P)"
echo "conda: $(which conda)"
echo "Python: $(which python)"
python --version
conda_prefix=/opt/conda/envs/env-feds-dask
source activate "${conda_prefix}"
conda install -c conda-forge scalene -y -p "${conda_prefix}"
echo "Starting algorithm in subshell"
(
pushd "$basedir"
{ # try
  echo "Running in directory: $(pwd -P)"
scalene --cli --no-browser --reduced-profile --html --column-width 180 --outfile "${output_dir}/profile.html" --- FireRun.py $1
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

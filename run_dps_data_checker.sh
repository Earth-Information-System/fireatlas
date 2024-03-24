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
pip install -r env-feds-2024-dps-requirements.txt
echo "Starting algorithm in subshell"
(
pushd "$basedir"
{ # try
  echo "Running in directory: $(pwd -P)"
  scalene --cli --no-browser --reduced-profile --html --column-width 180 --outfile "${output_dir}/profile.html" --- FireRunDataUpdateChecker.py
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

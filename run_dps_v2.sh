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
echo "Starting algorithm in subshell"
(
pushd "$basedir"
{ # try
  echo "Running in directory: $(pwd -P)"
  export_flag="$5"
  if [[ -n "$export_flag" && ( "$export_flag" == "true" ) ]]; then
    # these need to be double quoted to pass along strings correctly formatted
    python FireRunByRegion.py --regnm=$1 --bbox="$2" --tst="$3" --ted="$4" --export
  elif [[ -n "$export_flag" && ( "$export_flag" == "false" ) ]]; then
    python FireRunByRegion.py --regnm=$1 --bbox="$2" --tst="$3" --ted="$4"
  else
    # no export flag passed
    python FireRunByRegion.py --regnm=$1 --bbox="$2" --tst="$3" --ted="$4"
  fi
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

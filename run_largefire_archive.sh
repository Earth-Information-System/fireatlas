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
echo "Starting algorithm in subshell"
(
pushd "$basedir"
{ # try
  echo "Running in directory: $(pwd -P)"
  python combine_largefire.py -s $2 -e $3 -p -x --folder-name $1
  popd
  echo "Copying log to special output dir"
  cp "$basedir/running.log" ./output
} || { # catch
  popd
  echo "Copying log to special output dir"
  cp "$basedir/running.log" ./output
}
)
echo "Done!"

exit

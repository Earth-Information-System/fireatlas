#!/bin/bash
set -eo pipefail
mkdir output
wrkdir=$(pwd -P)
basedir=$( cd "$(dirname "$0")"; pwd -P )
echo "Basedir: $basedir"
echo "Initial working directory: $(pwd -P)"
echo "conda: $(which conda)"
echo "Python: $(which python)"
python --version
source activate /opt/conda/envs/vanilla
echo "Starting algorithm in subshell"
(
pushd "$basedir"
{ # try
  echo "Running in directory: $(pwd -P)"
  # python3 FireRunByRegionAndT.py --regnm="CaliTestRun" --bbox="[-125,36,-117,42]"
  scalene --cli --no-browser --reduced-profile --html --column-width 180 \
      --outfile "${wrkdir}/output/profile.html" --- FireRunByRegionAndT.py --regnm=$1 --t="$2"
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

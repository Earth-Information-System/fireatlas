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
source activate /opt/conda/envs/fire
echo "Starting algorithm in subshell"
(
pushd "$basedir"
{ # try
  echo "Running in directory: $(pwd -P)"
  # python3 FireRunFireForward.py --regnm="CaliTestRun"  --tst="[2023,6,1,\"AM\"]" --ted="[2023,9,1,\"AM\"]"
  scalene --cli --no-browser --reduced-profile --html --column-width 180 \
      --outfile "${output}/profile.html" --- FireRunFireForward.py --regnm=$1 --tst="$3" --ted="$4"
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

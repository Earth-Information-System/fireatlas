#!/bin/bash
set -eo pipefail
output_dir=${PWD}/output
mkdir "${output_dir}"
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
  CURRENT_YEAR=$(date +'%Y')
  scalene --cli --no-browser --reduced-profile --html --column-width 180 \
    --outfile "${output_dir}/profile.html" --- ../fireatlas/combine_largefire.py -s $CURRENT_YEAR -e $CURRENT_YEAR -x --nrt --folder-name $1
  popd
  echo "Copying log to special output dir"
  cp "$basedir/running.log" "$output_dir"
} || { # catch
  popd
  echo "Copying log to special output dir"
  cp "$basedir/running.log" "$output_dir"
}
)
echo "Done!"

exit

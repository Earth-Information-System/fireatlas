#!/bin/bash
# set -euxo pipefail
set -eo pipefail

# base dir set
basedir=$( cd "$(dirname "$0")"; pwd -P )
echo "Basedir: $basedir"
echo "Initial working directory: $(pwd -P)"

# install micromamba
curl micro.mamba.pm/install.sh | bash
echo $MAMBA_EXE
# set variables
export MAMBA_EXE="$basedir/micromamba" # path to micromamba binary
export MY_MAMBA_ENV="$basedir/micromamba/envs/rio-tiler-new"
echo $PATH
echo "micromamba version: $($MAMBA_EXE --version)"

# flag - may need to modify order of arguments
"$MAMBA_EXE" create -f "$basedir/rio-tiler-new.yml" -p "$MY_MAMBA_ENV"

echo "Starting algorithm in subshell"
(
pushd "$basedir"
{ # try
  echo "Running in directory: $(pwd -P)"
  "$MAMBA_EXE" run -p "$MY_MAMBA_ENV" python -u -c "import FireRun; FireRun.CreekSamplerun()"
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

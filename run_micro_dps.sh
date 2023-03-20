#!/bin/bash
# set -euxo pipefail
set -eo pipefail

# base dir set
basedir=$( cd "$(dirname "$0")"; pwd -P )
echo "Basedir: $basedir"
echo "Initial working directory: $(pwd -P)"

# automatic install micromamba
# curl micro.mamba.pm/install.sh | bash

# micromamba binary dir
MICROMAMBA_BIN_DIR="$basedir/bin"
# create that dir if it doesn't exist
# mkdir "$MICROMAMBA_BIN_DIR"
# set mambda micro exe 
MICROMAMBA_EXE=$(pwd -P)/bin/micromamba
echo $MICROMAMBA_EXE
# Download and extract binary to that location
curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj bin/micromamba

MY_MAMBA_ENV="$basedir/micromamba/envs/rio-tiler-new"
# MY_MAMBA_ENV="$MICROMAMBA_BIN_DIR/envs/rio-tiler-new"

echo $PATH
echo "micromamba version: $($MICROMAMBA_EXE --version)"

# flag - may need to modify order of arguments
"$MICROMAMBA_EXE" create -f "$basedir/rio-tiler-new.yml" -p "$MY_MAMBA_ENV"

echo "Starting algorithm in subshell"
(
pushd "$basedir"
{ # try
  echo "Running in directory: $(pwd -P)"
  "$MICROMAMBA_EXE" run -p "$MY_MAMBA_ENV" python -u -c "import FireRun; FireRun.CreekSamplerun()"
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

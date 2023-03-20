#!/bin/bash
# set -euxo pipefail
set -eo pipefail

# base dir set
basedir=$( cd "$(dirname "$0")"; pwd -P )
echo "Basedir: $basedir"
echo "Initial working directory: $(pwd -P)"

# install micromamba
curl micro.mamba.pm/install.sh | bash
# set variables
MAMBA_EXE="$basedir/micromamba" # path to micromamba binary
MY_MAMBA_ENV="$basedir/micromamba/envs/rio-tiler-new"

# flag - may need to modify order of arguments
"$MAMBA_EXE" create -p "$MY_MAMBA_ENV" "$basedir/rio-tiler-new.yml"

# creek sample run with exec
"$MAMBA_EXE" run -p "$MY_MAMBA_ENV" python -u -c "import FireRun; FireRun.CreekSamplerun()"

echo "Done!"

exit
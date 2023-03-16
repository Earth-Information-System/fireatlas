#!/bin/bash
# set -euxo pipefail
set -eo pipefail

basedir=$( cd "$(dirname "$0")"; pwd -P )
echo "Basedir: $basedir"
echo "Initial working directory: $(pwd -P)"

# install micromamba
curl micro.mamba.pm/install.sh | bash
export MAMBA_ROOT_PREFIX="/projects/micromamba"  # optional, defaults to ~/micromamba
eval "$(/projects/.local/bin/micromamba shell hook -s posix)"
echo "micromamba version: $(micromamba --version)"
echo "conda: $(which conda)"


[ -d "/projects/micromamba/envs" ] && echo "basic mamba envs ready: /projects/micromamba/envs"
[ -d "/projects/micromamba/envs/rio-tiler-new" ] && echo "rio-tiler-new already exists at path: /projects/micromamba/envs/rio-tiler-new"


if [ ! -d "/projects/micromamba/envs/rio-tiler-new" ]; then
    micromamba create -f "/projects/rio-tiler-new.yml"
fi

# conda prefix check
# FLAG - CONTRADICTION: conda claiming path DNE...

echo $CONDA_PREFIX
cd -
echo "$PWD"

if [ $CONDA_PREFIX != "/projects/micromamba/envs/rio-tiler-new" ]; then
    # activate env
    echo "activating rio-tiler-new..."
    micromamba activate rio-tiler-new
fi 

echo "$PWD"
cd -
cd -
[$CONDA_PREFIX != "/projects/micromamba/envs/rio-tiler-new"] && echo "Env activ failed"
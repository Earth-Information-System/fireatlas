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
    micromamba create -f "/projects/rio-tiler-new.yml" -p "/projects/micromamba/envs/rio-tiler-new"
fi

# conda prefix check
echo "CONDA_PREFIX: $CONDA_PREFIX"
echo "$PWD"
echo "Python: $(which python)"
python --version

micromamba run -p /projects/micromamba/envs/rio-tiler-new python -u -c "import FireRun; FireRun.CreekSamplerun()"
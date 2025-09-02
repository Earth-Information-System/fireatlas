#!/bin/bash
set -eo pipefail 

copy_s3_object() {
    local from_path="$1"
    local to_path="$2"
    if ! aws s3 cp "$from_path" "$to_path" >/dev/null 2>&1; then
        # log the error quietly, do not stop the script if fails
        echo "Copy failed from $from_path to $to_path, continuing..." >&2
    else
        echo "Copy succeeded from $from_path to $to_path"
    fi
}


run_id=$1

echo "Running script with run_id: $run_id"

output_dir=${PWD}/output
mkdir "${output_dir}"

basedir=$( cd "$(dirname "$0")"; pwd -P )
echo "Basedir: $basedir"
echo "Initial working directory: $(pwd -P)"
echo "conda: $(which conda)"
echo "Python: $(which python)"

python --version
source activate fire_env
conda list | grep s3fs

handle_exit() {
  popd
  echo "Copying log to special output dir"
  cp "$basedir/../running.log" "$output_dir"
  # force the calling process to know we've encountered an error put this in DPS failed state
  exit 128
}

trap 'handle_exit' EXIT

pushd "$basedir"
echo "Running in directory: $(pwd -P)"
# we now secretly look for s3://maap-ops-workspace/shared/gsfc_landslides/FEDSpreprocessed/<regnm>/.env
# and copy it locally to ../fireatlas/.env so that pydantic can pick up our overrides
copy_s3_object "s3://maap-ops-workspace/shared/gsfc_landslides/FEDSpreprocessed/${run_id}/.env" ../fireatlas/.env
copy_s3_object "s3://maap-ops-workspace/shared/zbecker/FEDSstaging/FEDSinput/run_definitions/${run_id}/run_config.yaml" ../fireatlas/run_config.yaml
ls -lah ../fireatlas/

python FireRunArchiveCoordinator.py "$run_id"

popd
echo "Copying log to special output dir"
cp "$basedir/../running.log" "$output_dir"

# unset trap since we are successful and send exit
trap - EXIT
echo "Done!"
exit 0






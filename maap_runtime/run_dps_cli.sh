#!/bin/bash
set -eo pipefail

# force all datetimes to be UTC
export TZ="Etc/UTC"

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

# check if the number of arguments is less than 5 (minimum required)
if [ $# -lt 5 ]; then
    echo "Usage: $0 regnm bbox tst ted --flag"
    echo ""
    echo "Usage: all arguments are passed positionally..."
    echo ""
    echo "Flags: --data-update, --preprocess-region, --preprocess-region-t, --fire-forward, --coordinate-all, --coordinate-all-no-veda-copy"
    echo ""
    echo 'Example: bash run_dps_cli.sh ShastaTrinity [-126,57,-116,37] [2023,1,1,"AM"] [2023,12,31,"PM"] --data-update'
    exit 1
fi

# required arguments
regnm=$1
bbox=$2
tst=$3
ted=$4

# shift the first four arguments to access the flags
shift 4

selected_flag=""

# Parse the rest of the arguments as flags
for arg in "$@"
do
    case $arg in
        --data-update)
            if [[ -n $selected_flag ]]; then
                echo "Error: More than one flag specified"
                exit 1
            fi
            selected_flag="data-update"
            ;;
        --preprocess-region)
            if [[ -n $selected_flag ]]; then
                echo "Error: More than one flag specified"
                exit 1
            fi
            selected_flag="preprocess-region"
            ;;
        --preprocess-region-t)
            if [[ -n $selected_flag ]]; then
                echo "Error: More than one flag specified"
                exit 1
            fi
            selected_flag="preprocess-region-t"
            ;;
        --fire-forward)
            if [[ -n $selected_flag ]]; then
                echo "Error: More than one flag specified"
                exit 1
            fi
            selected_flag="fire-forward"
            ;;
        --coordinate-all)
            if [[ -n $selected_flag ]]; then
                echo "Error: More than one flag specified"
                exit 1
            fi
            selected_flag="coordinate-all"
            ;;
        --coordinate-all-no-veda-copy)
            if [[ -n $selected_flag ]]; then
                echo "Error: More than one flag specified"
                exit 1
            fi
            selected_flag="coordinate-all-no-veda-copy"
            ;;
        *)
            echo "Unknown flag: $arg"
            exit 1
            ;;
    esac
done

# Check if a flag was selected
if [[ -z $selected_flag ]]; then
    echo "Error: No flag specified"
    exit 1
fi


echo "Running script with:"
echo "regnm: $regnm"
echo "bbox: $bbox"
echo "tst: $tst"
echo "ted: $ted"
echo "flag: $selected_flag"

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
copy_s3_object "s3://maap-ops-workspace/shared/gsfc_landslides/FEDSpreprocessed/${regnm}/.env" ../fireatlas/.env
ls -lah ../fireatlas/

if [[ $selected_flag == "data-update" ]]; then
  #scalene --cli --no-browser --reduced-profile --html --column-width 180 --outfile "${output_dir}/profile.html" --- FireRunDataUpdateChecker.py
  python FireRunDataUpdateChecker.py
elif [[ $selected_flag == "preprocess-region" ]]; then
  #scalene --cli --no-browser --reduced-profile --html --column-width 180 \
  #  --outfile "${output_dir}/profile.html" --- FireRunPreprocessRegion.py --regnm=$regnm --bbox="$bbox"
  python FireRunPreprocessRegion.py --regnm=$regnm --bbox="$bbox"
elif [[ $selected_flag == "preprocess-region-t" ]]; then
  #scalene --cli --no-browser --reduced-profile --html --column-width 180 \
  #  --outfile "${output_dir}/profile.html" --- FireRunByRegionAndT.py --regnm=$regnm --tst="$tst" --ted="$ted"
  python FireRunByRegionAndT.py --regnm=$regnm --tst="$tst" --ted="$ted"
elif [[ $selected_flag == "fire-forward" ]]; then
  #scalene --cli --no-browser --reduced-profile --html --column-width 180 \
  #  --outfile "${output_dir}/profile.html" --- FireRunFireForward.py --regnm=$regnm --tst="$tst" --ted="$ted"
  python FireRunFireForward.py --regnm=$regnm --tst="$tst" --ted="$ted"
elif [[ $selected_flag == "coordinate-all" ]]; then
  #scalene --cli --no-browser --reduced-profile --html --column-width 180 \
  #  --outfile "${output_dir}/profile.html" --- ../fireatlas/FireRunDaskCoordinator.py --regnm=$regnm --bbox="$bbox" --tst="$tst" --ted="$ted"
  python ../fireatlas/FireRunDaskCoordinator.py --regnm=$regnm --bbox="$bbox" --tst="$tst" --ted="$ted"
elif [[ $selected_flag == "coordinate-all-no-veda-copy" ]]; then
  #scalene --cli --no-browser --reduced-profile --html --column-width 180 \
  #  --outfile "${output_dir}/profile.html" --- ../fireatlas/FireRunDaskCoordinator.py --regnm=$regnm --bbox="$bbox" --tst="$tst" --ted="$ted"
  python ../fireatlas/FireRunDaskCoordinator.py --regnm=$regnm --bbox="$bbox" --tst="$tst" --ted="$ted" --no-veda-copy
fi

popd
echo "Copying log to special output dir"
cp "$basedir/../running.log" "$output_dir"

# unset trap since we are successful and send exit
trap - EXIT
echo "Done!"
exit 0


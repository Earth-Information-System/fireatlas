#!/bin/bash
set -eo pipefail

# check if the number of arguments is less than 4 (minimum required)
if [ $# -lt 5 ]; then
    echo "Usage: $0 regnm bbox tst ted --flag"
    echo ""
    echo "Usage: all arguments are passed positionally..."
    echo ""
    echo "Flags: --data-update, --preprocess-region, --preprocess-region-t, --fire-forward, --coordinate-all"
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

# Your script logic here
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
source activate /opt/conda/envs/vanilla

(
pushd "$basedir"
{ # try
  echo "Running in directory: $(pwd -P)"

  if [[ $selected_flag == "data-update" ]]; then
    scalene --cli --no-browser --reduced-profile --html --column-width 180 --outfile "${output_dir}/profile.html" --- FireRunDataUpdateChecker.py
  elif [[ $selected_flag == "preprocess-region" ]]; then
    scalene --cli --no-browser --reduced-profile --html --column-width 180 \
      --outfile "${output_dir}/profile.html" --- FireRunPreprocessRegion.py --regnm=$regnm --bbox="$bbox"
  elif [[ $selected_flag == "preprocess-region-t" ]]; then
    scalene --cli --no-browser --reduced-profile --html --column-width 180 \
      --outfile "${output_dir}/profile.html" --- FireRunByRegionAndT.py --regnm=$regnm --t="$tst"
  elif [[ $selected_flag == "fire-forward" ]]; then
    scalene --cli --no-browser --reduced-profile --html --column-width 180 \
      --outfile "${output_dir}/profile.html" --- FireRunFireForward.py --regnm=$regnm --tst="$tst" --ted="$ted"
  elif [[ $selected_flag == "coordinate-all" ]]; then
    scalene --cli --no-browser --reduced-profile --html --column-width 180 \
      --outfile "${output_dir}/profile.html" --- ../fireatlas/FireRunDaskCoordinator.py --regnm=$regnm --bbox="$bbox" --tst="$tst" --ted="$ted"
  fi

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

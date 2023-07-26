TARGETS=$(cat << EOF
./combine_largefire_archive/algorithm_config.yaml
./combine_largefire_nrt/algorithm_config.yaml
./conus_nrt/algorithm_config.yaml
./date_update_checker/algorithm_config.yaml
EOF)

for target in $TARGETS; do
  python3 register_algorithm.py $target
done
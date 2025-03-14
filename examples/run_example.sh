#!/bin/bash
tumour_id=LTXSIM001
script_dir=$(dirname "$(realpath "$0")")
input_directory="${script_dir}/example_cohort/input"
output_directory="${script_dir}/example_cohort/output/${tumour_id}"
segments_directory="${input_directory}/${tumour_id}/segments"
input_files=$(ls "${segments_directory}" | tr "\n" " ")

python3 "${script_dir}/../alpaca/run_alpaca.py" \
    --input_data_directory "${input_directory}" \
    --input_files ${input_files} \
    --output_directory "${output_directory}" \
    --ci_table_name ''

# Merge all the output files into a single file
awk 'FNR==1 && NR!=1 { next; } { print; }' $(ls "${output_directory}"/*.csv) > "${output_directory}/final_${tumour_id}.csv"

# Remove segment files
find "${output_directory}" -type f -name "*.csv" ! -name "final*" -exec rm {} +
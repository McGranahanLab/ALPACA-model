#!/bin/bash
tumour_id=LTX0530
input_tumour_directory="example_cohort/input/${tumour_id}"
output_directory="example_cohort/output/${tumour_id}"

echo "Tumour ID: ${tumour_id}"
# run alpaca:
alpaca \
    --input_tumour_directory "${input_tumour_directory}" \
    --output_directory "${output_directory}"

# get cn change to ancestor:
get_cn_change_to_ancestor \
    --output_path "${output_directory}" \
    --tumour_df_path "${output_directory}/final_${tumour_id}.csv" \
    --tree_path "${input_tumour_directory}/tree_paths.json"
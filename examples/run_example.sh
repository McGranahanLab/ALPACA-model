#!/bin/bash
tumour_id=LTX0000-Tumour1
input_tumour_directory="examples/example_cohort/input/${tumour_id}"
output_directory="examples/example_cohort/output/${tumour_id}"

echo "Tumour ID: ${tumour_id}"

# convert CONIPHER and Refphase outputs to ALPACA input:

refphase_rData="${input_tumour_directory}/${tumour_id}-refphase-results.RData"
CONIPHER_tree_object="${input_tumour_directory}/${tumour_id}.tree.RDS"
output_dir="${input_tumour_directory}"

alpaca input-conversion \
 --tumour_id $tumour_id \
 --refphase_rData $refphase_rData \
 --CONIPHER_tree_object $CONIPHER_tree_object \
 --output_dir $output_dir

# run alpaca:
alpaca run \
    --input_tumour_directory "${input_tumour_directory}" \
    --output_directory "${output_directory}"

# get cn change to ancestor:
alpaca ancestor-delta \
    --output_directory "${output_directory}" \
    --tumour_df_path "${output_directory}/ALPACA_output_${tumour_id}.csv" \
    --tree_path "${input_tumour_directory}/tree_paths.json"

# calculate clone copy number diversity:
alpaca ccd \
    --output_directory "${output_directory}" \
    --alpaca_output_path "${output_directory}/ALPACA_output_${tumour_id}.csv"

#!/bin/bash
tumour_id=LTX0000-Tumour1
input_tumour_directory="tests/x_chromosome/input/${tumour_id}"
output_directory="tests/x_chromosome/output/${tumour_id}"


refphase_rData="${input_tumour_directory}/${tumour_id}-refphase-results.RData"
CONIPHER_tree_object="${input_tumour_directory}/${tumour_id}.tree.RDS"
conversion_output_dir="${input_tumour_directory}"

alpaca input-conversion \
 --tumour_id $tumour_id \
 --refphase_rData $refphase_rData \
 --CONIPHER_tree_object $CONIPHER_tree_object \
 --output_dir $conversion_output_dir


# run alpaca:
alpaca run \
    --input_tumour_directory "${input_tumour_directory}" \
    --output_directory "${output_directory}" \
    --overwrite_output 0 \
    --solver pyomo \
    --pyomo_solver glpk \
    --debug
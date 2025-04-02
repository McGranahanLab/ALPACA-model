#!/bin/bash
tumour_id=LTX0000-Tumour1
output_dir="tests/test_output"
CONIPHER_tree_object="examples/example_cohort/input/${tumour_id}/${tumour_id}.tree.RDS"
SCRIPT_DIR="alpaca/scripts/submodules/alpaca_input_formatting"

Rscript "${SCRIPT_DIR}/convert_conipher_output/convert_conipher_output.R" \
    --CONIPHER_tree_object $CONIPHER_tree_object \
    --output_dir $output_dir
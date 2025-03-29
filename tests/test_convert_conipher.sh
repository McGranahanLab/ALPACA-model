#!/bin/bash


output_dir="tests/test_output"
CONIPHER_tree_object="examples/example_cohort/input/LTX0530/LTX0530.tree.RDS"
SCRIPT_DIR="alpaca/scripts/submodules/alpaca_input_formatting"
Rscript "${SCRIPT_DIR}/convert_conipher_output/convert_conipher_output.R" \
    --CONIPHER_tree_object $CONIPHER_tree_object \
    --output_dir $output_dir
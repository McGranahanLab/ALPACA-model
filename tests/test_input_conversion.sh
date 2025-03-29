#!/bin/bash


output_dir="tests/test_output"
CONIPHER_tree_object="examples/example_cohort/input/LTX0530/LTX0530.tree.RDS"
refphase_rData="examples/example_cohort/input/LTX0530/LTX0530-refphase-results.RData"
tumour_id='LTX0530'

SCRIPT_DIR="alpaca/scripts/submodules/alpaca_input_formatting"
bash "${SCRIPT_DIR}/input_conversion.sh" \
    --refphase_rData $refphase_rData \
    --tumour_id $tumour_id \
    --CONIPHER_tree_object $CONIPHER_tree_object \
    --output_dir $output_dir
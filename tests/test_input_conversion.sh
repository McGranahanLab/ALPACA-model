#!/bin/bash
tumour_id=LTX0000-Tumour1
output_dir="tests/test_output"
CONIPHER_tree_object="examples/example_cohort/input/${tumour_id}/${tumour_id}.tree.RDS"
refphase_rData="examples/example_cohort/input/${tumour_id}/${tumour_id}-refphase-results.RData"
SCRIPT_DIR="alpaca/scripts/submodules/alpaca_input_formatting"

bash "${SCRIPT_DIR}/input_conversion.sh" \
    --refphase_rData $refphase_rData \
    --tumour_id $tumour_id \
    --CONIPHER_tree_object $CONIPHER_tree_object \
    --output_dir $output_dir
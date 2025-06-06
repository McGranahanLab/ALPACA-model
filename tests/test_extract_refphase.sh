#!/bin/bash
tumour_id=LTX0000-Tumour1
output_dir="tests/test_output"
refphase_rData="examples/example_cohort/input/${tumour_id}/${tumour_id}-refphase-results.RData"
SCRIPT_DIR="alpaca/scripts/submodules/alpaca_input_formatting"

Rscript "${SCRIPT_DIR}/convert_refphase_output/extract_rephase_data.R" \
    --refphase_rData $refphase_rData \
    --output_dir $output_dir

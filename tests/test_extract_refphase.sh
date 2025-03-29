#!/bin/bash

output_dir="tests/test_output"
refphase_rData="examples/example_cohort/input/LTX0530/LTX0530-refphase-results.RData"
SCRIPT_DIR="alpaca/scripts/submodules/alpaca_input_formatting"


Rscript "${SCRIPT_DIR}/convert_refphase_output/extract_rephase_data.R" \
    --refphase_rData $refphase_rData \
    --output_dir $output_dir

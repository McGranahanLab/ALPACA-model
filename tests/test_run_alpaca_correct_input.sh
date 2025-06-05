#!/bin/bash
tumour_id=LTX0000-Tumour1
input_tumour_directory="examples/example_cohort/input/${tumour_id}"
output_directory="tests/correct_input/output/${tumour_id}"

# run alpaca:
alpaca run \
    --input_tumour_directory "${input_tumour_directory}" \
    --output_directory "${output_directory}" \
    --overwrite_output 0
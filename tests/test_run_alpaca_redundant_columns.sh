#!/bin/bash
tumour_id=LTX0000-Tumour1
input_tumour_directory="tests/redundant_columns/input/${tumour_id}"
output_directory="tests/redundant_columns/output/${tumour_id}"

# run alpaca:
alpaca run \
    --input_tumour_directory "${input_tumour_directory}" \
    --output_directory "${output_directory}"
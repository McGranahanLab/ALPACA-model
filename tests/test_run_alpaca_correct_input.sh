#!/bin/bash
tumour_id=LTX0000-Tumour1
input_tumour_directory="test/correct_input/input/${tumour_id}"
output_directory="test/correct_input/output/${tumour_id}"

# run alpaca:
alpaca run \
    --input_tumour_directory "${input_tumour_directory}" \
    --output_directory "${output_directory}"
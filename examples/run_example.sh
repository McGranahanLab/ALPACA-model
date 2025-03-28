#!/bin/bash
tumour_id=LTX0530
input_tumour_directory="example_cohort/input/${tumour_id}"
output_directory="example_cohort/output/${tumour_id}"

echo "Tumour ID: ${tumour_id}"

alpaca \
    --input_tumour_directory "${input_tumour_directory}" \
    --output_directory "${output_directory}"

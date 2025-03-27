#!/bin/bash
tumour_id=LTXSIM039
input_tumour_directory="../example_cohort/input/${tumour_id}"
output_directory="../example_cohort/output/${tumour_id}"

echo "Tumour ID: ${tumour_id}"

python3 -m alpaca \
    --input_tumour_directory "${input_tumour_directory}" \
    --output_directory "${output_directory}" \
    --ci_table_name ''



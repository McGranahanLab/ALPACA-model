#!/bin/bash
tumour_id=LTX0000-Tumour1
SCRIPT_DIR="alpaca/scripts/submodules/alpaca_input_formatting"
output_dir="tests/test_output"
refphase_segments_path="tests/test_output/phased_segs_test.tsv" # just two segments to speed up testing
refphase_snps_path="tests/test_output/phased_snps.tsv"
refphase_purity_ploidy_path="tests/test_output/purity_ploidy.tsv"

python3 "${SCRIPT_DIR}/convert_refphase_output/convert_refphase.py" \
    --tumour_id $tumour_id \
    --output_dir $output_dir \
    --refphase_segments $refphase_segments_path \
    --refphase_snps $refphase_snps_path \
    --refphase_purity_ploidy $refphase_purity_ploidy_path

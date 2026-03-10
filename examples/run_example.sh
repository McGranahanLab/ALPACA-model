#!/bin/bash
set -euo pipefail
tumour_id=LTX0000-Tumour1
input_tumour_directory="examples/example_cohort/input/${tumour_id}"
output_directory="examples/example_cohort/output/${tumour_id}"

echo "Tumour ID: ${tumour_id}"

# convert CONIPHER and Refphase outputs to ALPACA input:

refphase_rData="${input_tumour_directory}/${tumour_id}-refphase-results.RData"
CONIPHER_tree_object="${input_tumour_directory}/${tumour_id}.tree.RDS"
conversion_output_dir="${input_tumour_directory}"
CONIPHER_tree_index=1
heterozygous_SNPs_threshold=0

alpaca input-conversion \
 --tumour_id $tumour_id \
 --copy_number_tool refphase \
 --refphase_rData $refphase_rData \
 --CONIPHER_tree_object $CONIPHER_tree_object \
 --CONIPHER_tree_index $CONIPHER_tree_index \
 --output_dir $conversion_output_dir \
 --heterozygous_SNPs_threshold $heterozygous_SNPs_threshold \
 --recalculate_updated_cns 0 \
 --recalculate_not_updated_cns 0 \
 --recalculate_reference_cns 1
 
# run alpaca:
alpaca run \
    --input_tumour_directory "${input_tumour_directory}" \
    --output_directory "${output_directory}" \
    --plot_output_mode pdf \
    --genome_build hg19

# get cn change to ancestor:
alpaca ancestor-delta \
    --output_directory "${output_directory}" \
    --tumour_df_path "${output_directory}/ALPACA_output_${tumour_id}.csv" \
    --tree_path "${input_tumour_directory}/tree_paths.json"

# calculate clone copy number diversity:
alpaca ccd \
    --output_directory "${output_directory}" \
    --alpaca_output_path "${output_directory}/ALPACA_output_${tumour_id}.csv"

# calculate WGD ratios:
alpaca wgd \
    --output_directory "${output_directory}" \
    --alpaca_output_path "${output_directory}/ALPACA_output_${tumour_id}.csv" \
    --tree_path "${input_tumour_directory}/tree_paths.json"

alpaca plot-tumour \
    --input_directory "${input_tumour_directory}" \
    --output_directory "${output_directory}" \
    --alpaca_output_path "${output_directory}/ALPACA_output_${tumour_id}.csv" \
    --plot_output_mode notebook \
    --notebook_name "example_notebook" \
    --heatmap_palette "magma" \
    --genome_build hg19

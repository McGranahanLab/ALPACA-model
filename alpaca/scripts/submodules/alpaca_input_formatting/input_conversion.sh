#!/usr/bin/env bash

# Fail fast: exit on error, unset variables are errors, propagate failures in pipelines
set -euo pipefail
IFS=$'\n\t'

# Print a helpful message and exit when any command fails
trap 'code=$?; cmd="${BASH_COMMAND:-unknown}"; echo "ERROR: command failed: \"$cmd\" (exit code $code)" >&2; exit $code' ERR

# Ensure script is run with bash (not sh)
if [ -z "${BASH_VERSION:-}" ]; then
    echo "This script requires bash. Run with 'bash $0' or make executable and run './$0'" >&2
    exit 2
fi
# default arguments

usage() {
    echo "Usage: $0 --tumour_id TUMOUR_ID --refphase_rData RDATA_PATH --CONIPHER_tree_object TREE_OBJECT_PATH --output_dir OUTPUT_DIR"
    echo
    echo "Arguments:"
    echo "  --tumour_id              Tumour ID (required)"
    echo "  --refphase_rData         Path to refphase .RData file (required)"
    echo "  --CONIPHER_tree_object   Path to CONIPHER tree object .RDS file (required)"
    echo "  --CONIPHER_tree_index    Selected CONIPHER tree index (optional, default: 1)"
    echo "  --heterozygous_SNPs_threshold  Optional int threshold passed to convert_refphase.py - default value is 5"
    echo "  --ci_value               Confidence interval level (float, default: 0.5)"
    echo "  --n_bootstrap            Number of bootstrap iterations used by conversion steps while calculating confidence intervals (int, default: 100)"
    echo "  --recalculate_not_updated_cns  If set to 1, forces recalculation of copy numbers for segments that were not updated by Refphase (default: 0)"
    echo "  --recalculate_updated_cns      If set to 1, forces recalculation of copy numbers for segments that were updated by Refphase (default: 0)"
    echo "  --recalculate_reference_cns    If set to 1, forces recalculation of copy numbers for reference segments (default: 0)"
    echo "  --debug                  If set, writes a detailed debug report for convert_refphase.py (default: off)"
    echo "  --output_dir             Output directory (required)"
    echo "  --help                   Display this help message"
    exit 1
}
# default arguments
CONIPHER_tree_index=1
heterozygous_SNPs_threshold=5
ci_value=0.5
n_bootstrap=100
recalculate_not_updated_cns=0
recalculate_updated_cns=0
recalculate_reference_cns=0
debug=0
# Parse named arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --tumour_id) tumour_id="$2"; shift ;;
        --refphase_rData) refphase_rData="$2"; shift ;;
        --CONIPHER_tree_object) CONIPHER_tree_object="$2"; shift ;;
        --CONIPHER_tree_index) CONIPHER_tree_index="$2"; shift ;;
        --heterozygous_SNPs_threshold) heterozygous_SNPs_threshold="$2"; shift ;;
        --ci_value) ci_value="$2"; shift ;;
        --n_bootstrap) n_bootstrap="$2"; shift ;;
        --recalculate_not_updated_cns) recalculate_not_updated_cns="$2"; shift ;;
        --recalculate_updated_cns) recalculate_updated_cns="$2"; shift ;;
        --recalculate_reference_cns) recalculate_reference_cns="$2"; shift ;;
        --debug) debug=1 ;;
        --output_dir) output_dir="$2"; shift ;;
        --help) usage ;;
        *) echo "Unknown parameter passed: $1"; usage ;;
    esac
    shift
done

# Check if required arguments are provided
if [ -z "$tumour_id" ] || [ -z "$refphase_rData" ] || [ -z "$CONIPHER_tree_object" ] || [ -z "$output_dir" ]; then
    echo "Error: All arguments are required"
    usage
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"


echo "Tumour ID: ${tumour_id}"
echo "Extract TSV files from RData"
echo "make output directory: $output_dir"
mkdir -p $output_dir

# exit if refphase_rData does not exist:
if [ ! -f $refphase_rData ]; then
    echo "Error: refphase_rData (${refphase_rData}) file does not exist" >&2
    exit 1
fi

if ! command -v python3 &> /dev/null; then
    echo "Error: python3 is not available in the environment." >&2
    exit 1
fi
echo "===================================="
echo "Extracting data from REFPHASE output"
python3 "${SCRIPT_DIR}/convert_refphase_output/extract_rephase_data.py" \
    --refphase_rData $refphase_rData \
    --output_dir $output_dir
if [ $? -ne 0 ]; then
    echo "Python extract_rephase_data.py failed" >&2
    exit 1
fi

refphase_segments_path="${output_dir}/phased_segs.tsv"
refphase_snps_path="${output_dir}/phased_snps.tsv"
refphase_purity_ploidy_path="${output_dir}/purity_ploidy.tsv"
echo "===================================="
echo "Extracting data from CONIPHER output"
python3 "${SCRIPT_DIR}/convert_conipher_output/convert_conipher_output.py" \
    --CONIPHER_tree_object $CONIPHER_tree_object \
    --CONIPHER_tree_index $CONIPHER_tree_index \
    --output_dir $output_dir
if [ $? -ne 0 ]; then
    echo "Python convert_conipher_output.py failed" >&2
    exit 1
fi
echo "===================================="
echo "Converting REFPHASE output"
convert_refphase_cmd=(
    python3 "${SCRIPT_DIR}/convert_refphase_output/convert_refphase.py"
    --tumour_id "${tumour_id}"
    --output_dir "${output_dir}"
    --refphase_segments "${refphase_segments_path}"
    --refphase_snps "${refphase_snps_path}"
    --refphase_purity_ploidy "${refphase_purity_ploidy_path}"
    --conipher_cp_table "${output_dir}/cp_table.csv"
    --heterozygous_SNPs_threshold "${heterozygous_SNPs_threshold}"
    --ci_value "${ci_value}"
    --n_bootstrap "${n_bootstrap}"
    --recalculate_not_updated_cns "${recalculate_not_updated_cns}"
    --recalculate_updated_cns "${recalculate_updated_cns}"
    --recalculate_reference_cns "${recalculate_reference_cns}"
)
if [ "${debug}" -eq 1 ]; then
    convert_refphase_cmd+=(--debug)
fi
"${convert_refphase_cmd[@]}"

if [ $? -ne 0 ]; then
    echo "Python convert_refphase.py failed" >&2
    exit 1
fi

# Ensure that the final output (cp_table.csv and ALPACA_input_table.csv) contain
# the same set of samples.

python3 "${SCRIPT_DIR}/check_final_outputs.py" \
    --output_dir "${output_dir}"


# Write a small conversion report summarising the arguments used and the run time
report_path="${output_dir%/}/conversion_report.txt"
sample_filter_report_path="${output_dir%/}/sample_filter_report.txt"
cat > "$report_path" <<-REPORT
Conversion report
=================
Date: $(date -u +"%Y-%m-%dT%H:%M:%SZ")
Script: ${SCRIPT_DIR}/convert_refphase_output/convert_refphase.py

Arguments and values:
    tumour_id: ${tumour_id}
    refphase_rData: ${refphase_rData}
    CONIPHER_tree_object: ${CONIPHER_tree_object}
    CONIPHER_tree_index: ${CONIPHER_tree_index}
    heterozygous_SNPs_threshold: ${heterozygous_SNPs_threshold}
    ci_value: ${ci_value}
    n_bootstrap: ${n_bootstrap}
    recalculate_not_updated_cns: ${recalculate_not_updated_cns}
    recalculate_updated_cns: ${recalculate_updated_cns}
    recalculate_reference_cns: ${recalculate_reference_cns}
    debug: ${debug}
    output_dir: ${output_dir}

REPORT

if [ -f "$sample_filter_report_path" ]; then
    {
        echo
        echo "Sample filtering report:"
        sed 's/^/    /' "$sample_filter_report_path"
    } >> "$report_path"
fi

echo "Wrote conversion report to $report_path"

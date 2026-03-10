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
    echo "Usage: $0 --tumour_id TUMOUR_ID --CONIPHER_tree_object TREE_OBJECT_PATH --output_dir OUTPUT_DIR [--copy_number_tool refphase|battenberg] [...]"
    echo
    echo "Arguments:"
    echo "  --tumour_id              Tumour ID (required)"
    echo "  --copy_number_tool       Copy number source tool: refphase|battenberg (optional, default: refphase)"
    echo "  --chromosome             Optional chromosome filter (e.g. 1, chr1, X). If set, only this chromosome is processed."
    echo "  --refphase_rData         Path to refphase .RData file (required when --copy_number_tool refphase)"
    echo "  --battenberg_inventory   Path to Battenberg inventory file (recommended for --copy_number_tool battenberg)"
    echo "  --battenberg_input_dir   Path to Battenberg directory for auto-discovery (optional for --copy_number_tool battenberg)"
    echo "  --CONIPHER_tree_object   Path to CONIPHER tree object .RDS file (required)"
    echo "  --CONIPHER_tree_index    Selected CONIPHER tree index (optional, default: 1)"
    echo "  --heterozygous_SNPs_threshold  Optional int threshold passed to convert_refphase.py - default value is 5"
    echo "  --ci_value               Confidence interval level (float, default: 0.5)"
    echo "  --n_bootstrap            Number of bootstrap iterations used by conversion steps while calculating confidence intervals (int, default: 100)"
    echo "  --recalculate_not_updated_cns  If set to 1, forces recalculation of copy numbers for segments that were not updated by Refphase (default: 0)"
    echo "  --recalculate_updated_cns      If set to 1, forces recalculation of copy numbers for segments that were updated by Refphase (default: 0)"
    echo "  --recalculate_reference_cns    If set to 1, forces recalculation of copy numbers for reference segments (default: 0)"
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
copy_number_tool="refphase"
chromosome=""
battenberg_inventory=""
battenberg_input_dir=""
# Parse named arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --tumour_id) tumour_id="$2"; shift ;;
        --copy_number_tool) copy_number_tool="$2"; shift ;;
        --chromosome) chromosome="$2"; shift ;;
        --refphase_rData) refphase_rData="$2"; shift ;;
        --battenberg_inventory) battenberg_inventory="$2"; shift ;;
        --battenberg_input_dir) battenberg_input_dir="$2"; shift ;;
        --CONIPHER_tree_object) CONIPHER_tree_object="$2"; shift ;;
        --CONIPHER_tree_index) CONIPHER_tree_index="$2"; shift ;;
        --heterozygous_SNPs_threshold) heterozygous_SNPs_threshold="$2"; shift ;;
        --ci_value) ci_value="$2"; shift ;;
        --n_bootstrap) n_bootstrap="$2"; shift ;;
        --recalculate_not_updated_cns) recalculate_not_updated_cns="$2"; shift ;;
        --recalculate_updated_cns) recalculate_updated_cns="$2"; shift ;;
        --recalculate_reference_cns) recalculate_reference_cns="$2"; shift ;;
        --output_dir) output_dir="$2"; shift ;;
        --help) usage ;;
        *) echo "Unknown parameter passed: $1"; usage ;;
    esac
    shift
done

copy_number_tool="$(echo "${copy_number_tool}" | tr '[:upper:]' '[:lower:]')"
if [[ "${copy_number_tool}" != "refphase" && "${copy_number_tool}" != "battenberg" ]]; then
    echo "Error: --copy_number_tool must be one of: refphase, battenberg" >&2
    usage
fi

# Check if required arguments are provided
if [ -z "${tumour_id:-}" ] || [ -z "${CONIPHER_tree_object:-}" ] || [ -z "${output_dir:-}" ]; then
    echo "Error: All arguments are required"
    usage
fi

if [ "${copy_number_tool}" = "refphase" ] && [ -z "${refphase_rData:-}" ]; then
    echo "Error: --refphase_rData is required when --copy_number_tool refphase" >&2
    usage
fi

if [ "${copy_number_tool}" = "battenberg" ] && [ -z "${battenberg_inventory:-}" ] && [ -z "${battenberg_input_dir:-}" ]; then
    echo "Error: for --copy_number_tool battenberg, provide --battenberg_inventory and/or --battenberg_input_dir" >&2
    usage
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"


echo "Tumour ID: ${tumour_id}"
echo "Copy number tool: ${copy_number_tool}"
echo "make output directory: $output_dir"
mkdir -p $output_dir

if [ "${copy_number_tool}" = "refphase" ]; then
    # exit if refphase_rData does not exist:
    if [ ! -f "${refphase_rData}" ]; then
        echo "Error: refphase_rData (${refphase_rData}) file does not exist" >&2
        exit 1
    fi
fi

if ! command -v python3 &> /dev/null; then
    echo "Error: python3 is not available in the environment." >&2
    exit 1
fi
echo "===================================="
if [ "${copy_number_tool}" = "refphase" ]; then
    echo "Extracting data from REFPHASE output"
    refphase_extract_cmd=(
        python3 "${SCRIPT_DIR}/convert_refphase_output/extract_rephase_data.py"
        --refphase_rData "$refphase_rData"
        --output_dir "$output_dir"
    )
    if [ -n "${chromosome:-}" ]; then
        refphase_extract_cmd+=(--chromosome "$chromosome")
    fi
    "${refphase_extract_cmd[@]}"
    if [ $? -ne 0 ]; then
        echo "Python extract_rephase_data.py failed" >&2
        exit 1
    fi
else
    echo "Extracting data from BATTENBERG output"
    battenberg_cmd=(
        python3 "${SCRIPT_DIR}/convert_battenberg_output/extract_battenberg_data.py"
        --tumour_id "$tumour_id"
        --output_dir "$output_dir"
    )
    if [ -n "${chromosome:-}" ]; then
        battenberg_cmd+=(--chromosome "$chromosome")
    fi
    if [ -n "${battenberg_inventory:-}" ]; then
        battenberg_cmd+=(--battenberg_inventory "$battenberg_inventory")
    fi
    if [ -n "${battenberg_input_dir:-}" ]; then
        battenberg_cmd+=(--battenberg_input_dir "$battenberg_input_dir")
    fi
    "${battenberg_cmd[@]}"
    if [ $? -ne 0 ]; then
        echo "Python extract_battenberg_data.py failed" >&2
        exit 1
    fi
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
echo "Converting copy number output (${copy_number_tool})"
convert_cmd=(
    python3 "${SCRIPT_DIR}/convert_refphase_output/convert_refphase.py"
    --tumour_id "$tumour_id"
    --output_dir "$output_dir"
    --copy_number_tool "$copy_number_tool"
    --refphase_segments "$refphase_segments_path"
    --refphase_snps "$refphase_snps_path"
    --refphase_purity_ploidy "$refphase_purity_ploidy_path"
    --conipher_cp_table "${output_dir}/cp_table.csv"
    --heterozygous_SNPs_threshold "${heterozygous_SNPs_threshold}"
    --ci_value "${ci_value}"
    --n_bootstrap "${n_bootstrap}"
    --recalculate_not_updated_cns "${recalculate_not_updated_cns}"
    --recalculate_updated_cns "${recalculate_updated_cns}"
    --recalculate_reference_cns "${recalculate_reference_cns}"
)
if [ -n "${chromosome:-}" ]; then
    convert_cmd+=(--chromosome "$chromosome")
fi
"${convert_cmd[@]}"

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
    copy_number_tool: ${copy_number_tool}
    chromosome: ${chromosome:-}
    refphase_rData: ${refphase_rData:-}
    battenberg_inventory: ${battenberg_inventory:-}
    battenberg_input_dir: ${battenberg_input_dir:-}
    CONIPHER_tree_object: ${CONIPHER_tree_object}
    CONIPHER_tree_index: ${CONIPHER_tree_index}
    heterozygous_SNPs_threshold: ${heterozygous_SNPs_threshold}
    ci_value: ${ci_value}
    n_bootstrap: ${n_bootstrap}
    recalculate_not_updated_cns: ${recalculate_not_updated_cns}
    recalculate_updated_cns: ${recalculate_updated_cns}
    recalculate_reference_cns: ${recalculate_reference_cns}
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

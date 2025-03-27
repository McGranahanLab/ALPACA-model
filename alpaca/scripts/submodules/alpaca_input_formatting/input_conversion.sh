

# default arguments


usage() {
    echo "Usage: $0 --tumour_id TUMOUR_ID --refphase_rData RDATA_PATH --CONIPHER_tree_object TREE_OBJECT_PATH --output_dir OUTPUT_DIR"
    echo
    echo "Arguments:"
    echo "  --tumour_id              Tumour ID (required)"
    echo "  --refphase_rData         Path to refphase .RData file (required)"
    echo "  --CONIPHER_tree_object   Path to CONIPHER tree object .RDS file (required)"
    echo "  --output_dir             Output directory (required)"
    echo "  --help                   Display this help message"
    exit 1
}

# Parse named arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --tumour_id) tumour_id="$2"; shift ;;
        --refphase_rData) refphase_rData="$2"; shift ;;
        --CONIPHER_tree_object) CONIPHER_tree_object="$2"; shift ;;
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
    echo "Error: refphase_rData (${refphase_rData}) file does not exist"
    exit 1
fi

echo "Extracting data from REFPHASE output"
Rscript "${SCRIPT_DIR}/convert_refphase_output/extract_rephase_data.R" \
    --refphase_rData $refphase_rData \
    --output_dir $output_dir

refphase_segments_path="${output_dir}/phased_segs.tsv"
refphase_snps_path="${output_dir}/phased_snps.tsv"
refphase_purity_ploidy_path="${output_dir}/purity_ploidy.tsv"

echo "Converting REFPHASE output"
python3 "${SCRIPT_DIR}/convert_refphase_output/convert_refphase.py" \
    --tumour_id $tumour_id \
    --output_dir $output_dir \
    --refphase_segments $refphase_segments_path \
    --refphase_snps $refphase_snps_path \
    --refphase_purity_ploidy $refphase_purity_ploidy_path

echo "Extracting data from CONIPHER output"
Rscript "${SCRIPT_DIR}/convert_conipher_output/convert_conipher_output.R" \
    --CONIPHER_tree_object $CONIPHER_tree_object \
    --output_dir $output_dir
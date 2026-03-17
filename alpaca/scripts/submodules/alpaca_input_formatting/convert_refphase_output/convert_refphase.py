import pandas as pd
import argparse
import os
import re
import numpy as np
from functions import calculate_confidence_intervals, get_consensus_segmentation, calibrate_battenberg_cns_and_cis

# arguments
parser = argparse.ArgumentParser(
    description="Calculate confidence intervals from copy-number tool output"
)
parser.add_argument(
    "--tumour_id", type=str, help="Unique identifier for the tumour", required=True
)
parser.add_argument("--output_dir", type=str, help="Output directory", required=True)
parser.add_argument(
    "--chromosome",
    type=str,
    required=False,
    help="Optional chromosome filter (e.g. 1, chr1, X). If set, only this chromosome is processed.",
)
parser.add_argument(
    "--copy_number_tool",
    type=str,
    choices=["refphase", "battenberg"],
    default="refphase",
    help="Copy-number input source. Default is refphase for backwards compatibility.",
)
parser.add_argument(
    "--refphase_segments",
    type=str,
    help="Location of intermediate segments file",
    required=True,
)
parser.add_argument(
    "--refphase_snps", type=str, help="Location of intermediate SNPs file", required=True
)
parser.add_argument(
    "--refphase_purity_ploidy",
    type=str,
    help="Location of intermediate purity ploidy file",
    required=True,
)

parser.add_argument(
    "--conipher_cp_table",
    type=str,
    help="Path to CONIPHER cp_table CSV (required)",
    required=True,
)

# options
parser.add_argument(
    "--heterozygous_SNPs_threshold",
    type=int,
    default=5,
    help="Minimum number of heterozygous SNPs to consider a segment. Segments with fewer heterozygous SNPs will be discarded.",
)
parser.add_argument("--ci_value", type=float, help="Confidence interval value.")
parser.add_argument("--n_bootstrap", type=int, help="Number of bootstrap samples.")
parser.add_argument(
    "--recalculate_not_updated_cns",
    type=int,
    choices=[0, 1],
    default=0,
    help="Refphase updates copy-numbers for segments where allelic imbalance is detected. \
        The remaining segments inherit the copy-number of their parent ASCAT segment. \
        When calculating confidence intervals for these non-updated segments, two behaviours are possible. \
        If set to 1, we will recalculate confidence intervals and fractional copy-numbers for these segments using BAF and LOGr of the subset of SNPs\
        assigned  to the Refphase segment in questions. Otherwise, we will first center the SNPs around the original ASCAT copy-numbers, and then calculate\
        confidence intervals. The rationale for such behaviour is that in the second case, there is not enough evidence to divert from the null\
        (i.e. ASCAT solution), but the uncertainty in the copy-number estimate should still be captured and should be lower compared to the entire\
        parent ASCAT segment",
)
parser.add_argument(
    "--recalculate_updated_cns",
    type=int,
    choices=[0, 1],
    default=0,
    help="Refphase updates copy-numbers for segments where allelic imbalance is detected. \
        While doing so, it uses ASCAT equations to calculate CNS based on BAF, LOG, purity, ploidy etc. \
        Since we are using the same data and equations to caclculate confidence intervals, we can also re-calculate the original copy number as well.\
        However, for many segments, such recalculated copy number differs slightly from the value provided by the refphase. If this argument is 0, \
        instead of calculating the copy number, we will just calculate the intervals and center them around the original refphase provided value",
)

parser.add_argument(
    "--recalculate_reference_cns",
    type=int,
    choices=[0, 1],
    default=0,
    help="Recalculates the copy-number for segments marked as 'is_reference' True in Refphase. \
        Default refphase behaviour is to recalculate and then round these copy numbers to nearest integers. \
        Setting this option to '1' will trigger recalculation without the rounding, i.e. leaving the copy number for these segments in fractional state",
)

parser.add_argument(
    "--split_segments",
    type=int,
    choices=[0, 1],
    default=0,
    help="Split input into separate files for each segment. Useful for parallel processing.",
)

# placeholder for future argument - currently always set to True
BATTENBERG_RECALIBRATE_CI = True

args = parser.parse_args()
tumour_id = args.tumour_id
output_dir = args.output_dir
copy_number_tool = args.copy_number_tool
chromosome = args.chromosome
ci_value = args.ci_value
n_bootstrap = args.n_bootstrap
recalculate_not_updated_cns = bool(args.recalculate_not_updated_cns)
recalculate_updated_cns = bool(args.recalculate_updated_cns)
recalculate_reference_cns = bool(args.recalculate_reference_cns)
split_segments = bool(args.split_segments)

if copy_number_tool == "battenberg":
    # Battenberg path always recalculates copy numbers for all segments.
    recalculate_not_updated_cns = True
    recalculate_updated_cns = False
    recalculate_reference_cns = False
    print(
        "Battenberg mode enabled: forcing recalculate_not_updated_cns=1, "
        "recalculate_updated_cns=0, recalculate_reference_cns=0."
    )

# create output directory:
os.makedirs(output_dir, exist_ok=True)
# read data
refphase_segments = pd.read_csv(args.refphase_segments, sep="\t")
refphase_snps = pd.read_csv(args.refphase_snps, sep="\t")
refphase_purity_ploidy = pd.read_csv(args.refphase_purity_ploidy, sep="\t")
cp_table = pd.read_csv(args.conipher_cp_table, index_col="clone")
conipher_samples = cp_table.columns
# rename columns:
refphase_segments = refphase_segments.rename(
    columns={
        "group_name": "sample",
        "seqnames": "chr",
        "patient_tumour": "tumour_id",
    }
)
refphase_snps = refphase_snps.rename(
    columns={
        "group_name": "sample",
        "seqnames": "chr",
        "patient_tumour": "tumour_id",
    }
)

if "was_cn_updated" not in refphase_segments.columns:
    refphase_segments["was_cn_updated"] = False
if "is_reference" not in refphase_segments.columns:
    refphase_segments["is_reference"] = False
if "phasing" not in refphase_snps.columns:
    refphase_snps["phasing"] = "a"


# sanitize chromosome names:
def _sanitize_chr_names(df):
    if pd.api.types.is_numeric_dtype(df["chr"]):
        return
    # Extract the actual chromosome identifier (letters/numbers only, ignoring prefixes)
    df["chr"] = df["chr"].str.extract(r"([0-9]+|[XYM][Tt]?)", flags=re.IGNORECASE)[0]
    # Map to numbers
    chr_map = {
        "X": "23",
        "x": "23",
        "Y": "24",
        "y": "24",
        "MT": "25",
        "Mt": "25",
        "mt": "25",
        "M": "25",
        "m": "25",
    }
    df["chr"] = df["chr"].replace(chr_map)
    df["chr"] = pd.to_numeric(df["chr"], errors="coerce")


def _normalize_chr_value(chr_value):
    s = str(chr_value)
    extracted = re.findall(r"([0-9]+|[XYM][Tt]?)", s, flags=re.IGNORECASE)
    if not extracted:
        raise ValueError(f"Could not parse chromosome value: {chr_value}")
    token = extracted[0].upper()
    token = {"X": "23", "Y": "24", "MT": "25", "M": "25"}.get(token, token)
    return int(token)


_sanitize_chr_names(refphase_segments)
_sanitize_chr_names(refphase_snps)

if chromosome is not None:
    target_chr = _normalize_chr_value(chromosome)
    refphase_segments = refphase_segments[refphase_segments["chr"] == target_chr].copy()
    refphase_snps = refphase_snps[refphase_snps["chr"] == target_chr].copy()
    if refphase_segments.empty:
        raise ValueError(
            f"No segment rows remain after --chromosome filter ({chromosome})."
        )
    if refphase_snps.empty:
        raise ValueError(
            f"No SNP rows remain after --chromosome filter ({chromosome})."
        )

# create segment column by combining chromosome, start and end:
refphase_segments["segment"] = (
    refphase_segments["chr"].astype(str)
    + "_"
    + refphase_segments["start"].astype(str)
    + "_"
    + refphase_segments["end"].astype(str)
)

# remove segments with fewer than args.heterozygous_SNPs_threshold heterozygous SNPs:
# use this to use total number accross all samples:
# SNP_count = refphase_segments.groupby('segment')['heterozygous_SNP_number'].sum().sort_values()
# segments_above_threshold = SNP_count[SNP_count>args.heterozygous_SNPs_threshold]

# use this to use number of heterozygous SNPs in each sample:
if copy_number_tool == "refphase":
    refphase_segments = refphase_segments.groupby("segment").filter(
        lambda x: (x["heterozygous_SNP_number"] >= args.heterozygous_SNPs_threshold).all()
    )
else:
    # Battenberg segments are sample-specific; thresholding should be applied per-row.
    refphase_segments = refphase_segments[
        refphase_segments["heterozygous_SNP_number"] >= args.heterozygous_SNPs_threshold
    ]
# NB: if segments with 0 SNPs are kept, the copy number for such segments will not be recalculated with calculate_confidence_intervals
# and confidence interval values will be set to equal the input copy number values.

# calculate confidence intervals:
print(f"Calculating confidence intervals for {tumour_id}")
# assign SNPS to segments:
snps_with_segments = refphase_snps.merge(
    refphase_segments,
    left_on=["sample", "chr"],
    right_on=["sample", "chr"],
    how="inner",
)

snps_with_segments = snps_with_segments[
    (snps_with_segments["pos"] >= snps_with_segments["start"])
    & (snps_with_segments["pos"] <= snps_with_segments["end"])
]
# add purity and ploidy information
snps_with_segments_purity_ploidy = snps_with_segments.merge(
    refphase_purity_ploidy, left_on="sample", right_on="sample_id", how="inner"
)

# estimate the confidence intervals:
confidence_intervals = (
    snps_with_segments_purity_ploidy.groupby(["segment", "sample"])
    .apply(
        calculate_confidence_intervals,
        ci_value=ci_value,
        n_bootstrap=n_bootstrap,
        recalculate_not_updated_cns=recalculate_not_updated_cns,
        recalculate_updated_cns=recalculate_updated_cns,
        recalculate_reference_cns=recalculate_reference_cns,
    )
    .reset_index()
    .drop(columns=["level_2"])
)
# add 0 SNP segments if such segments not filtereted out, i.e when args.heterozygous_SNPs_threshold=0
if args.heterozygous_SNPs_threshold == 0:
    # add 0 SNP segments if such segments not filtereted out, i.e when args.heterozygous_SNPs_threshold=0
    # for such segments, we set confidence intervals to 0.5 to reflect low certainty
    CI_span = 0.5
    zero_snp_segments = refphase_segments[
        refphase_segments.heterozygous_SNP_number == 0
    ]
    zero_snp_segments = zero_snp_segments[["segment", "sample", "cn_a", "cn_b"]]
    CI_half = CI_span / 2.0
    zero_snp_segments["lower_CI_A"] = (zero_snp_segments["cn_a"] - CI_half).clip(
        lower=0
    )
    zero_snp_segments["upper_CI_A"] = zero_snp_segments["cn_a"] + CI_half
    zero_snp_segments["lower_CI_B"] = (zero_snp_segments["cn_b"] - CI_half).clip(
        lower=0
    )
    zero_snp_segments["upper_CI_B"] = zero_snp_segments["cn_b"] + CI_half
    zero_snp_segments.rename(columns={"cn_a": "cpnA", "cn_b": "cpnB"}, inplace=True)
    confidence_intervals = pd.concat(
        [confidence_intervals, zero_snp_segments], ignore_index=True
    )
confidence_intervals["chr"] = confidence_intervals["segment"].apply(
    lambda x: int(x.split("_")[0])
)
confidence_intervals["start"] = confidence_intervals["segment"].apply(
    lambda x: int(x.split("_")[1])
)
LOW_SNP_THRESHOLD = 10
# we don't trust confidence intervals for segments with very low number of heterozygous SNPs, 
# ensure such segments have confidence intervals of at least 0.5 (i.e. +/- 0.25 around the original copy number value):
low_snp_ci_span = 0.5
low_snp_ci_half = low_snp_ci_span / 2.0
low_snp_segments = pd.MultiIndex.from_frame(
    refphase_segments.loc[
        refphase_segments["heterozygous_SNP_number"] <= LOW_SNP_THRESHOLD, ["segment", "sample"]
    ].drop_duplicates()
)
if len(low_snp_segments) > 0:
    low_snp_mask = pd.MultiIndex.from_frame(
        confidence_intervals[["segment", "sample"]]
    ).isin(low_snp_segments)
    confidence_intervals.loc[low_snp_mask, "lower_CI_A"] = (
        confidence_intervals.loc[low_snp_mask, "cpnA"] - low_snp_ci_half
    ).clip(lower=0)
    confidence_intervals.loc[low_snp_mask, "upper_CI_A"] = (
        confidence_intervals.loc[low_snp_mask, "cpnA"] + low_snp_ci_half
    )
    confidence_intervals.loc[low_snp_mask, "lower_CI_B"] = (
        confidence_intervals.loc[low_snp_mask, "cpnB"] - low_snp_ci_half
    ).clip(lower=0)
    confidence_intervals.loc[low_snp_mask, "upper_CI_B"] = (
        confidence_intervals.loc[low_snp_mask, "cpnB"] + low_snp_ci_half
    )

confidence_intervals = confidence_intervals.sort_values(by=["sample", "chr", "start"])
confidence_intervals.drop(columns=["chr", "start"], inplace=True)
if copy_number_tool == "battenberg" and BATTENBERG_RECALIBRATE_CI:    
    confidence_intervals = calibrate_battenberg_cns_and_cis(confidence_intervals, refphase_segments)

ci_table = confidence_intervals.merge(refphase_segments)[
    [
        "segment",
        "sample",
        "cn_a",
        "cn_b",
        "cpnA",
        "cpnB",
        "lower_CI_A",
        "upper_CI_A",
        "lower_CI_B",
        "upper_CI_B",
        "was_cn_updated",
        *(
            ["is_clonal"]
            if "is_clonal" in refphase_segments.columns
            else []
        ),
    ]
].drop_duplicates()

ci_table["tumour_id"] = tumour_id
ci_table["ci_value"] = ci_value
for allele in ["A", "B"]:
    assert all(
        ci_table[f"cpn{allele}"] >= ci_table[f"lower_CI_{allele}"]
    ), f"cpn{allele} >= lower_CI_{allele}"
    assert all(
        ci_table[f"cpn{allele}"] <= ci_table[f"upper_CI_{allele}"]
    ), f"cpn{allele} <= upper_CI_{allele}"

# ensure all input segments are present in output:
input_segments = refphase_segments["segment"].unique()
output_segments = ci_table["segment"].unique()
missing_segments = set(input_segments) - set(output_segments)
assert (
    len(missing_segments) == 0
), f"Some input segments are missing in the output: {missing_segments}"

# keep only samples present in CONIPHER cp_table:
ci_table = ci_table[ci_table["sample"].isin(conipher_samples)]
ci_table = get_consensus_segmentation(ci_table)
alpaca_input = ci_table.copy()
ci_table.drop(columns=["cn_a", "cn_b", "cpnA", "cpnB", "was_cn_updated"], inplace=True)
ci_table.to_csv(f"{output_dir}/ci_table.csv", index=False)
print(f"{tumour_id} done")

# keep only relevant columns:
print(f"Creating ALPACA input table for {tumour_id}")
alpaca_input = alpaca_input[
    [
        "tumour_id",
        "sample",
        "segment",
        "cpnA",
        "cpnB",
        *(["is_clonal"] if "is_clonal" in alpaca_input.columns else []),
    ]
]
# write to file:
alpaca_input.to_csv(f"{output_dir}/ALPACA_input_table.csv", index=False)

# split input into separate files for each segment to faciliate parallel processing:
if split_segments:
    output_dir_segments = f"{output_dir}/segments"
    os.makedirs(output_dir_segments, exist_ok=True)
    for segment in alpaca_input["segment"].unique():
        alpaca_input[alpaca_input["segment"] == segment].to_csv(
            f"{output_dir_segments}/ALPACA_input_table_{tumour_id}_{segment}.csv",
            index=False,
        )

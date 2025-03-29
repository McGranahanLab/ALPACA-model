import pandas as pd
import argparse
import os
from functions import calculate_confidence_intervals

# arguments
parser = argparse.ArgumentParser(
    description="Calculate confidence intervals from refphase output"
)
parser.add_argument(
    "--tumour_id", type=str, help="Unique identifier for the tumour", required=True
)
parser.add_argument("--output_dir", type=str, help="Output directory", required=True)
parser.add_argument(
    "--refphase_segments",
    type=str,
    help="Location of refphase segments file",
    required=True,
)
parser.add_argument(
    "--refphase_snps", type=str, help="Location of refphase snps file", required=True
)
parser.add_argument(
    "--refphase_purity_ploidy",
    type=str,
    help="Location of refphase purity ploidy file",
    required=True,
)

# options
parser.add_argument(
    "--heterozygous_SNPs_threshold",
    type=int,
    default=5,
    help="Minimum number of heterozygous SNPs to consider a segment. Segments with fewer heterozygous SNPs will be discarded.",
)
parser.add_argument(
    "--ci_value", type=float, default=0.5, help="Confidence interval value."
)
parser.add_argument(
    "--n_bootstrap", type=int, default=100, help="Number of bootstrap samples."
)
parser.add_argument(
    "--split_segments", type=bool, default=False, help="Split input into separate files for each segment. Useful for parallel processing."
)


args = parser.parse_args()
tumour_id = args.tumour_id
output_dir = args.output_dir
ci_value = args.ci_value
n_bootstrap = args.n_bootstrap
split_segments = args.split_segments
# create output directory:
os.makedirs(output_dir, exist_ok=True)

# read data
refphase_segments = pd.read_csv(args.refphase_segments, sep="\t")
refphase_snps = pd.read_csv(args.refphase_snps, sep="\t")
refphase_purity_ploidy = pd.read_csv(args.refphase_purity_ploidy, sep="\t")


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
refphase_segments = refphase_segments.groupby("segment").filter(
    lambda x: (x["heterozygous_SNP_number"] >= args.heterozygous_SNPs_threshold).all()
)


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
    .apply(calculate_confidence_intervals, ci_value=ci_value, n_bootstrap=n_bootstrap)
    .reset_index()
)

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
    ]
].drop_duplicates()

ci_table["tumour_id"] = tumour_id
ci_table["ci_value"] = ci_value

for allele in ["A", "B"]:
    assert all(
        ci_table[f"cpn{allele}"] >= ci_table[f"lower_CI_{allele}"]
    ), f"cpn{allele} >= lower_CI_{allele}"
    assert all(
        ci_table[f"cpn{allele}"] < ci_table[f"upper_CI_{allele}"]
    ), f"cpn{allele} < upper_CI_{allele}"

alpaca_input = ci_table.copy()
ci_table.drop(columns=["cn_a", "cn_b","cpnA", "cpnB", "was_cn_updated"], inplace=True)
ci_table.to_csv(f"{output_dir}/ci_table.csv", index=False)
print(f"{tumour_id} done")

# keep only relevant columns:
print(f"Creating ALPACA input table for {tumour_id}")
alpaca_input = alpaca_input[["tumour_id", "sample", "segment", "cpnA", "cpnB"]]
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

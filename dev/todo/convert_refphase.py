import pandas as pd
import argparse
import os
from functions import calculate_confidence_intervals


    """TODO
    compare recalculated fractinal copy number with the orignal one - the difference seems substantial in some cases
    for example:
    LTX0530:
    segment	sample	cn_a	cn_b	cpnA	cpnB	lower_CI_A	upper_CI_A	lower_CI_B	upper_CI_B	was_cn_updated	tumour_id	ci_value
    10_1411078_1417878	sample_1	3.83696176246305	2.94485073828175	3.5882281799743800	2.753486417616070	3.2097693907961200	3.987499336506610	2.4623174151579900	3.060667510200020	TRUE	LTX0530	0.5
    10_1411078_1417878	sample_2	7.0	2.0	4.843061397927750	2.2535936523772300	4.278604415383480	5.2191565066586100	1.9856149329683400	2.4321466488506800	FALSE	LTX0530	0.5
    
    also, check the segmetns which are neither balanced nor imbalanced:
    
    def calculate_confidence_intervals(seg_sample_df, ci_value=0.95, n_bootstrap=1000):
    # balanced = (seg_sample_df["phasing"] == "none").all()
    # imbalanced = (seg_sample_df["phasing"] != "none").all()
    # if balanced + imbalanced == 0:
    # TODO log this error
    # print("Segment is neither balanced nor imbalanced")
    # raise ValueError
    baf_a = seg_sample_df.query('phasing == "a"')["baf"].mean()
    baf_b = seg_sample_df.query('phasing == "b"')["baf"].mean()
    if math.isnan(baf_a) and math.isnan(baf_b):
        baf_a, baf_b = 0.5, 0.5
    if math.isnan(baf_a):
        baf_a = 1 - baf_b
    if math.isnan(baf_b):
        baf_b = 1 - baf_a
    bafs = {"A": baf_a, "B": baf_b}
    cis = {"A": {}, "B": {}}
    cn_frac = {}
    for allele in ["A", "B"]:
        cn_frac[allele] = max(0, calculate_cn(seg_sample_df, bafs[allele]))
        bootstrap_values = []
        for i in range(n_bootstrap):
            bootstrap_sample_df = bootstrap_sample(seg_sample_df, i)
            bootstrap_value = calculate_cn(bootstrap_sample_df, bafs[allele])
            bootstrap_values.append(bootstrap_value)
        # remove nans from the bootstrap values:
        bootstrap_values = [x for x in bootstrap_values if not np.isnan(x)]
        one_tail_percentile_value = (1 - ci_value) / 2 * 100
        lower_bound = np.percentile(bootstrap_values, one_tail_percentile_value)
        upper_bound = np.percentile(bootstrap_values, 100 - one_tail_percentile_value)
        lower_CI = max(lower_bound, 0)
        upper_CI = max(upper_bound, 0.001)
        cis[allele] = {"lower_CI": lower_CI, "upper_CI": upper_CI}

    return pd.DataFrame(
        {
            "cpnA": cn_frac["A"],
            "lower_CI_A": cis["A"]["lower_CI"],
            "upper_CI_A": cis["A"]["upper_CI"],
            "cpnB": cn_frac["B"],
            "lower_CI_B": cis["B"]["lower_CI"],
            "upper_CI_B": cis["B"]["upper_CI"],
        },
        index=[0],
    )
    
    """
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
    "--ci_value", type=float, default=0.5, help="Confidence interval value"
)
parser.add_argument(
    "--n_bootstrap", type=int, default=100, help="Number of bootstrap samples"
)

args = parser.parse_args()
tumour_id = args.tumour_id
output_dir = args.output_dir
ci_value = args.ci_value
n_bootstrap = args.n_bootstrap
# create output directory:
output_dir_segments = f"{output_dir}/segments"
os.makedirs(output_dir_segments, exist_ok=True)

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

# TODO after testing, remove redundant columns"
# ci_table.drop(columns=["cn_a", "cn_b","cpnA", "cpnB", "was_cn_updated"], inplace=True)
ci_table.to_csv(f"{output_dir}/ci_table.csv", index=False)
print(f"{tumour_id} done")

# keep only relevant columns:
print(f"Creating ALPACA input table for {tumour_id}")
alpaca_input = ci_table[["tumour_id", "sample", "segment", "cpnA", "cpnB"]].copy()
# write to file:

alpaca_input.to_csv(f"{output_dir}/ALPACA_input_table.csv", index=False)

# split input into separate files for each segment to faciliate parallel processing:
for segment in alpaca_input["segment"].unique():
    alpaca_input[alpaca_input["segment"] == segment].to_csv(
        f"{output_dir_segments}/ALPACA_input_table_{tumour_id}_{segment}.csv",
        index=False,
    )

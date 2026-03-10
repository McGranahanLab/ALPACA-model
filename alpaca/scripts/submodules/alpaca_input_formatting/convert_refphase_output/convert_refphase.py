import argparse
import os
import re
from functions import calculate_confidence_intervals, DebugReporter
import pandas as pd


# arguments:
def parse_args():
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
    parser.add_argument(
        "--conipher_cp_table",
        type=str,
        help="Path to CONIPHER cp_table CSV (required)",
        required=True,
    )
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
    parser.add_argument(
        "--debug",
        action="store_true",
        help="If set, writes a detailed per-step debug report with dataframe previews.",
    )
    return parser.parse_args()


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


def main():
    args = parse_args()
    tumour_id = args.tumour_id
    output_dir = args.output_dir
    ci_value = args.ci_value
    n_bootstrap = args.n_bootstrap
    recalculate_not_updated_cns = bool(args.recalculate_not_updated_cns)
    recalculate_updated_cns = bool(args.recalculate_updated_cns)
    recalculate_reference_cns = bool(args.recalculate_reference_cns)
    split_segments = bool(args.split_segments)

    os.makedirs(output_dir, exist_ok=True)
    debug_reporter = DebugReporter(enabled=args.debug, output_dir=output_dir)
    if args.debug:
        print(f"Debug mode enabled. Writing detailed report to {debug_reporter.path}")

    try:
        debug_reporter.section(
            "Step 1: Parse arguments",
            "Capture all options passed to convert_refphase.py before data loading starts.",
        )
        for key, value in vars(args).items():
            debug_reporter.note(f"{key}: {value}")

        debug_reporter.section(
            "Step 2: Load input dataframes",
            "Load Refphase segments/SNPs/purity-ploidy and CONIPHER cp_table from disk.",
        )
        refphase_segments = pd.read_csv(args.refphase_segments, sep="\t")
        refphase_snps = pd.read_csv(args.refphase_snps, sep="\t")
        refphase_purity_ploidy = pd.read_csv(args.refphase_purity_ploidy, sep="\t")
        cp_table = pd.read_csv(args.conipher_cp_table, index_col="clone")
        conipher_samples = cp_table.columns
        debug_reporter.dataframe("refphase_segments (raw)", refphase_segments)
        debug_reporter.dataframe("refphase_snps (raw)", refphase_snps)
        debug_reporter.dataframe("refphase_purity_ploidy (raw)", refphase_purity_ploidy)
        debug_reporter.dataframe("cp_table (raw)", cp_table)
        debug_reporter.note(
            f"Number of CONIPHER samples retained for final filtering: {len(conipher_samples)}"
        )

        debug_reporter.section(
            "Step 3: Standardize key column names",
            "Rename refphase column names to the standard fields used by downstream joins.",
        )
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
        debug_reporter.dataframe("refphase_segments (renamed)", refphase_segments)
        debug_reporter.dataframe("refphase_snps (renamed)", refphase_snps)

        debug_reporter.section(
            "Step 4: Sanitize chromosome naming",
            "Convert chromosome labels to numeric coding shared across inputs.",
        )
        _sanitize_chr_names(refphase_segments)
        _sanitize_chr_names(refphase_snps)
        debug_reporter.dataframe("refphase_segments (chromosomes sanitized)", refphase_segments)
        debug_reporter.dataframe("refphase_snps (chromosomes sanitized)", refphase_snps)

        debug_reporter.section(
            "Step 5: Build segment IDs and filter low-support segments",
            "Create segment IDs from chr/start/end and remove segments below the heterozygous SNP threshold.",
        )
        refphase_segments["segment"] = (
            refphase_segments["chr"].astype(str)
            + "_"
            + refphase_segments["start"].astype(str)
            + "_"
            + refphase_segments["end"].astype(str)
        )
        debug_reporter.dataframe("refphase_segments (with segment id)", refphase_segments)
        segments_before_filter = refphase_segments["segment"].nunique()
        refphase_segments = refphase_segments.groupby("segment").filter(
            lambda x: (x["heterozygous_SNP_number"] >= args.heterozygous_SNPs_threshold).all()
        )
        segments_after_filter = refphase_segments["segment"].nunique()
        debug_reporter.note(
            f"Segments before threshold filtering: {segments_before_filter}; after filtering: {segments_after_filter}; threshold: {args.heterozygous_SNPs_threshold}"
        )
        debug_reporter.dataframe("refphase_segments (post threshold filter)", refphase_segments)

        debug_reporter.section(
            "Step 6: Assign SNPs to segments",
            "Join SNPs to segments by sample and chromosome, then keep SNPs within segment start/end boundaries.",
        )
        snps_with_segments = refphase_snps.merge(
            refphase_segments,
            left_on=["sample", "chr"],
            right_on=["sample", "chr"],
            how="inner",
        )
        debug_reporter.dataframe("snps_with_segments (merged)", snps_with_segments)
        snps_with_segments = snps_with_segments[
            (snps_with_segments["pos"] >= snps_with_segments["start"])
            & (snps_with_segments["pos"] <= snps_with_segments["end"])
        ]
        debug_reporter.dataframe("snps_with_segments (position filtered)", snps_with_segments)

        debug_reporter.section(
            "Step 7: Add purity/ploidy and compute confidence intervals",
            "Merge per-sample purity/ploidy metadata, then run bootstrap CI estimation per segment/sample pair.",
        )
        snps_with_segments_purity_ploidy = snps_with_segments.merge(
            refphase_purity_ploidy, left_on="sample", right_on="sample_id", how="inner"
        )
        debug_reporter.dataframe(
            "snps_with_segments_purity_ploidy", snps_with_segments_purity_ploidy
        )

        print(f"Calculating confidence intervals for {tumour_id}")
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
        debug_reporter.dataframe("confidence_intervals (initial)", confidence_intervals)

        debug_reporter.section(
            "Step 8: Handle segments with zero heterozygous SNPs",
            "If threshold is zero, generate fallback CI bounds for zero-SNP segments and append them.",
        )
        if args.heterozygous_SNPs_threshold == 0:
            ci_span = 0.5
            zero_snp_segments = refphase_segments[refphase_segments.heterozygous_SNP_number == 0]
            zero_snp_segments = zero_snp_segments[["segment", "sample", "cn_a", "cn_b"]]
            ci_half = ci_span / 2.0
            zero_snp_segments["lower_CI_A"] = (zero_snp_segments["cn_a"] - ci_half).clip(lower=0)
            zero_snp_segments["upper_CI_A"] = zero_snp_segments["cn_a"] + ci_half
            zero_snp_segments["lower_CI_B"] = (zero_snp_segments["cn_b"] - ci_half).clip(lower=0)
            zero_snp_segments["upper_CI_B"] = zero_snp_segments["cn_b"] + ci_half
            zero_snp_segments.rename(columns={"cn_a": "cpnA", "cn_b": "cpnB"}, inplace=True)
            debug_reporter.dataframe("zero_snp_segments (generated fallback rows)", zero_snp_segments)
            confidence_intervals = pd.concat(
                [confidence_intervals, zero_snp_segments], ignore_index=True
            )
            debug_reporter.note(f"Appended {len(zero_snp_segments)} zero-SNP segment rows.")
        else:
            debug_reporter.note(
                "Skipped zero-SNP fallback because heterozygous_SNPs_threshold is not 0."
            )
        debug_reporter.dataframe("confidence_intervals (after zero-SNP step)", confidence_intervals)

        debug_reporter.section(
            "Step 9: Build CI table and validate",
            "Sort confidence intervals, merge back segment metadata, and validate CI bounds.",
        )
        confidence_intervals["chr"] = confidence_intervals["segment"].apply(
            lambda x: int(x.split("_")[0])
        )
        confidence_intervals["start"] = confidence_intervals["segment"].apply(
            lambda x: int(x.split("_")[1])
        )
        confidence_intervals = confidence_intervals.sort_values(by=["sample", "chr", "start"])
        confidence_intervals.drop(columns=["chr", "start"], inplace=True)
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
        debug_reporter.dataframe("ci_table (pre-validation)", ci_table)
        for allele in ["A", "B"]:
            assert all(
                ci_table[f"cpn{allele}"] >= ci_table[f"lower_CI_{allele}"]
            ), f"cpn{allele} >= lower_CI_{allele}"
            assert all(
                ci_table[f"cpn{allele}"] <= ci_table[f"upper_CI_{allele}"]
            ), f"cpn{allele} <= upper_CI_{allele}"
        debug_reporter.note("CI bound assertions passed for both alleles.")

        input_segments = refphase_segments["segment"].unique()
        output_segments = ci_table["segment"].unique()
        missing_segments = set(input_segments) - set(output_segments)
        assert (
            len(missing_segments) == 0
        ), f"Some input segments are missing in the output: {missing_segments}"
        debug_reporter.note("All input segments are present in ci_table output.")

        debug_reporter.section(
            "Step 10: Filter samples and write final output files",
            "Keep only samples present in CONIPHER cp_table, then write ci_table.csv and ALPACA_input_table.csv.",
        )
        ci_table = ci_table[ci_table["sample"].isin(conipher_samples)]
        debug_reporter.dataframe("ci_table (post CONIPHER sample filter)", ci_table)

        alpaca_input = ci_table.copy()
        ci_table.drop(columns=["cn_a", "cn_b", "cpnA", "cpnB", "was_cn_updated"], inplace=True)
        debug_reporter.dataframe("ci_table (final columns)", ci_table)
        ci_table.to_csv(f"{output_dir}/ci_table.csv", index=False)
        print(f"{tumour_id} done")

        print(f"Creating ALPACA input table for {tumour_id}")
        alpaca_input = alpaca_input[["tumour_id", "sample", "segment", "cpnA", "cpnB"]]
        debug_reporter.dataframe("alpaca_input (final)", alpaca_input)
        alpaca_input.to_csv(f"{output_dir}/ALPACA_input_table.csv", index=False)

        debug_reporter.section(
            "Step 11: Optional segment splitting",
            "Optionally split ALPACA input into one file per segment for parallel downstream processing.",
        )
        if split_segments:
            output_dir_segments = f"{output_dir}/segments"
            os.makedirs(output_dir_segments, exist_ok=True)
            for segment in alpaca_input["segment"].unique():
                alpaca_input[alpaca_input["segment"] == segment].to_csv(
                    f"{output_dir_segments}/ALPACA_input_table_{tumour_id}_{segment}.csv",
                    index=False,
                )
            debug_reporter.note(
                f"Created {alpaca_input['segment'].nunique()} per-segment files in {output_dir_segments}."
            )
        else:
            debug_reporter.note("Segment splitting disabled (split_segments=0).")

        if args.debug:
            debug_reporter.note("\nDebug report completed successfully.")
    finally:
        debug_reporter.close()


if __name__ == "__main__":
    main()

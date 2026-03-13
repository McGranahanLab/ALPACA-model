import pandas as pd
import numpy as np
import math


def estimate_cn_ascat(baf, logr, purity, ploidy, logr_compaction=1.0):
    cn = (
        purity
        - 1
        + baf * 2 ** (logr / logr_compaction) * ((1 - purity) * 2 + purity * ploidy)
    ) / purity
    return cn


def bootstrap_sample(data, i):
    return data.sample(frac=1, replace=True, random_state=i)


def calculate_final_value_cn_tot(seg_sample_df, logr_shift=0, logr_scale=1):
    baf = 1
    mean_logr = seg_sample_df["logr"].mean()
    purity = seg_sample_df["purity"].unique()[0]
    ploidy = seg_sample_df["ploidy"].unique()[0]
    final_value = estimate_cn_ascat(
        baf, logr_shift + (mean_logr * logr_scale), purity, ploidy
    )
    return final_value


def calculate_confidence_intervals_logr(seg_sample_df, ci_value, n_bootstrap):
    cn_tot = seg_sample_df["cn_a"] + seg_sample_df["cn_b"]
    assert len(cn_tot.unique()) == 1
    cn_tot = cn_tot.unique()[0]
    bootstrap_values = []
    for i in range(n_bootstrap):
        bootstrap_sample_df = bootstrap_sample(seg_sample_df, i)
        bootstrap_value = calculate_final_value_cn_tot(bootstrap_sample_df)
        bootstrap_values.append(bootstrap_value)
    # remove nans from the bootstrap values:
    bootstrap_values = [x for x in bootstrap_values if not np.isnan(x)]
    one_tail_percentile_value = (1 - ci_value) / 2 * 100
    lower_bound = np.percentile(bootstrap_values, one_tail_percentile_value)
    upper_bound = np.percentile(bootstrap_values, 100 - one_tail_percentile_value)
    if cn_tot > 0:
        a_frac = seg_sample_df["cn_a"].values[0] / cn_tot
        b_frac = seg_sample_df["cn_b"].values[0] / cn_tot
    else:
        # to avoid division by zero in homozygous deletions:
        a_frac, b_frac = 0.5, 0.5
    # apply the ratio to bounds:
    # prevent negative values
    lower_CI_A = max(lower_bound * a_frac, 0)
    lower_CI_B = max(lower_bound * b_frac, 0)
    # ci range cannot be 0, this ensures that minimum span will be 0.001
    upper_CI_A = max(upper_bound * a_frac, 0.001)
    upper_CI_B = max(upper_bound * b_frac, 0.001)

    return pd.DataFrame(
        {
            "lower_CI_A": lower_CI_A,
            "upper_CI_A": upper_CI_A,
            "lower_CI_B": lower_CI_B,
            "upper_CI_B": upper_CI_B,
        },
        index=[0],
    )


def calculate_cn(seg_sample_df, baf, logr_shift=0, logr_scale=1):
    mean_logr = seg_sample_df["logr"].mean()
    purity = seg_sample_df["purity"].unique()[0]
    ploidy = seg_sample_df["ploidy"].unique()[0]
    final_value = estimate_cn_ascat(
        baf, logr_shift + (mean_logr * logr_scale), purity, ploidy
    )
    return final_value


def recalculate_updated_cns(seg_sample_df):
    """
    Placeholder for refphase reference-segment copy-number recalculation.
    """
    return None


def _estimate_phased_bafs(seg_sample_df):
    baf_a = seg_sample_df.query('phasing == "a"')["baf"].mean()
    baf_b = seg_sample_df.query('phasing == "b"')["baf"].mean()
    if math.isnan(baf_a) and math.isnan(baf_b):
        baf_a, baf_b = 0.5, 0.5
    if math.isnan(baf_a):
        baf_a = 1 - baf_b
    if math.isnan(baf_b):
        baf_b = 1 - baf_a
    return {"A": baf_a, "B": baf_b}


def calculate_confidence_intervals(seg_sample_df, ci_value, n_bootstrap, recalculate_not_updated_cns, recalculate_updated_cns, recalculate_reference_cns):
    cis = {"A": {}, "B": {}}
    cn_frac = {}
    refphase_updated_cns = seg_sample_df.was_cn_updated.unique()[0]
    refphase_reference_segment = seg_sample_df.is_reference.unique()[0]
    bootstrap_values = {"A": [], "B": []}
    for i in range(n_bootstrap):
        bootstrap_sample_df = bootstrap_sample(seg_sample_df, i)
        bafs = _estimate_phased_bafs(bootstrap_sample_df)
        for allele in ["A", "B"]:
            bootstrap_value = calculate_cn(bootstrap_sample_df, bafs[allele])
            if not np.isnan(bootstrap_value):
                bootstrap_values[allele].append(bootstrap_value)
    for allele in ["A", "B"]:
        one_tail_percentile_value = (1 - ci_value) / 2 * 100
        lower_bound = np.percentile(
            bootstrap_values[allele], one_tail_percentile_value
        )
        upper_bound = np.percentile(
            bootstrap_values[allele], 100 - one_tail_percentile_value
        )
        lower_CI = max(lower_bound, 0)
        upper_CI = max(upper_bound, 0.001)
        cis[allele] = {"lower_CI": lower_CI, "upper_CI": upper_CI}
        # We belive that current implementation of Refphase contains a bug and segments marked as 'is_reference' should also be marked as 'was_cn_updated'.
        # Therefore, we include an option here to fix this and recalculate reference segment copy number without rounding:
        if refphase_reference_segment:
            if recalculate_reference_cns:
                cn_frac[allele] = (lower_CI + upper_CI) / 2
                continue
        if refphase_updated_cns:
            # if segment was updated by refphase
            if recalculate_updated_cns:
                # but we want to recalculate the copy numbers based on BAF and LogR
                # keep the calculated CIs and set fractional copy number as mean of the intervals
                # alternatively, recalculate the copy number, but with few bootstraps it might fall outsie the cofidence intervals on some occasions:
                # cn_frac[allele] = max(0, calculate_cn(seg_sample_df, bafs[allele]))
                cn_frac[allele] = (lower_CI + upper_CI) / 2
            else:
                # but we want to use refphase updated copy number and center the CIs around it:
                half_ci_span = (upper_CI - lower_CI) / 2
                refphase_cns = seg_sample_df[f'cn_{allele.lower()}'].unique()[0]  # value updated by refphase
                lower_CI = refphase_cns - half_ci_span
                upper_CI = refphase_cns + half_ci_span
                cis[allele] = {"lower_CI": lower_CI, "upper_CI": upper_CI}
                cn_frac[allele] = refphase_cns
        else:
            # if segment was NOT updated by refphase
            if recalculate_not_updated_cns:
                # but we want to recalculate the copy numbers based on BAF and LogR:
                # keep the calculated CIs and set fractional copy number as mean of the intervals
                # alternatively, recalculate the copy number, but with few bootstraps it might fall outsie the cofidence intervals on some occasions:
                # cn_frac[allele] = max(0, calculate_cn(seg_sample_df, bafs[allele]))
                cn_frac[allele] = (lower_CI + upper_CI) / 2
            else:
                # but we want to use refphase copy number and center the CIs around it:
                half_ci_span = (upper_CI - lower_CI) / 2
                refphase_cns = seg_sample_df[f'cn_{allele.lower()}'].unique()[0]
                lower_CI = refphase_cns - half_ci_span
                upper_CI = refphase_cns + half_ci_span
                cis[allele] = {"lower_CI": lower_CI, "upper_CI": upper_CI}
                cn_frac[allele] = refphase_cns 
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


def get_consensus_segmentation(ci_table):
    # skip if all samples already have exact same segments:
    common_segments = set(ci_table["segment"].unique())
    if ci_table.groupby("sample")["segment"].apply(lambda x: set(x.unique()) == common_segments).all():
        return ci_table
    
    original_columns = ci_table.columns.tolist()
    segment_parts = ci_table["segment"].str.split("_", expand=True)
    ci_table = ci_table.copy()
    ci_table["chr"] = segment_parts[0].astype(np.int64)
    ci_table["start"] = segment_parts[1].astype(np.int64)
    ci_table["end"] = segment_parts[2].astype(np.int64)

    coverage_before = (
        (ci_table["end"] - ci_table["start"] + 1)
        .groupby(ci_table["sample"])
        .sum()
        .to_dict()
    )

    samples = ci_table["sample"].drop_duplicates().tolist()
    n_samples = len(samples)
    consensus_parts = []

    for chr_value, chr_df in ci_table.groupby("chr", sort=True):
        print('Processing chromosome:', chr_value)
        chr_df = chr_df.sort_values(["sample", "start", "end"]).reset_index(drop=True)
        sample_chr_groups = {sample: df for sample, df in chr_df.groupby("sample", sort=False)}
        breakpoints = np.unique(
            np.concatenate(
                [
                    chr_df["start"].to_numpy(dtype=np.int64),
                    chr_df["end"].to_numpy(dtype=np.int64) + 1,
                ]
            )
        )
        if breakpoints.size < 2:
            continue
        interval_starts = breakpoints[:-1]
        interval_ends = breakpoints[1:] - 1
        n_intervals = interval_starts.size
        coverage_counts = np.zeros(n_intervals, dtype=np.int32)
        sample_maps = {}
        for sample in samples:
            sample_df = sample_chr_groups.get(sample)
            if sample_df is None or sample_df.empty:
                idx = np.full(n_intervals, -1, dtype=np.int64)
                valid = np.zeros(n_intervals, dtype=bool)
            else:
                sample_df = sample_df.sort_values("start").reset_index(drop=True)
                seg_starts = sample_df["start"].to_numpy(dtype=np.int64)
                seg_ends = sample_df["end"].to_numpy(dtype=np.int64)
                idx = np.searchsorted(seg_starts, interval_starts, side="right") - 1
                valid = idx >= 0
                if valid.any():
                    valid_idx = idx[valid]
                    valid[valid] &= interval_ends[valid] <= seg_ends[valid_idx]
            coverage_counts += valid
            sample_maps[sample] = (sample_df, idx)

        keep_mask = coverage_counts == n_samples
        if not keep_mask.any():
            continue

        kept_starts = interval_starts[keep_mask]
        kept_ends = interval_ends[keep_mask]
        for sample in samples:
            sample_df, idx = sample_maps[sample]
            sample_consensus = sample_df.iloc[idx[keep_mask]].copy()
            sample_consensus["chr"] = chr_value
            sample_consensus["start"] = kept_starts
            sample_consensus["end"] = kept_ends
            sample_consensus["segment"] = (
                sample_consensus["chr"].astype(str)
                + "_"
                + sample_consensus["start"].astype(str)
                + "_"
                + sample_consensus["end"].astype(str)
            )
            consensus_parts.append(sample_consensus)
    if consensus_parts:
        ci_table = pd.concat(consensus_parts, ignore_index=True)
    else:
        ci_table = ci_table.iloc[0:0].copy()

    ci_table = ci_table.sort_values(["chr", "start", "sample"]).reset_index(drop=True)
    ci_table = ci_table.drop_duplicates(subset=["sample", "segment"], keep="first")

    coverage_after = (
        (ci_table["end"] - ci_table["start"] + 1)
        .groupby(ci_table["sample"])
        .sum()
        .to_dict()
    )
    print("Coverage report (bp) before/after consensus segmentation:")
    for sample in samples:
        before_bp = int(coverage_before.get(sample, 0))
        after_bp = int(coverage_after.get(sample, 0))
        dropped_bp = before_bp - after_bp
        retained_pct = (100.0 * after_bp / before_bp) if before_bp > 0 else 0.0
        print(
            f"{sample}: before={before_bp}, after={after_bp}, "
            f"dropped={dropped_bp}, retained={retained_pct:.2f}%"
        )

    ci_table = ci_table.drop(columns=["chr", "start", "end"])
    ci_table = ci_table[original_columns]
    return ci_table


def calibrate_battenberg_cns_and_cis(confidence_intervals, refphase_segments):
    # recalculated confidence intervals and copy numbers for battenberg are sometimes different from the original ones. 
    # To ensure consistency, update copy numbers and ci values to match the ones in the original files and emit an alert if the difference is large.
    confidence_intervals = confidence_intervals.merge(
        refphase_segments[["segment", "sample", "cntot"]],
        on=["segment", "sample"],
        how="left"
    )
    assert confidence_intervals["cntot"].notnull().all(), "cntot values are missing for some segments after merge. Please check the input files and merging keys."
    confidence_intervals['AB'] = confidence_intervals['cpnA'] + confidence_intervals['cpnB']
    confidence_intervals['baf'] = confidence_intervals['cpnB'] / confidence_intervals['AB']
    confidence_intervals['cn_diff'] = confidence_intervals['AB'] - confidence_intervals['cntot']
    confidence_intervals['AB_cal'] = confidence_intervals['AB'] + confidence_intervals['cn_diff']
    confidence_intervals['cpnA_cal'] = confidence_intervals['AB_cal'] * (1 - confidence_intervals['baf'])
    confidence_intervals['cpnB_cal'] = confidence_intervals['AB_cal'] * confidence_intervals['baf']
    confidence_intervals['cn_diff_A'] = confidence_intervals['cpnA_cal'] - confidence_intervals['cpnA']
    confidence_intervals['cn_diff_B'] = confidence_intervals['cpnB_cal'] - confidence_intervals['cpnB']
    confidence_intervals['lower_CI_A_cal'] = confidence_intervals['lower_CI_A'] + confidence_intervals['cn_diff_A']
    confidence_intervals['upper_CI_A_cal'] = confidence_intervals['upper_CI_A'] + confidence_intervals['cn_diff_A']
    confidence_intervals['lower_CI_B_cal'] = confidence_intervals['lower_CI_B'] + confidence_intervals['cn_diff_B']
    confidence_intervals['upper_CI_B_cal'] = confidence_intervals['upper_CI_B'] + confidence_intervals['cn_diff_B']
    diff_A = (confidence_intervals['cpnA_cal'] - confidence_intervals['cpnA']).abs().mean()
    diff_B = (confidence_intervals['cpnB_cal'] - confidence_intervals['cpnB']).abs().mean()
    if max(diff_A, diff_B) > 0.5:
        print(f"Warning: Average difference between original and recalculated copy numbers is {max(diff_A, diff_B):.3f}, which is above the threshold of 0.5. Please check the input files and recalculation logic.")
    
    # ensure lower bound is not negative:
    confidence_intervals['lower_CI_A_cal'] = confidence_intervals['lower_CI_A_cal'].clip(lower=0)
    confidence_intervals['lower_CI_B_cal'] = confidence_intervals['lower_CI_B_cal'].clip(lower=0)
    
    confidence_intervals = confidence_intervals.drop(
        columns=[
            "cntot",
            "AB",
            "baf",
            "cn_diff",
            "AB_cal",
            "cpnA",
            "cpnB",
            "cn_diff_A",
            "cn_diff_B",
            "lower_CI_A",
            "upper_CI_A",
            "lower_CI_B",
            "upper_CI_B",
        ]).rename(columns={
            "cpnA_cal": "cpnA",
            "cpnB_cal": "cpnB",
            "lower_CI_A_cal": "lower_CI_A",
            "upper_CI_A_cal": "upper_CI_A",
            "lower_CI_B_cal": "lower_CI_B",
            "upper_CI_B_cal": "upper_CI_B",
        })
    return confidence_intervals

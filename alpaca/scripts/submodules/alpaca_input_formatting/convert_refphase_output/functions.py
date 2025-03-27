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


def calculate_confidence_intervals_logr(seg_sample_df, ci_value=0.95, n_bootstrap=1000):
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
    ## prevent negative values
    lower_CI_A = max(lower_bound * a_frac, 0)
    lower_CI_B = max(lower_bound * b_frac, 0)
    ## ci range cannot be 0, this ensures that minimum span will be 0.001
    upper_CI_A = max(upper_bound * a_frac, 0.001)
    upper_CI_B = max(upper_bound * b_frac, 0.001)

    return pd.DataFrame(
        {
            f"lower_CI_A": lower_CI_A,
            f"upper_CI_A": upper_CI_A,
            f"lower_CI_B": lower_CI_B,
            f"upper_CI_B": upper_CI_B,
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


def calculate_confidence_intervals(seg_sample_df, ci_value=0.95, n_bootstrap=1000):
    balanced = (seg_sample_df["phasing"] == "none").all()
    imbalanced = (seg_sample_df["phasing"] != "none").all()
    if balanced + imbalanced == 0:
        # TODO log this error
        print("Segment is neither balanced nor imbalanced")
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
            f"cpnA": cn_frac["A"],
            f"lower_CI_A": cis["A"]["lower_CI"],
            f"upper_CI_A": cis["A"]["upper_CI"],
            f"cpnB": cn_frac["B"],
            f"lower_CI_B": cis["B"]["lower_CI"],
            f"upper_CI_B": cis["B"]["upper_CI"],
        },
        index=[0],
    )

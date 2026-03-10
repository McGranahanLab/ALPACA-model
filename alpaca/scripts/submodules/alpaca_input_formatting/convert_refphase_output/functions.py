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


def calculate_confidence_intervals(seg_sample_df, ci_value, n_bootstrap, recalculate_not_updated_cns, recalculate_updated_cns, recalculate_reference_cns):
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
    refphase_updated_cns = seg_sample_df.was_cn_updated.unique()[0]
    refphase_reference_segment = seg_sample_df.is_reference.unique()[0]
    for allele in ["A", "B"]:
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


class DebugReporter:
    def __init__(self, enabled, output_dir):
        self.enabled = enabled
        self.path = None
        self._handle = None
        if self.enabled:
            self.path = os.path.join(output_dir, "convert_refphase_debug_report.txt")
            self._handle = open(self.path, "w", encoding="utf-8")
            self._write(f"convert_refphase debug report\nGenerated: {datetime.utcnow().isoformat()}Z\n")

    def _write(self, text):
        if not self.enabled:
            return
        self._handle.write(f"{text}\n")
        self._handle.flush()

    def section(self, title, description):
        if not self.enabled:
            return
        self._write(f"\n{'=' * 90}\n{title}\n{'=' * 90}")
        self._write(description)

    def note(self, message):
        self._write(message)

    def dataframe(self, name, df):
        if not self.enabled:
            return
        self._write(f"\nDataFrame: {name}")
        if not isinstance(df, pd.DataFrame):
            self._write(f"  Not a pandas DataFrame (found {type(df)}).")
            return
        self._write(f"  shape: {df.shape}")
        self._write(f"  columns: {list(df.columns)}")
        self._write("  head(5):")
        if df.empty:
            self._write("  <empty dataframe>")
        else:
            self._write(df.head(5).to_string(index=True))

    def close(self):
        if self._handle is not None:
            self._handle.close()
            self._handle = None

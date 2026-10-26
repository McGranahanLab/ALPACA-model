"""Detection, solving and reporting for degenerate tumours.

A tumour is degenerate when its clone tree has a single clone, when it was
sequenced as a single sample, or both. Clone and sample counts come from the
clone proportions table (rows are the tree's clones, columns are samples).

- single clone (any number of samples): the ILP and the elbow search are
  skipped; each allele gets the integer that first minimises the number of
  confidence-interval misses and then the summed distance to the observations.
- single sample, several clones: the ILP runs with heterogeneity restricted
  (clones are forced clonal where the sample's CI admits a single integer) and
  each segment is checked for an equally good alternative solution.
"""

import glob
import math
import os
from datetime import datetime
from typing import Dict, Iterable, List, Optional, Tuple

import pandas as pd

from alpaca.utils import count_integers, get_length_from_name, read_dataframe, save_dataframe

STANDARD = "standard"
SINGLE_CLONE = "single_clone"
SINGLE_SAMPLE = "single_sample"
SINGLE_CLONE_AND_SINGLE_SAMPLE = "single_clone_and_single_sample"
EDGE_CASE_CATEGORIES = (SINGLE_CLONE, SINGLE_SAMPLE, SINGLE_CLONE_AND_SINGLE_SAMPLE)

# per-segment flags
NO_COMMON_INTEGER = "no_common_integer"
CI_CONTAINS_NO_INTEGER = "ci_contains_no_integer"
ROUNDING_TIE = "rounding_tie"
AMBIGUOUS = "ambiguous"
AMBIGUITY_CHECK_FAILED = "ambiguity_check_failed"

REPORT_SUFFIX = "_edge_case_report"
REPORT_COLUMNS = [
    "tumour_id",
    "segment",
    "category",
    "n_clones",
    "n_samples",
    "solver_path",
    "restrict_heterogeneity",
    "flags",
    "details",
]

FLAG_DESCRIPTIONS = {
    NO_COMMON_INTEGER: "the samples' confidence intervals share no common integer, "
    "so no single clonal state fits every sample (possible missing subclone, "
    "purity/ploidy misfit or segmentation problem)",
    CI_CONTAINS_NO_INTEGER: "the sample's confidence interval contains no integer, "
    "so a monoclonal sample cannot be fitted exactly (possible missing subclone or "
    "purity/ploidy misfit)",
    ROUNDING_TIE: "two integer states fit the observations equally well; the one "
    "closer to diploid was chosen",
    AMBIGUOUS: "another clone copy-number assignment fits the sample equally well "
    "at the selected complexity",
    AMBIGUITY_CHECK_FAILED: "the ambiguity check could not be completed",
}

CATEGORY_DESCRIPTIONS = {
    SINGLE_CLONE: "The clone tree contains a single clone. Every sample is "
    "monoclonal, so each segment was solved directly (no ILP, no elbow search) by "
    "choosing, per allele, the integer that minimises confidence-interval misses "
    "and then the summed distance to the observed copy numbers.",
    SINGLE_SAMPLE: "The tumour has a single sample but several clones. Subclonal "
    "copy numbers are weakly identifiable from one sample (two observations per "
    "segment for two unknowns per clone), so heterogeneity was restricted where "
    "enabled and each segment was checked for equally good alternative solutions.",
    SINGLE_CLONE_AND_SINGLE_SAMPLE: "The tumour has a single clone and a single "
    "sample. Each segment reduces to rounding the observed copy number; it was "
    "solved directly (no ILP, no elbow search).",
}


def classify_tumour(cp_table: pd.DataFrame) -> str:
    """Classify a tumour from its clone proportions table (clones x samples)."""
    single_clone = len(cp_table.index) == 1
    single_sample = len(cp_table.columns) == 1
    if single_clone and single_sample:
        return SINGLE_CLONE_AND_SINGLE_SAMPLE
    if single_clone:
        return SINGLE_CLONE
    if single_sample:
        return SINGLE_SAMPLE
    return STANDARD


def read_cp_table(tumour_dir: str) -> pd.DataFrame:
    return pd.read_csv(os.path.join(tumour_dir, "cp_table.csv"), index_col="clone")


def classify_tumour_dir(tumour_dir: str) -> str:
    return classify_tumour(read_cp_table(tumour_dir))


# ---------------------------------------------------------------------------
# single-clone direct solution
# ---------------------------------------------------------------------------


def homozygous_deletion_permitted(
    observed: Dict[str, Dict[str, float]], segment: str, model_config: dict
) -> bool:
    """Whether a lone clone (present in every sample) may have CN 0 on both alleles.

    Mirrors the homozygous-deletion constraint of the ILP backends.
    """
    if not model_config.get("limit_homozygous_deletions_threshold_flag", True):
        return True
    threshold = float(model_config.get("homozygous_deletion_threshold", 1))
    strong_threshold = float(model_config.get("homo_del_strong_evidence_threshold", 0.5))
    size_limit = int(model_config.get("homo_del_size_limit", 5 * 10**7))
    samples = observed["A"].keys()
    low = [s for s in samples if observed["A"][s] < threshold and observed["B"][s] < threshold]
    strong = [
        s
        for s in samples
        if observed["A"][s] < strong_threshold and observed["B"][s] < strong_threshold
    ]
    if strong:
        return True
    return get_length_from_name(segment) < size_limit and bool(low)


def _allele_scores(values: List[float], lowers: List[float], uppers: List[float], cn: int):
    ci_misses = sum(1 for lo, up in zip(lowers, uppers) if cn < lo or cn > up)
    distance = sum(abs(cn - v) for v in values)
    return ci_misses, distance


def solve_single_clone(
    input_table: pd.DataFrame,
    ci_table: pd.DataFrame,
    segment: str,
    model_config: dict,
) -> Tuple[Tuple[int, int], Dict[str, float], List[str], str]:
    """Pick integer (A, B) for a single clone directly.

    Candidates are ranked by (CI misses, summed distance) over both alleles and
    then by the number of events from diploid; homozygous deletion is excluded
    unless permitted. Returns ((A, B), scores, flags, details).
    """
    ci = ci_table[ci_table.segment == segment].set_index("sample")
    samples = sorted(input_table["sample"].unique())
    observed = {
        allele: dict(zip(input_table["sample"], input_table[f"cpn{allele}"]))
        for allele in ("A", "B")
    }
    bounds = {
        allele: (
            [float(ci.loc[s, f"lower_CI_{allele}"]) for s in samples],
            [float(ci.loc[s, f"upper_CI_{allele}"]) for s in samples],
        )
        for allele in ("A", "B")
    }
    per_allele = {}
    for allele in ("A", "B"):
        values = [float(observed[allele][s]) for s in samples]
        lowers, uppers = bounds[allele]
        top = int(math.ceil(max(values + uppers))) + 1
        per_allele[allele] = {
            cn: _allele_scores(values, lowers, uppers, cn) for cn in range(0, top + 1)
        }
    homdel_ok = homozygous_deletion_permitted(observed, segment, model_config)
    candidates = []
    for a, (ci_a, d_a) in per_allele["A"].items():
        for b, (ci_b, d_b) in per_allele["B"].items():
            if a == 0 and b == 0 and not homdel_ok:
                continue
            fit = (ci_a + ci_b, round(d_a + d_b, 6))
            events = abs(a - 1) + abs(b - 1)
            candidates.append((fit, events, a, b))
    candidates.sort()
    best_fit, events, a, b = candidates[0]
    flags = []
    details = []
    for allele, cn in (("A", a), ("B", b)):
        lowers, uppers = bounds[allele]
        if len(samples) > 1:
            if max(math.ceil(lo) for lo in lowers) > min(math.floor(up) for up in uppers):
                flags.append(f"{NO_COMMON_INTEGER}_{allele}")
        elif count_integers(lowers[0], uppers[0]) == 0:
            flags.append(f"{CI_CONTAINS_NO_INTEGER}_{allele}")
        fits = {c: (miss, round(dist, 6)) for c, (miss, dist) in per_allele[allele].items()}
        best = min(fits.values())
        tied = [c for c, fit in fits.items() if fit == best]
        chosen = f"chose {cn}"
        if len(tied) > 1:
            flags.append(f"{ROUNDING_TIE}_{allele}")
            chosen += f" (tied states {tied})"
        obs = ", ".join(
            f"{s}={observed[allele][s]:.2f} [{lo:.2f}, {up:.2f}]"
            for s, lo, up in zip(samples, lowers, uppers)
        )
        details.append(f"{allele}: observed {obs}; {chosen}")
    scores = {"CI_score": best_fit[0], "D_score": round(best_fit[1], 3), "complexity": events}
    return (a, b), scores, flags, "; ".join(details)


# ---------------------------------------------------------------------------
# per-segment and tumour-level reports
# ---------------------------------------------------------------------------


def make_segment_report(
    tumour_id: str,
    segment: str,
    category: str,
    n_clones: int,
    n_samples: int,
    solver_path: str,
    restrict_heterogeneity: bool,
    flags: Iterable[str],
    details: str = "",
) -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "tumour_id": tumour_id,
                "segment": segment,
                "category": category,
                "n_clones": n_clones,
                "n_samples": n_samples,
                "solver_path": solver_path,
                "restrict_heterogeneity": bool(restrict_heterogeneity),
                "flags": ";".join(flags),
                "details": details,
            }
        ],
        columns=REPORT_COLUMNS,
    )


def segment_report_name(tumour_id: str, segment: str, fmt: str) -> str:
    return f"{tumour_id}_{segment}{REPORT_SUFFIX}.{fmt}"


def find_segment_reports(reports_dir: str, tumour_id: Optional[str] = None) -> List[str]:
    prefix = f"{tumour_id}_" if tumour_id else ""
    files = []
    for ext in ("csv", "parquet"):
        files += glob.glob(os.path.join(reports_dir, f"{prefix}*{REPORT_SUFFIX}.{ext}"))
    return sorted(files)


def combine_segment_reports(files: Iterable[str]) -> pd.DataFrame:
    frames = [read_dataframe(f) for f in files]
    frames = [f for f in frames if not f.empty]
    if not frames:
        return pd.DataFrame(columns=REPORT_COLUMNS)
    combined = pd.concat(frames, ignore_index=True)
    combined["flags"] = combined["flags"].fillna("").astype(str)
    combined["details"] = combined["details"].fillna("").astype(str)
    return combined.sort_values(["tumour_id", "segment"]).reset_index(drop=True)


def _flag_base(flag: str) -> str:
    """Strip the allele suffix from per-allele flags, e.g. rounding_tie_A."""
    return flag[:-2] if flag.endswith(("_A", "_B")) else flag


def warning_file_name(category: str) -> str:
    return f"warning_{category}.txt"


def render_warning(
    tumour_id: str,
    category: str,
    n_clones: int,
    n_samples: int,
    segment_reports: pd.DataFrame,
    cohort_status: Optional[str] = None,
) -> str:
    lines = [
        f"WARNING: tumour {tumour_id} is a {category.replace('_', ' ')} tumour",
        f"clones: {n_clones}, samples: {n_samples}",
        f"generated: {datetime.now().isoformat(timespec='seconds')}",
        "",
        CATEGORY_DESCRIPTIONS[category],
    ]
    if cohort_status:
        lines += ["", f"Cohort status: {cohort_status}"]
    flagged = segment_reports[segment_reports["flags"] != ""]
    lines += [
        "",
        f"Segments solved: {len(segment_reports)}",
        f"Segments flagged: {len(flagged)}",
    ]
    if not flagged.empty:
        used = sorted({_flag_base(f) for fl in flagged["flags"] for f in fl.split(";") if f})
        lines += ["", "Flags:"]
        lines += [f"  {f}: {FLAG_DESCRIPTIONS.get(f, '')}" for f in used]
        lines += ["", "Flagged segments:"]
        for _, row in flagged.iterrows():
            lines.append(f"  {row['segment']}\t{row['flags']}")
            if row["details"]:
                lines.append(f"      {row['details']}")
    return "\n".join(lines) + "\n"


def load_segment_reports(reports: str, tumour_id: Optional[str] = None) -> pd.DataFrame:
    """Read segment reports from a directory of per-segment files or one combined table."""
    if os.path.isdir(reports):
        combined = combine_segment_reports(find_segment_reports(reports, tumour_id))
    elif os.path.isfile(reports):
        combined = combine_segment_reports([reports])
    else:
        combined = pd.DataFrame(columns=REPORT_COLUMNS)
    if tumour_id is not None:
        combined = combined[combined["tumour_id"].astype(str) == str(tumour_id)]
    return combined.reset_index(drop=True)


def write_tumour_edge_case_outputs(
    tumour_dir: str,
    reports: str,
    output_dir: str,
    tumour_id: Optional[str] = None,
    output_format: str = "csv",
    cohort_status: Optional[str] = None,
) -> Optional[str]:
    """Write warning_<category>.txt and the combined segment report for one tumour.

    ``reports`` is a directory of per-segment reports or one combined report table.
    Returns the category, or None for a standard tumour (nothing is written).
    """
    cp_table = read_cp_table(tumour_dir)
    category = classify_tumour(cp_table)
    if category == STANDARD:
        return None
    tumour_id = tumour_id or os.path.basename(os.path.normpath(tumour_dir))
    segment_reports = load_segment_reports(reports, tumour_id)
    os.makedirs(output_dir, exist_ok=True)
    text = render_warning(
        tumour_id,
        category,
        len(cp_table.index),
        len(cp_table.columns),
        segment_reports,
        cohort_status=cohort_status,
    )
    with open(os.path.join(output_dir, warning_file_name(category)), "w") as fh:
        fh.write(text)
    save_dataframe(
        segment_reports,
        output_dir,
        f"edge_case_report.{output_format}",
        output_format=output_format,
    )
    return category

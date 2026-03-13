#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path
import re
from typing import Dict, Iterable, List

import numpy as np
import pandas as pd

BATTENBERG_PLOIDY_SOURCE_COLUMN = "psi"
BATTENBERG_SOLUTION_ID = 'A'


def _sanitize_chr_series(chr_series: pd.Series) -> pd.Series:
    if pd.api.types.is_numeric_dtype(chr_series):
        return chr_series
    out = chr_series.astype(str).str.extract(
        r"([0-9]+|[XYM][Tt]?)", flags=re.IGNORECASE
    )[0]
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
    out = out.replace(chr_map)
    return pd.to_numeric(out, errors="coerce")


def _normalize_chr_value(chr_value: str) -> int:
    extracted = re.findall(r"([0-9]+|[XYM][Tt]?)", str(chr_value), flags=re.IGNORECASE)
    if not extracted:
        raise ValueError(f"Could not parse chromosome value: {chr_value}")
    token = extracted[0].upper()
    token = {"X": "23", "Y": "24", "MT": "25", "M": "25"}.get(token, token)
    return int(token)


def _extract_sample_from_logr_segmented(path: Path) -> str:
    name = path.name
    m = re.match(r"^(?P<sample>.+)\.logRsegmented\.txt(?:\.gz)?$", name)
    if m:
        return m.group("sample")
    if ".logRsegmented" in name:
        return name.split(".logRsegmented", 1)[0]
    raise ValueError(f"Could not derive sample name from file name: {path}")


def _resolve_path(path_value: str, base_dir: Path) -> Path:
    path = Path(path_value).expanduser()
    if not path.is_absolute():
        path = (base_dir / path).resolve()
    return path


def _choose_path(
    candidates: Iterable[Path],
    *,
    label: str,
    sample: str,
    preferred_substring: str | None = None,
) -> Path:
    candidates_sorted = sorted(set(candidates))
    if not candidates_sorted:
        raise FileNotFoundError(
            f"No files found for sample '{sample}' ({label})."
        )
    if preferred_substring:
        preferred = [
            p
            for p in candidates_sorted
            if preferred_substring.lower() in p.name.lower()
        ]
        if preferred:
            return preferred[0]
    return candidates_sorted[0]


def _find_col(df: pd.DataFrame, candidates: List[str], required: bool = True) -> str | None:
    lower_to_orig = {str(col).lower(): col for col in df.columns}
    for candidate in candidates:
        if candidate.lower() in lower_to_orig:
            return lower_to_orig[candidate.lower()]
    if required:
        raise KeyError(
            f"Required column not found. Expected one of {candidates}, found {list(df.columns)}"
        )
    return None


def _read_table_with_optional_header(
    path: Path,
    *,
    sep: str,
    required_col_groups: List[List[str]],
    fallback_first_columns: List[str],
    min_columns: int = 3,
) -> pd.DataFrame:
    """
    Read a table that may be missing a header line.
    If required columns are not detected, assume headerless format where:
    col1=chromosome, col2=position, col3=value (or equivalent).
    """
    df = pd.read_csv(path, sep=sep, compression="infer")
    has_required_columns = all(
        _find_col(df, candidates, required=False) is not None
        for candidates in required_col_groups
    )
    if has_required_columns:
        return df

    # Fall back to headerless interpretation.
    df_no_header = pd.read_csv(path, sep=sep, compression="infer", header=None)
    if df_no_header.shape[1] < min_columns:
        raise ValueError(
            f"Header not detected and table has too few columns in {path}. "
            f"Expected at least {min_columns}, found {df_no_header.shape[1]}."
        )
    colnames = list(fallback_first_columns)
    if len(colnames) > df_no_header.shape[1]:
        raise ValueError(
            f"Fallback column definition ({len(colnames)}) exceeds table width "
            f"({df_no_header.shape[1]}) for {path}."
        )
    colnames.extend(
        [f"extra_col_{i}" for i in range(len(colnames), df_no_header.shape[1])]
    )
    df_no_header.columns = colnames
    return df_no_header


def _read_inventory(path: Path, tumour_id: str | None = None) -> List[Dict[str, Path | str]]:
    inventory = pd.read_csv(path, sep=None, engine="python")
    base_dir = path.parent

    tumour_col = _find_col(
        inventory, ["tumour_id", "tumor_id", "case_id"], required=False
    )
    if tumour_col is not None and tumour_id:
        inventory = inventory[
            inventory[tumour_col].astype(str) == str(tumour_id)
        ].copy()
        if inventory.empty:
            raise ValueError(
                f"No rows found in inventory for tumour_id='{tumour_id}'."
            )

    sample_col = _find_col(inventory, ["sample", "sample_name"], required=False)
    logr_col = _find_col(
        inventory,
        ["logr_segmented_path", "logr_segmented", "logRsegmented_path", "logRsegmented"],
    )
    mutant_col = _find_col(
        inventory,
        [
            "mutant_logr_path",
            "mutant_logr",
            "mutantLogR_path",
            "mutantLogR_gcCorrected_path",
        ],
    )
    het_col = _find_col(
        inventory,
        [
            "baf_segmented_path",
            "BAFsegmented_path",
            "heterozygous_baf_path",
            "heterozygous_baf",
            "heterozygousMutBAFs_path",
            "heterozygousMutBAFs_haplotyped_path",
        ],
    )
    purity_col = _find_col(
        inventory,
        ["purity_ploidy_path", "purity_ploidy", "battenberg_purity_ploidy_path"],
    )
    subclones_col = _find_col(
        inventory,
        [
            "subclones_path",
            "subclones_path",
            "battenberg_subclones_path",
            "battenberg_subclones_path",
            "battenberg_solution_path",
        ],
        required=False,
    )

    required_cols = [logr_col, mutant_col, het_col, purity_col]
    for col in required_cols:
        if inventory[col].isnull().any():
            raise ValueError(f"Inventory column '{col}' contains empty values.")

    resolved_rows: List[Dict[str, Path | str]] = []
    for _, row in inventory.iterrows():
        logr_path = _resolve_path(str(row[logr_col]), base_dir)
        sample = (
            str(row[sample_col]).strip()
            if sample_col is not None and str(row[sample_col]).strip()
            else _extract_sample_from_logr_segmented(logr_path)
        )
        if (
            subclones_col is not None
            and pd.notna(row[subclones_col])
            and str(row[subclones_col]).strip()
        ):
            subclones_path = _resolve_path(str(row[subclones_col]), base_dir)
        else:
            subclones_path = _choose_path(
                base_dir.rglob(f"{sample}*subclones*.txt*"),
                label="default Battenberg subclones",
                sample=sample,
                preferred_substring="default",
            )
        resolved_rows.append(
            {
                "sample": sample,
                "logr_segmented_path": logr_path,
                "mutant_logr_path": _resolve_path(str(row[mutant_col]), base_dir),
                "heterozygous_baf_path": _resolve_path(str(row[het_col]), base_dir),
                "purity_ploidy_path": _resolve_path(str(row[purity_col]), base_dir),
                "subclones_path": subclones_path,
            }
        )
    return resolved_rows


def _discover_inventory(input_dir: Path) -> List[Dict[str, Path | str]]:
    sample_rows: List[Dict[str, Path | str]] = []
    logr_segmented_files = sorted(input_dir.rglob("*.logRsegmented.txt*"))
    if not logr_segmented_files:
        raise FileNotFoundError(
            f"No Battenberg logRsegmented files found in {input_dir}."
        )
    seen_samples = set()
    for logr_segmented_path in logr_segmented_files:
        sample = _extract_sample_from_logr_segmented(logr_segmented_path)
        if sample in seen_samples:
            raise ValueError(
                f"Sample '{sample}' resolved from multiple logRsegmented files."
            )
        seen_samples.add(sample)

        mutant_candidates = input_dir.rglob(f"{sample}_mutantLogR_gcCorrected.tab*")
        het_candidates = input_dir.rglob(f"{sample}*BAFsegmented.txt*")
        purity_candidates = input_dir.rglob(f"{sample}_battenbergA*purity_ploidy.txt*")
        subclones_candidates = input_dir.rglob(f"{sample}*subclones*.txt*")

        sample_rows.append(
            {
                "sample": sample,
                "logr_segmented_path": logr_segmented_path,
                "mutant_logr_path": _choose_path(
                    mutant_candidates, label="mutant LogR", sample=sample
                ),
                "heterozygous_baf_path": _choose_path(
                    het_candidates, label="heterozygous BAF", sample=sample
                ),
                "purity_ploidy_path": _choose_path(
                    purity_candidates,
                    label="purity/ploidy",
                    sample=sample,
                    preferred_substring="default",
                ),
                "subclones_path": _choose_path(
                    subclones_candidates,
                    label="default Battenberg subclones",
                    sample=sample,
                    preferred_substring="default",
                ),
            }
        )
    return sample_rows


def _read_logr_segmented(path: Path, chromosome_filter: int | None = None) -> pd.DataFrame:
    df = _read_table_with_optional_header(
        path,
        sep="\t",
        required_col_groups=[
            ["chr", "chromosome", "seqnames"],
            ["pos", "position"],
            ["segmented_logr", "segmented_logR"],
        ],
        fallback_first_columns=["chr", "pos", "segmented_logr"],
        min_columns=3,
    )
    chr_col = _find_col(df, ["chr", "chromosome", "seqnames"])
    pos_col = _find_col(df, ["pos", "position"])
    logr_col = _find_col(df, ["segmented_logr", "segmented_logR"])
    out = df[[chr_col, pos_col, logr_col]].rename(
        columns={chr_col: "chr", pos_col: "pos", logr_col: "segmented_logr"}
    )
    out["chr"] = _sanitize_chr_series(out["chr"])
    out["pos"] = pd.to_numeric(out["pos"], errors="coerce")
    out["segmented_logr"] = pd.to_numeric(out["segmented_logr"], errors="coerce")
    out = out.dropna(subset=["chr", "pos", "segmented_logr"]).copy()
    out["chr"] = out["chr"].astype(int)
    out["pos"] = out["pos"].astype(int)
    if chromosome_filter is not None:
        out = out[out["chr"] == chromosome_filter].copy()
    out = out.sort_values(["chr", "pos"]).reset_index(drop=True)
    return out


def _read_mutant_logr(path: Path, sample: str, chromosome_filter: int | None = None) -> pd.DataFrame:
    df = _read_table_with_optional_header(
        path,
        sep="\t",
        required_col_groups=[
            ["chromosome", "chr", "seqnames"],
            ["position", "pos"],
        ],
        fallback_first_columns=["chromosome", "position", "value"],
        min_columns=3,
    )
    chr_col = _find_col(df, ["chromosome", "chr", "seqnames"])
    pos_col = _find_col(df, ["position", "pos"])
    value_candidates = [c for c in df.columns if c not in {chr_col, pos_col}]
    if not value_candidates:
        raise ValueError(f"No mutant LogR value column found in {path}")
    value_col = sample if sample in df.columns else value_candidates[0]
    out = df[[chr_col, pos_col, value_col]].rename(
        columns={chr_col: "chr", pos_col: "pos", value_col: "logr"}
    )
    out["chr"] = _sanitize_chr_series(out["chr"])
    out["pos"] = pd.to_numeric(out["pos"], errors="coerce")
    out["logr"] = pd.to_numeric(out["logr"], errors="coerce")
    out = out.dropna(subset=["chr", "pos", "logr"]).copy()
    out["chr"] = out["chr"].astype(int)
    out["pos"] = out["pos"].astype(int)
    if chromosome_filter is not None:
        out = out[out["chr"] == chromosome_filter].copy()
    return out


def _read_het_baf(path: Path, sample: str, chromosome_filter: int | None = None) -> pd.DataFrame:
    df = _read_table_with_optional_header(
        path,
        sep="\t",
        required_col_groups=[
            ["chromosome", "chr", "seqnames"],
            ["position", "pos"],
        ],
        fallback_first_columns=["chromosome", "position", "value"],
        min_columns=3,
    )
    chr_col = _find_col(df, ["chromosome", "chr", "seqnames"])
    pos_col = _find_col(df, ["position", "pos"])
    # Use the segmented BAF column explicitly to avoid ambiguity with similarly named columns.
    value_col = _find_col(df, ["BAFseg"], required=False)
    if value_col is None:
        # headerless fallback: the 3rd column is loaded as "value"
        value_col = _find_col(df, ["value"], required=False)
    if value_col is None:
        raise ValueError(
            f"Could not find BAF value column in {path}. Expected 'BAFseg' "
            "or third-column fallback ('value') for headerless files."
        )
    out = df[[chr_col, pos_col, value_col]].rename(
        columns={chr_col: "chr", pos_col: "pos", value_col: "baf"}
    )
    out["chr"] = _sanitize_chr_series(out["chr"])
    out["pos"] = pd.to_numeric(out["pos"], errors="coerce")
    out["baf"] = pd.to_numeric(out["baf"], errors="coerce")
    out = out.dropna(subset=["chr", "pos", "baf"]).copy()
    out["chr"] = out["chr"].astype(int)
    out["pos"] = out["pos"].astype(int)
    if chromosome_filter is not None:
        out = out[out["chr"] == chromosome_filter].copy()
    return out


def _read_purity_ploidy(path: Path, sample: str) -> Dict[str, float | str]:
    df = pd.read_csv(path, sep=r"\s+", compression="infer")
    purity_col = _find_col(df, ["purity"])
    _find_col(df, ["ploidy"])
    if BATTENBERG_PLOIDY_SOURCE_COLUMN not in df.columns:
        raise KeyError(
            f"Expected '{BATTENBERG_PLOIDY_SOURCE_COLUMN}' column in {path}, "
            f"found {list(df.columns)}"
        )
    if df.empty:
        raise ValueError(f"Purity/ploidy file is empty: {path}")
    row = df.iloc[0]
    return {
        "sample_id": sample,
        "purity": float(row[purity_col]),
        "ploidy": float(row[BATTENBERG_PLOIDY_SOURCE_COLUMN]),
    }


def _read_subclones(path: Path, chromosome_filter: int | None = None) -> pd.DataFrame:
    df = _read_table_with_optional_header(
        path,
        sep="\t",
        required_col_groups=[
            ["chr", "chromosome", "seqnames"],
            ["startpos", "start"],
            ["endpos", "end"],
            ["cntot", "cn_tot", "ntot"],
            [f"frac1_{BATTENBERG_SOLUTION_ID}"]
        ],
        fallback_first_columns=["chr", "startpos", "endpos", "cntot"],
        min_columns=4,
    )
    chr_col = _find_col(df, ["chr", "chromosome", "seqnames"])
    start_col = _find_col(df, ["startpos", "start"])
    end_col = _find_col(df, ["endpos", "end"])
    cntot_col = _find_col(df, ["cntot", "cn_tot", "ntot"])
    frac_col = _find_col(df, [f"frac1_{BATTENBERG_SOLUTION_ID}"])
    out = df[[chr_col, start_col, end_col, cntot_col, frac_col]].rename(
        columns={chr_col: "chr", start_col: "start", end_col: "end", cntot_col: "cntot"}
    )
    out["chr"] = _sanitize_chr_series(out["chr"])
    out["start"] = pd.to_numeric(out["start"], errors="coerce")
    out["end"] = pd.to_numeric(out["end"], errors="coerce")
    out["cntot"] = pd.to_numeric(out["cntot"], errors="coerce")
    out = out.dropna(subset=["chr", "start", "end", "cntot"]).copy()
    out["chr"] = out["chr"].astype(int)
    out["start"] = out["start"].astype(int)
    out["end"] = out["end"].astype(int)
    if chromosome_filter is not None:
        out = out[out["chr"] == chromosome_filter].copy()
    out = out.drop_duplicates(subset=["chr", "start", "end"], keep="first")
    return out


def _estimate_cn_tot(logr: float, purity: float, ploidy: float) -> float:
    cn_tot = (
        purity
        - 1
        + (2 ** logr) * ((1 - purity) * 2 + purity * ploidy)
    ) / purity
    return max(float(cn_tot), 0.0)


def _count_snps_per_segment(snps: pd.DataFrame, segments: pd.DataFrame) -> pd.Series:
    if snps.empty or segments.empty:
        return pd.Series(dtype="int64")
    joined = snps.merge(
        segments[["segment", "chr", "start", "end"]],
        on="chr",
        how="inner",
    )
    joined = joined[(joined["pos"] >= joined["start"]) & (joined["pos"] <= joined["end"])]
    return joined.groupby("segment").size()


def _segments_from_logr(logr_segmented: pd.DataFrame) -> pd.DataFrame:
    out = logr_segmented.copy()
    # Battenberg segmenting rule: exact segmented_logr value + genomic adjacency.
    out["segment_group"] = (
        (out["chr"] != out["chr"].shift(1))
        | (out["segmented_logr"] != out["segmented_logr"].shift(1))
    ).cumsum()
    segments = (
        out.groupby(["chr", "segment_group"], as_index=False)
        .agg(
            start=("pos", "min"),
            end=("pos", "max"),
            segmented_logr=("segmented_logr", "first"),
            segment_snp_number=("pos", "size"),
        )
        .drop(columns=["segment_group"])
    )
    segments["segment"] = (
        segments["chr"].astype(str)
        + "_"
        + segments["start"].astype(str)
        + "_"
        + segments["end"].astype(str)
    )
    return segments


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Extract intermediate ALPACA tables from Battenberg sample-level files."
    )
    parser.add_argument("--tumour_id", required=False, type=str)
    parser.add_argument(
        "--chromosome",
        required=False,
        type=str,
        help="Optional chromosome filter (e.g. 1, chr1, X). If set, only this chromosome is kept.",
    )
    parser.add_argument("--output_dir", required=True, type=str)
    parser.add_argument("--battenberg_inventory", required=False, type=str)
    parser.add_argument("--battenberg_input_dir", required=False, type=str)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    output_dir = Path(args.output_dir).expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    chromosome_filter = (
        _normalize_chr_value(args.chromosome) if args.chromosome is not None else None
    )

    if not args.battenberg_inventory and not args.battenberg_input_dir:
        raise ValueError(
            "Provide --battenberg_inventory or --battenberg_input_dir."
        )

    if args.battenberg_inventory:
        inventory_rows = _read_inventory(
            Path(args.battenberg_inventory).expanduser().resolve(),
            tumour_id=args.tumour_id,
        )
    else:
        inventory_rows = _discover_inventory(
            Path(args.battenberg_input_dir).expanduser().resolve()
        )

    phased_segments_all = []
    phased_snps_all = []
    purity_ploidy_rows = []
    resolved_inventory_rows = []

    for row in inventory_rows:
        sample = str(row["sample"])
        logr_segmented_path = Path(row["logr_segmented_path"])
        mutant_logr_path = Path(row["mutant_logr_path"])
        heterozygous_baf_path = Path(row["heterozygous_baf_path"])
        purity_ploidy_path = Path(row["purity_ploidy_path"])
        subclones_path = Path(row["subclones_path"])
        for path in [
            logr_segmented_path,
            mutant_logr_path,
            heterozygous_baf_path,
            purity_ploidy_path,
            subclones_path,
        ]:
            if not path.exists():
                raise FileNotFoundError(f"Missing required input file for {sample}: {path}")

        resolved_inventory_rows.append(
            {
                "sample": sample,
                "logr_segmented_path": str(logr_segmented_path),
                "mutant_logr_path": str(mutant_logr_path),
                "heterozygous_baf_path": str(heterozygous_baf_path),
                "purity_ploidy_path": str(purity_ploidy_path),
                "subclones_path": str(subclones_path),
            }
        )

        logr_segmented = _read_logr_segmented(logr_segmented_path, chromosome_filter)
        mutant_logr = _read_mutant_logr(mutant_logr_path, sample, chromosome_filter)
        het_baf = _read_het_baf(heterozygous_baf_path, sample, chromosome_filter)
        subclones = _read_subclones(subclones_path, chromosome_filter)
        purity_ploidy = _read_purity_ploidy(purity_ploidy_path, sample)
        if logr_segmented.empty:
            raise ValueError(
                f"No logRsegmented rows remain for sample {sample} after chromosome filtering."
            )
        if mutant_logr.empty:
            raise ValueError(
                f"No mutant LogR rows remain for sample {sample} after chromosome filtering."
            )
        if het_baf.empty:
            raise ValueError(
                f"No heterozygous BAF rows remain for sample {sample} after chromosome filtering."
            )
        if subclones.empty:
            raise ValueError(
                f"No default Battenberg subclones rows remain for sample {sample} "
                "after chromosome filtering."
            )

        snps = het_baf.merge(mutant_logr, on=["chr", "pos"], how="inner")
        if snps.empty:
            raise ValueError(
                f"No overlapping heterozygous BAF and mutant LogR SNPs found for sample {sample}."
            )
        snps["phasing"] = "b"
        snps["group_name"] = sample
        snps["seqnames"] = snps["chr"]
        snps["strand"] = "*"
        snps["germline_zygosity"] = "het"
        snps = snps[
            [
                "group_name",
                "seqnames",
                "strand",
                "pos",
                "baf",
                "logr",
                "germline_zygosity",
                "phasing",
            ]
        ].copy()

        segments = _segments_from_logr(logr_segmented)
        ## retain original cntot value, needed downstream:
        # keep also clonal/subclonal assignment
        segments = segments.merge(
            subclones,
            on=["chr", "start", "end"],
            how="left",
            validate="m:1",
        )
        # subclones containes merged segments (larger than derived from logr)
        # for such segments find matching parent:
        cols_to_fix = ["cntot", f"frac1_{BATTENBERG_SOLUTION_ID}"]
        rows_to_fix = segments[cols_to_fix].isna().any(axis=1)
        if rows_to_fix.any():
            missing_segments = segments.loc[
                rows_to_fix, ["segment", "chr", "start", "end"]
            ].copy()
            missing_segments["_segment_index"] = missing_segments.index
            parent_lookup = subclones[["chr", "start", "end", *cols_to_fix]].rename(
                columns={
                    "start": "parent_start",
                    "end": "parent_end",
                    **{col: f"parent_{col}" for col in cols_to_fix},
                }
            )
            parent_matches = missing_segments.merge(parent_lookup, on="chr", how="left")
            parent_matches = parent_matches[
                (parent_matches["parent_start"] <= parent_matches["start"])
                & (parent_matches["parent_end"] >= parent_matches["end"])
            ]

            parent_counts = parent_matches["_segment_index"].value_counts()
            multiple_matches = parent_counts[parent_counts > 1]
            if not multiple_matches.empty:
                conflict = segments.loc[multiple_matches.index[0]]
                raise ValueError(
                    f"Multiple matching parent segments found in {subclones_path} "
                    f"for segment {conflict['segment']} "
                    f"(chr{conflict['chr']}:{conflict['start']}-{conflict['end']}) "
                    f"in sample {sample}."
                )
            if not parent_matches.empty:
                parent_matches = parent_matches.drop_duplicates("_segment_index").set_index(
                    "_segment_index"
                )
                for col in cols_to_fix:
                    segments.loc[parent_matches.index, col] = segments.loc[
                        parent_matches.index, col
                    ].fillna(parent_matches[f"parent_{col}"])
        missing_cntot = segments["cntot"].isna().sum()
        if missing_cntot:
            raise ValueError(
                f"Could not map 'cntot' for {missing_cntot} segment(s) in sample {sample} "
                f"using coordinates from {subclones_path}."
            )
        segment_counts = _count_snps_per_segment(
            snps.rename(columns={"seqnames": "chr"}), segments
        )
        segments["heterozygous_SNP_number"] = (
            segments["segment"].map(segment_counts).fillna(0).astype(int)
        )

        purity = float(purity_ploidy["purity"])
        ploidy = float(purity_ploidy["ploidy"])
        if purity <= 0:
            raise ValueError(f"Purity must be > 0 for sample {sample}, got {purity}.")
        segments["cn_tot"] = np.nan  # we don't use _estimate_cn_tot here because we just carry over original Battenber cntot values
        
        # These dummy fileds are introduced only to match the refphase data format
        segments["cn_a"] = np.nan
        segments["cn_b"] = np.nan

        segments["group_name"] = sample
        if args.tumour_id:
            segments["patient_tumour"] = args.tumour_id
        segments["seqnames"] = segments["chr"]
        segments["strand"] = "*"
        segments["width"] = segments["end"] - segments["start"] + 1
        segments["is_LOH"] = False
        segments["mirrored_vs_ref"] = False
        segments["any_ai"] = False
        segments["diptest_pvalue"] = np.nan
        segments["ai_pvalue"] = np.nan
        segments["effect_size"] = np.nan
        segments["is_ai"] = False
        segments["is_reference"] = False
        segments["was_cn_updated"] = False
        segments["homozygous_SNP_number"] = 0
        segments["cn_a_integer"] = np.nan
        segments["cn_b_integer"] = np.nan
        segments['is_clonal'] = segments[f"frac1_{BATTENBERG_SOLUTION_ID}"] == 1
        breakpoint()
        segments = segments[
            [
                "group_name",
                "seqnames",
                "strand",
                "start",
                "end",
                "width",
                "cn_a",
                "cn_b",
                "cn_tot",
                "cntot",
                "is_LOH",
                "mirrored_vs_ref",
                "any_ai",
                "diptest_pvalue",
                "ai_pvalue",
                "effect_size",
                "is_ai",
                "is_reference",
                "was_cn_updated",
                "heterozygous_SNP_number",
                "homozygous_SNP_number",
                "cn_a_integer",
                "cn_b_integer",
                "is_clonal",
            ]
        ].copy()

        phased_segments_all.append(segments)
        phased_snps_all.append(snps)
        purity_ploidy_rows.append(purity_ploidy)

    phased_segments = pd.concat(phased_segments_all, ignore_index=True)
    phased_snps = pd.concat(phased_snps_all, ignore_index=True)
    purity_ploidy_df = pd.DataFrame(purity_ploidy_rows)
    resolved_inventory = pd.DataFrame(resolved_inventory_rows)

    phased_segments.to_csv(output_dir / "phased_segs.tsv", sep="\t", index=False)
    phased_snps.to_csv(output_dir / "phased_snps.tsv", sep="\t", index=False)
    purity_ploidy_df.to_csv(output_dir / "purity_ploidy.tsv", sep="\t", index=False)
    resolved_inventory.to_csv(output_dir / "battenberg_inventory_resolved.tsv", sep="\t", index=False)
    print(f"Wrote Battenberg intermediate tables to {output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

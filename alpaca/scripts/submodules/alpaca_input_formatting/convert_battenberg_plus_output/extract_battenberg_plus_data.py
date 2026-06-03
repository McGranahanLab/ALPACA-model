#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path
import re

import numpy as np
import pandas as pd


def _read_table(path: Path) -> pd.DataFrame:
    # Try delimiter inference first, then common explicit fallbacks.
    for kwargs in (
        {"sep": None, "engine": "python"},
        {"sep": "\t"},
        {"sep": ","},
    ):
        try:
            return pd.read_csv(path, **kwargs)
        except Exception:
            continue
    raise ValueError(f"Could not parse table file: {path}")


def _normalize_chr_value(value) -> int:
    token = re.findall(r"([0-9]+|[XYM][Tt]?)", str(value), flags=re.IGNORECASE)
    if not token:
        raise ValueError(f"Could not parse chromosome value: {value}")
    chrom = token[0].upper()
    chrom = {"X": "23", "Y": "24", "MT": "25", "M": "25"}.get(chrom, chrom)
    return int(chrom)


def _standardize_fractional_columns(df: pd.DataFrame) -> pd.DataFrame:
    col_map = {c.upper(): c for c in df.columns}
    required = {
        "SAMPLE": "sample",
        "CHR": "chr",
        "STARTPOS": "startpos",
        "ENDPOS": "endpos",
        "COPY_NUMBER_A": "copy_number_a",
        "COPY_NUMBER_B": "copy_number_b",
    }
    missing = [src for src in required if src not in col_map]
    if missing:
        raise ValueError(
            "fractional_copy_number file is missing required columns: "
            + ", ".join(missing)
        )

    out = df.rename(columns={col_map[src]: dst for src, dst in required.items()}).copy()
    out["sample"] = out["sample"].astype(str)
    out["chr"] = out["chr"].map(_normalize_chr_value)
    out["startpos"] = pd.to_numeric(out["startpos"], errors="raise").astype(int)
    out["endpos"] = pd.to_numeric(out["endpos"], errors="raise").astype(int)
    out["copy_number_a"] = pd.to_numeric(out["copy_number_a"], errors="raise")
    out["copy_number_b"] = pd.to_numeric(out["copy_number_b"], errors="raise")
    return out


def _standardize_ci_columns(df: pd.DataFrame) -> pd.DataFrame:
    lower = {c.lower(): c for c in df.columns}
    required = {
        "chr": "chr",
        "startpos": "startpos",
        "endpos": "endpos",
        "nmajor": "nmajor",
        "nminor": "nminor",
        "nmajor_ci_lower": "nmajor_ci_lower",
        "nmajor_ci_upper": "nmajor_ci_upper",
        "nminor_ci_lower": "nminor_ci_lower",
        "nminor_ci_upper": "nminor_ci_upper",
    }
    missing = [src for src in required if src not in lower]
    if missing:
        raise ValueError(
            "battenberg_ci file is missing required columns: " + ", ".join(missing)
        )

    out = df.rename(columns={lower[src]: dst for src, dst in required.items()}).copy()
    out["chr"] = out["chr"].map(_normalize_chr_value)
    out["startpos"] = pd.to_numeric(out["startpos"], errors="raise").astype(int)
    out["endpos"] = pd.to_numeric(out["endpos"], errors="raise").astype(int)
    for col in [
        "nmajor",
        "nminor",
        "nmajor_ci_lower",
        "nmajor_ci_upper",
        "nminor_ci_lower",
        "nminor_ci_upper",
    ]:
        out[col] = pd.to_numeric(out[col], errors="raise")

    if "ntot" in lower:
        out["ntot"] = pd.to_numeric(out[lower["ntot"]], errors="coerce")
    else:
        out["ntot"] = out["nmajor"] + out["nminor"]

    return out


def _resolve_ci_file_paths(ci_files: list[str]) -> list[Path]:
    return [Path(p).expanduser().resolve() for p in ci_files if str(p).strip()]


def _infer_sample_from_path(path: Path) -> str:
    # Current battenberg_plus file naming uses '<sample>--<run_id>_..._extended.txt.gz'.
    if "--" in path.name:
        return path.name.split("--", 1)[0]
    stem = path.stem
    for suffix in (".txt", ".tsv", ".csv"):
        if stem.lower().endswith(suffix):
            stem = stem[: -len(suffix)]
    return stem


def _load_ci_tables(
    ci_files: list[str],
    ci_dir: str | None,
    fractional_samples: set[str],
) -> pd.DataFrame:
    ci_paths: list[Path] = []
    if ci_files:
        ci_paths.extend(_resolve_ci_file_paths(ci_files))

    if ci_dir:
        # Include all files with 'extended' substring in the filename.
        directory = Path(ci_dir).expanduser().resolve()
        ci_paths.extend(
            [
                p
                for p in directory.iterdir()
                if p.is_file() and "extended" in p.name.lower()
            ]
        )

    unique_paths = sorted({p.resolve() for p in ci_paths})

    if not unique_paths:
        raise ValueError(
            "No battenberg_plus CI files could be resolved for fractional samples."
        )

    tables = []
    for path in unique_paths:
        ci_df = _standardize_ci_columns(_read_table(path))
        sample = _infer_sample_from_path(path)
        ci_df["sample"] = sample
        tables.append(ci_df)

    out = pd.concat(tables, ignore_index=True)
    # Keep only requested samples; this avoids accidental ingestion of unrelated files in ci_dir.
    len_before = len(out)
    # Since R replaces '-' with '.', unify here:
    out["sample"] = out["sample"].str.replace("-", ".", regex=False)
    out = out[out["sample"].isin(fractional_samples)].copy()
    len_after = len(out)
    if len_before > 0 and len_after == 0:
        raise ValueError(
            "No rows from battenberg_plus CI tables matched samples in fractional copy-number table. "
            f"ci samples: {sorted(out['sample'].unique())}, "
            f"fractional samples: {sorted(fractional_samples)}"
        )
    return out


def _phase_ci_rows(merged: pd.DataFrame, tolerance: float) -> pd.DataFrame:
    direct = np.isclose(
        merged["nmajor"], merged["copy_number_a"], atol=tolerance
    ) & np.isclose(merged["nminor"], merged["copy_number_b"], atol=tolerance)
    flipped = np.isclose(
        merged["nmajor"], merged["copy_number_b"], atol=tolerance
    ) & np.isclose(merged["nminor"], merged["copy_number_a"], atol=tolerance)

    unresolved = ~(direct | flipped)
    if unresolved.any():
        preview = merged.loc[
            unresolved,
            [
                "sample",
                "chr",
                "startpos",
                "endpos",
                "copy_number_a",
                "copy_number_b",
                "nmajor",
                "nminor",
            ],
        ].head(10)
        raise ValueError(
            "Could not determine phasing for some battenberg_plus rows. "
            f"Mismatch count: {int(unresolved.sum())}. Preview:\n{preview.to_string(index=False)}"
        )

    both = direct & flipped
    direct = direct | both
    flipped = flipped & ~both

    out = merged.copy()
    out["cpnA"] = np.where(direct, out["nmajor"], out["nminor"])
    out["cpnB"] = np.where(direct, out["nminor"], out["nmajor"])
    out["lower_CI_A"] = np.where(direct, out["nmajor_ci_lower"], out["nminor_ci_lower"])
    out["upper_CI_A"] = np.where(direct, out["nmajor_ci_upper"], out["nminor_ci_upper"])
    out["lower_CI_B"] = np.where(direct, out["nminor_ci_lower"], out["nmajor_ci_lower"])
    out["upper_CI_B"] = np.where(direct, out["nminor_ci_upper"], out["nmajor_ci_upper"])
    out["was_flipped"] = flipped
    return out


def _write_outputs(output_dir: Path, phased_df: pd.DataFrame, ci_value: float) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)

    seg_out = phased_df[
        [
            "sample",
            "chr",
            "startpos",
            "endpos",
            "cpnA",
            "cpnB",
            "ntot",
        ]
    ].copy()
    seg_out = seg_out.rename(
        columns={
            "sample": "group_name",
            "chr": "seqnames",
            "startpos": "start",
            "endpos": "end",
            "cpnA": "cn_a",
            "cpnB": "cn_b",
            "ntot": "heterozygous_SNP_number",
        }
    )
    seg_out["was_cn_updated"] = False
    seg_out["is_reference"] = False
    seg_out["cntot"] = seg_out["cn_a"] + seg_out["cn_b"]

    ci_out = phased_df[
        [
            "sample",
            "chr",
            "startpos",
            "endpos",
            "cpnA",
            "cpnB",
            "lower_CI_A",
            "upper_CI_A",
            "lower_CI_B",
            "upper_CI_B",
        ]
    ].copy()
    ci_out["segment"] = (
        ci_out["chr"].astype(str)
        + "_"
        + ci_out["startpos"].astype(str)
        + "_"
        + ci_out["endpos"].astype(str)
    )
    ci_out["ci_value"] = ci_value
    ci_out = ci_out[
        [
            "segment",
            "sample",
            "cpnA",
            "cpnB",
            "lower_CI_A",
            "upper_CI_A",
            "lower_CI_B",
            "upper_CI_B",
            "ci_value",
        ]
    ]

    seg_out.to_csv(output_dir / "phased_segs.tsv", sep="\t", index=False)
    ci_out.to_csv(output_dir / "precomputed_ci.tsv", sep="\t", index=False)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Extract battenberg_plus outputs into ALPACA conversion intermediates"
    )
    parser.add_argument(
        "--tumour_id", type=str, required=True, help="Tumour identifier"
    )
    parser.add_argument(
        "--fractional_copy_number",
        type=str,
        required=True,
        help="Path to combined phased fractional copy-number table",
    )
    parser.add_argument(
        "--battenberg_ci_file",
        type=str,
        action="append",
        default=[],
        required=False,
        help="Path to a battenberg_plus CI file. Can be provided multiple times.",
    )
    parser.add_argument(
        "--battenberg_ci_dir",
        type=str,
        required=False,
        help="Optional directory containing one battenberg_plus CI file per sample",
    )
    parser.add_argument(
        "--chromosome",
        type=str,
        required=False,
        help="Optional chromosome filter (e.g. 1, chr1, X)",
    )
    parser.add_argument(
        "--match_tolerance",
        type=float,
        default=0.05,
        help="Absolute tolerance used when matching unphased major/minor copy numbers to phased A/B",
    )
    parser.add_argument(
        "--ci_value",
        type=float,
        default=0.5,
        help="Confidence level label stored in precomputed_ci.tsv",
    )
    parser.add_argument(
        "--output_dir", type=str, required=True, help="Output directory"
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if not args.battenberg_ci_file and not args.battenberg_ci_dir:
        raise ValueError(
            "Provide one or more --battenberg_ci_file and/or --battenberg_ci_dir for battenberg_plus conversion."
        )

    fractional = _standardize_fractional_columns(
        _read_table(Path(args.fractional_copy_number))
    )

    if args.chromosome is not None:
        target_chr = _normalize_chr_value(args.chromosome)
        fractional = fractional[fractional["chr"] == target_chr].copy()

    if fractional.empty:
        raise ValueError(
            "No rows remained in fractional copy-number table after filtering."
        )

    ci_df = _load_ci_tables(
        ci_files=args.battenberg_ci_file,
        ci_dir=args.battenberg_ci_dir,
        fractional_samples=set(fractional["sample"].astype(str).unique()),
    )
    if args.chromosome is not None:
        target_chr = _normalize_chr_value(args.chromosome)
        ci_df = ci_df[ci_df["chr"] == target_chr].copy()
    
    if ci_df.empty:
        raise ValueError(
            "No rows remained in battenberg_plus CI tables after filtering."
        )

    merge_cols = ["sample", "chr", "startpos", "endpos"]
    merged = ci_df.merge(
        fractional,
        on=merge_cols,
        how="inner",
        validate="many_to_one",
    )
    if merged.empty:
        raise ValueError(
            "No overlapping rows between battenberg_plus CI tables and fractional copy-number table. "
            "Check sample names and segment coordinates."
        )

    phased = _phase_ci_rows(merged, tolerance=args.match_tolerance)

    # Keep one row per sample/segment.
    phased = phased.sort_values(
        ["sample", "chr", "startpos", "endpos"]
    ).drop_duplicates(subset=["sample", "chr", "startpos", "endpos"], keep="first")

    _write_outputs(Path(args.output_dir), phased, ci_value=float(args.ci_value))
    print(f"Wrote battenberg_plus intermediates to {args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

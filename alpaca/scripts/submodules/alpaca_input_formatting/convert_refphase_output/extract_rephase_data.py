#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path
import sys
import re

import numpy as np
import pandas as pd

THIS_DIR = Path(__file__).resolve().parent
INPUT_FORMATTING_DIR = THIS_DIR.parent
if str(INPUT_FORMATTING_DIR) not in sys.path:
    sys.path.insert(0, str(INPUT_FORMATTING_DIR))

from r_object_io import get_field, read_rdata, to_dataframe


def _select_refphase_results_object(objects):
    for _, obj in objects.items():
        try:
            get_field(obj, "phased_segs")
            get_field(obj, "phased_snps")
            get_field(obj, "sample_data")
            return obj
        except Exception:
            continue
    return next(iter(objects.values()))


def _unwrap_singleton(value):
    if isinstance(value, np.ndarray):
        if value.shape == ():
            return _unwrap_singleton(value.item())
        if value.size == 1:
            return _unwrap_singleton(value.reshape(-1)[0])
    if isinstance(value, (list, tuple)) and len(value) == 1:
        return _unwrap_singleton(value[0])
    return value


def _flatten_sample_wide_object_df(df: pd.DataFrame, *, name: str) -> pd.DataFrame:
    # Some rdata decodings yield one-row tables where each column is a sample
    # and each cell is a GRanges/GPos-like object. Convert to long format.
    if df.empty or len(df) != 1:
        return df

    reconstructed = []
    for sample_name in df.columns:
        cell = _unwrap_singleton(df.iloc[0][sample_name])
        if isinstance(cell, str):
            return df
        try:
            cell_df = to_dataframe(cell, name=f"{name}${sample_name}")
        except Exception:
            return df
        if cell_df.empty:
            continue
        cell_df = cell_df.copy()
        cell_df["group_name"] = str(sample_name)
        reconstructed.append(cell_df)

    if not reconstructed:
        return df

    out = pd.concat(reconstructed, ignore_index=True, sort=False)
    ordered = ["group_name"] + [c for c in out.columns if c != "group_name"]
    return out.loc[:, ordered]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Extract TSV files from refphase .RData output"
    )
    parser.add_argument(
        "--refphase_rData", type=str, required=True, help="Path to refphase .RData file"
    )
    parser.add_argument("--output_dir", type=str, required=True, help="Output directory")
    parser.add_argument(
        "--chromosome",
        type=str,
        required=False,
        help="Optional chromosome filter (e.g. 1, chr1, X). If set, only this chromosome is kept.",
    )
    return parser.parse_args()


def _normalize_chr_series(series: pd.Series) -> pd.Series:
    if pd.api.types.is_numeric_dtype(series):
        return pd.to_numeric(series, errors="coerce")
    out = series.astype(str).str.extract(r"([0-9]+|[XYM][Tt]?)", flags=re.IGNORECASE)[0]
    out = out.replace({"X": "23", "x": "23", "Y": "24", "y": "24", "MT": "25", "M": "25", "mt": "25", "m": "25"})
    return pd.to_numeric(out, errors="coerce")


def _normalize_chr_value(chr_value: str) -> int:
    extracted = re.findall(r"([0-9]+|[XYM][Tt]?)", str(chr_value), flags=re.IGNORECASE)
    if not extracted:
        raise ValueError(f"Could not parse chromosome value: {chr_value}")
    token = extracted[0].upper()
    token = {"X": "23", "Y": "24", "MT": "25", "M": "25"}.get(token, token)
    return int(token)


def main() -> int:
    args = parse_args()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    print("Loading refphase results from .RData file")
    print(args.refphase_rData)
    objects = read_rdata(args.refphase_rData)
    if not objects:
        raise ValueError("No objects were found in the provided .RData file.")

    refphase_results = _select_refphase_results_object(objects)

    try:
        phased_segments = to_dataframe(
            get_field(refphase_results, "phased_segs"),
            name="refphase_results$phased_segs",
        )
        phased_segments = _flatten_sample_wide_object_df(
            phased_segments, name="refphase_results$phased_segs"
        )
        snps = to_dataframe(
            get_field(refphase_results, "phased_snps"),
            name="refphase_results$phased_snps",
        )
        snps = _flatten_sample_wide_object_df(snps, name="refphase_results$phased_snps")
        purity_ploidy = to_dataframe(
            get_field(refphase_results, "sample_data"),
            name="refphase_results$sample_data",
        )
    except Exception as exc:
        available = sorted(vars(refphase_results).keys()) if hasattr(
            refphase_results, "__dict__"
        ) else "unknown"
        raise RuntimeError(
            "Failed to convert refphase object to tables. "
            f"Top-level fields detected: {available}"
        ) from exc

    if "seqnames" not in snps.columns:
        raise RuntimeError(
            "Converted phased_snps table is missing 'seqnames'. "
            f"Columns detected: {list(snps.columns)}"
        )
    if "group_name" not in snps.columns:
        raise RuntimeError(
            "Converted phased_snps table is missing 'group_name'. "
            f"Columns detected: {list(snps.columns)}"
        )

    if args.chromosome is not None:
        target_chr = _normalize_chr_value(args.chromosome)
        phased_segments = phased_segments[
            _normalize_chr_series(phased_segments["seqnames"]) == target_chr
        ].copy()
        snps = snps[_normalize_chr_series(snps["seqnames"]) == target_chr].copy()
        if phased_segments.empty:
            raise ValueError(
                f"No phased segments remain after --chromosome filter ({args.chromosome})."
            )
        if snps.empty:
            raise ValueError(
                f"No phased SNPs remain after --chromosome filter ({args.chromosome})."
            )

    print("Writing phased segments, snps and purity ploidy to tsv files")
    print(str(output_dir))
    phased_segments.to_csv(output_dir / "phased_segs.tsv", sep="\t", index=False)
    snps.to_csv(output_dir / "phased_snps.tsv", sep="\t", index=False)
    purity_ploidy.to_csv(output_dir / "purity_ploidy.tsv", sep="\t", index=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

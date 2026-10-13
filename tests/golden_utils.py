"""Shared utilities for golden-file comparison in ALPACA tests."""

from pathlib import Path

import pandas as pd

from alpaca.utils import read_dataframe

# Columns whose values vary run-to-run and must be excluded from golden comparisons.
VOLATILE_COLS = ["run_time_seconds"]
# Additional columns that vary when using the Gurobi backend directly.
GUROBI_VOLATILE_COLS = VOLATILE_COLS + ["gurobi_time_CI", "gurobi_time_D"]


def assert_gurobi_log_is_present(log_path):
    """Assert that the Gurobi log file was created. When debugging solution provided by the users
    we are only interested in its feasibility and Gurobi metrics, so the output table is not created"""
    log_path = Path(log_path)
    assert log_path.exists(), (
        f"Gurobi log file not found: {log_path}\n"
        "This likely means that Gurobi did not run at all, so the provided solution was not even checked."
    )


def resolve_output_path(path):
    """Return the existing output file for ``path``, whichever format it was written in.

    ALPACA writes parquet by default and CSV with ``--output_format csv``. Tests refer to
    outputs by their CSV name; this returns the path as given if it exists, otherwise the
    same file with the other extension (.parquet <-> .csv).
    """
    path = Path(path)
    if path.exists():
        return path
    for ext in (".parquet", ".csv"):
        alt = path.with_suffix(ext)
        if alt.exists():
            return alt
    return path


def assert_csv_matches_golden(actual_path, golden_path, ignore_cols=None, update=False):
    """Compare an ALPACA output table (CSV or parquet) against a committed CSV golden file.

    Pass update=True (via the --update-golden flag) to regenerate the golden
    file instead of comparing. Volatile columns listed in ignore_cols are
    dropped before both writing and reading so they never appear in the golden.
    """
    ignore_cols = ignore_cols or []
    actual_path = resolve_output_path(actual_path)
    golden_path = Path(golden_path)

    actual_df = read_dataframe(actual_path)
    actual_df = actual_df.drop(columns=[c for c in ignore_cols if c in actual_df.columns])

    if update:
        golden_path.parent.mkdir(parents=True, exist_ok=True)
        actual_df.to_csv(golden_path, index=False)
        return

    assert golden_path.exists(), (
        f"Golden file not found: {golden_path}\n"
        "Run with --update-golden to generate it."
    )
    golden_df = pd.read_csv(golden_path)
    golden_df = golden_df.drop(columns=[c for c in ignore_cols if c in golden_df.columns])

    pd.testing.assert_frame_equal(
        actual_df.reset_index(drop=True),
        golden_df.reset_index(drop=True),
        check_like=True,
        atol=1e-6,
        rtol=0,
    )

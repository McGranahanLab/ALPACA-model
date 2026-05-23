"""Integration tests for alpaca input-conversion commands.

All tests in this module require Rscript and are automatically skipped when it
is absent from PATH.  They are tagged with ``requires_r`` so they can also be
excluded explicitly::

    pytest tests/ -m "not requires_r"

Golden files (where applicable) live in ``tests/<scenario>/expected/`` and are
updated with ``--update-golden``.
"""

import shutil
import subprocess
import sys

import pandas as pd
import pytest

from golden_utils import VOLATILE_COLS, assert_csv_matches_golden

TUMOUR_ID = "LTX0000-Tumour1"

_RSCRIPT_MISSING = shutil.which("Rscript") is None
_SKIP_R = pytest.mark.skipif(_RSCRIPT_MISSING, reason="Rscript not on PATH")


def _golden(repo_root, scenario):
    return repo_root / "tests" / scenario / "expected"


def _script_dir(repo_root):
    return repo_root / "alpaca" / "scripts" / "submodules" / "alpaca_input_formatting"


# ---------------------------------------------------------------------------
# Full pipeline conversion (uses alpaca input-conversion CLI)
# ---------------------------------------------------------------------------


@pytest.mark.integration
@pytest.mark.requires_r
@_SKIP_R
def test_input_conversion(tmp_path, repo_root, run_alpaca, update_golden):
    """Full CONIPHER + RefPhase → ALPACA format conversion."""
    example_input = repo_root / "examples" / "example_cohort" / "input" / TUMOUR_ID
    refphase_rdata = example_input / f"{TUMOUR_ID}-refphase-results.RData"
    conipher_tree = example_input / f"{TUMOUR_ID}.tree.RDS"

    run_alpaca(
        "input-conversion",
        "--tumour_id", TUMOUR_ID,
        "--refphase_rData", str(refphase_rdata),
        "--CONIPHER_tree_object", str(conipher_tree),
        "--CONIPHER_tree_index", "1",
        "--output_dir", str(tmp_path),
    )

    for expected_file in ("ALPACA_input_table.csv", "ci_table.csv", "cp_table.csv"):
        assert (tmp_path / expected_file).exists(), f"Missing output: {expected_file}"

    golden_dir = _golden(repo_root, "input_conversion")
    for csv_name in ("ALPACA_input_table.csv", "ci_table.csv", "cp_table.csv"):
        assert_csv_matches_golden(
            tmp_path / csv_name,
            golden_dir / csv_name,
            ignore_cols=VOLATILE_COLS,
            update=update_golden,
        )


@pytest.mark.integration
@pytest.mark.requires_r
@_SKIP_R
def test_input_conversion_missing_samples(tmp_path, repo_root):
    """When CONIPHER cp_table lists a sample not present in RefPhase output,
    ci_table samples should be a strict subset of cp_table samples."""
    input_dir = repo_root / "tests" / "missing_samples" / "input" / TUMOUR_ID
    refphase_rdata = input_dir / f"{TUMOUR_ID}-refphase-results.RData"
    script_dir = _script_dir(repo_root) / "convert_refphase_output"

    # Step 1: extract RefPhase RData → TSVs
    extract_script = script_dir / "extract_rephase_data.py"
    result = subprocess.run(
        [sys.executable, str(extract_script),
         "--refphase_rData", str(refphase_rdata),
         "--output_dir", str(tmp_path)],
        capture_output=True, text=True, cwd=repo_root,
    )
    assert result.returncode == 0, result.stderr

    # Step 2: build a cp_table with one extra sample absent from RefPhase,
    # simulating CONIPHER knowing about a sample that RefPhase lacks.
    base_cp = pd.read_csv(
        repo_root / "examples" / "example_cohort" / "input" / TUMOUR_ID / "cp_table.csv",
        index_col="clone",
    )
    base_cp["CONIPHER_EXTRA_SAMPLE"] = 0.0
    cp_table_path = tmp_path / "cp_table.csv"
    base_cp.to_csv(cp_table_path)

    # Step 3: convert TSVs to ALPACA ci_table using the inflated cp_table
    convert_script = script_dir / "convert_refphase.py"
    result = subprocess.run(
        [sys.executable, str(convert_script),
         "--tumour_id", TUMOUR_ID,
         "--output_dir", str(tmp_path),
         "--refphase_segments", str(tmp_path / "phased_segs.tsv"),
         "--refphase_snps", str(tmp_path / "phased_snps.tsv"),
         "--refphase_purity_ploidy", str(tmp_path / "purity_ploidy.tsv"),
         "--conipher_cp_table", str(cp_table_path)],
        capture_output=True, text=True, cwd=repo_root,
    )
    assert result.returncode == 0, result.stderr
    assert (tmp_path / "ci_table.csv").exists()

    cp_samples = set(base_cp.columns)
    ci_samples = set(pd.read_csv(tmp_path / "ci_table.csv")["sample"])
    assert cp_samples != ci_samples, (
        f"Expected sample mismatch between cp_table ({cp_samples}) "
        f"and ci_table ({ci_samples}), but they matched."
    )


# ---------------------------------------------------------------------------
# Individual converter scripts (invoked directly as Python modules)
# ---------------------------------------------------------------------------


@pytest.mark.integration
@pytest.mark.requires_r
@_SKIP_R
def test_convert_conipher(tmp_path, repo_root):
    """Convert a CONIPHER RDS tree object to tree_paths.json."""
    script = (
        _script_dir(repo_root)
        / "convert_conipher_output"
        / "convert_conipher_output.py"
    )
    conipher_tree = (
        repo_root / "examples" / "example_cohort" / "input" / TUMOUR_ID
        / f"{TUMOUR_ID}.tree.RDS"
    )
    result = subprocess.run(
        [sys.executable, str(script),
         "--CONIPHER_tree_object", str(conipher_tree),
         "--output_dir", str(tmp_path)],
        capture_output=True, text=True, cwd=repo_root,
    )
    assert result.returncode == 0, result.stderr
    assert (tmp_path / "tree_paths.json").exists()


@pytest.mark.integration
@pytest.mark.requires_r
@_SKIP_R
def test_extract_and_convert_refphase(tmp_path, repo_root, update_golden):
    """Extract RefPhase RData → TSVs, then convert TSVs → ALPACA CI table."""
    example_input = repo_root / "examples" / "example_cohort" / "input" / TUMOUR_ID
    refphase_rdata = example_input / f"{TUMOUR_ID}-refphase-results.RData"
    script_dir = _script_dir(repo_root) / "convert_refphase_output"

    # Step 1: extract RData to TSVs
    extract_script = script_dir / "extract_rephase_data.py"
    result = subprocess.run(
        [sys.executable, str(extract_script),
         "--refphase_rData", str(refphase_rdata),
         "--output_dir", str(tmp_path)],
        capture_output=True, text=True, cwd=repo_root,
    )
    assert result.returncode == 0, result.stderr
    for tsv in ("phased_segs.tsv", "phased_snps.tsv", "purity_ploidy.tsv"):
        assert (tmp_path / tsv).exists(), f"Extraction did not produce {tsv}"

    # Step 2: convert TSVs to ALPACA ci_table
    conipher_cp_table = example_input / "cp_table.csv"
    convert_script = script_dir / "convert_refphase.py"
    result = subprocess.run(
        [sys.executable, str(convert_script),
         "--tumour_id", TUMOUR_ID,
         "--output_dir", str(tmp_path),
         "--refphase_segments", str(tmp_path / "phased_segs.tsv"),
         "--refphase_snps", str(tmp_path / "phased_snps.tsv"),
         "--refphase_purity_ploidy", str(tmp_path / "purity_ploidy.tsv"),
         "--conipher_cp_table", str(conipher_cp_table)],
        capture_output=True, text=True, cwd=repo_root,
    )
    assert result.returncode == 0, result.stderr
    assert (tmp_path / "ci_table.csv").exists()

    golden_dir = _golden(repo_root, "refphase_conversion")
    assert_csv_matches_golden(
        tmp_path / "ci_table.csv",
        golden_dir / "ci_table.csv",
        ignore_cols=VOLATILE_COLS,
        update=update_golden,
    )


@pytest.mark.integration
@pytest.mark.requires_r
@_SKIP_R
def test_convert_ascat(tmp_path, repo_root, update_golden):
    """Convert ASCAT RData to ALPACA format.

    Skipped if the test dataset (in dev/, which is gitignored) is absent.
    """
    ascat_rdata = repo_root / "dev" / "ASCAT_test_data" / "patient_1" / "ASCAT_objects.Rdata"
    if not ascat_rdata.exists():
        pytest.skip(f"ASCAT test dataset not found: {ascat_rdata}")

    script_dir = _script_dir(repo_root) / "convert_ascat_output"

    # Step 1: unpack RData → TSVs using R helper
    r_helper = script_dir / "unpack_RDS_ascat_helper.R"
    result = subprocess.run(
        ["Rscript", str(r_helper),
         "--rdata", str(ascat_rdata),
         "--output_dir", str(tmp_path),
         "--out_prefix", TUMOUR_ID],
        capture_output=True, text=True, cwd=repo_root,
    )
    assert result.returncode == 0, result.stderr

    # Step 2: convert TSVs → ALPACA format
    convert_script = script_dir / "convert_ascat_output.py"
    result = subprocess.run(
        [sys.executable, str(convert_script),
         "--tumour_id", TUMOUR_ID,
         "--output_dir", str(tmp_path),
         "--segments", str(tmp_path / f"{TUMOUR_ID}_segments.tsv"),
         "--snps", str(tmp_path / f"{TUMOUR_ID}_snps.tsv"),
         "--purity_ploidy", str(tmp_path / f"{TUMOUR_ID}_purity_ploidy.tsv"),
         "--n_boot", "10",
         "--gamma", "1"],
        capture_output=True, text=True, cwd=repo_root,
    )
    assert result.returncode == 0, result.stderr
    assert (tmp_path / "ALPACA_input_table.csv").exists()

    golden_dir = _golden(repo_root, "ascat_conversion")
    assert_csv_matches_golden(
        tmp_path / "ALPACA_input_table.csv",
        golden_dir / "ALPACA_input_table.csv",
        ignore_cols=VOLATILE_COLS,
        update=update_golden,
    )

"""End-to-end runs of single clone / single sample tumours in tumour and segment mode."""

import pandas as pd
import pytest

from alpaca.utils import read_dataframe
from mock_inputs import (
    SINGLE_CLONE_PROPORTIONS,
    SINGLE_CLONE_TREE_PATHS,
    SINGLE_CLONE_TRUTH,
    TUMOUR_ID,
    make_mock_tumour_dir,
)

pytestmark = pytest.mark.integration

_PYOMO_SCIP = ("--solver", "pyomo", "--pyomo_solver", "scip", "--pyomo_solver_options", "seed=1 threads=1")
SEG1, SEG2 = list(SINGLE_CLONE_TRUTH)
SAMPLES3 = ("S1", "S2", "S3")


def _override_observations(tumour_dir, segment, observed, half_width=0.2):
    """Replace one segment's observed copy numbers ({sample: (cpnA, cpnB)}) and CIs."""
    paths = [tumour_dir / "ALPACA_input_table.csv", *(tumour_dir / "segments").glob(f"*_{segment}.csv")]
    for path in paths:
        df = pd.read_csv(path)
        for sample, (a, b) in observed.items():
            df.loc[(df["sample"] == sample) & (df["segment"] == segment), ["cpnA", "cpnB"]] = [a, b]
        df.to_csv(path, index=False)
    ci = pd.read_csv(tumour_dir / "ci_table.csv")
    for sample, (a, b) in observed.items():
        rows = (ci["sample"] == sample) & (ci["segment"] == segment)
        ci.loc[rows, ["lower_CI_A", "upper_CI_A", "lower_CI_B", "upper_CI_B"]] = [
            max(a - half_width, 0), a + half_width, max(b - half_width, 0), b + half_width
        ]
    ci.to_csv(tumour_dir / "ci_table.csv", index=False)


def _single_clone_dir(root, samples):
    return make_mock_tumour_dir(
        root,
        samples=samples,
        tree_paths=SINGLE_CLONE_TREE_PATHS,
        proportions=SINGLE_CLONE_PROPORTIONS,
        truth=SINGLE_CLONE_TRUTH,
    )


@pytest.fixture
def single_clone_dir(tmp_path):
    tumour_dir = _single_clone_dir(tmp_path / "in", SAMPLES3)
    # the samples disagree on allele A: no integer lies in every sample's CI
    _override_observations(tumour_dir, SEG2, {"S1": (2.0, 0.0), "S2": (3.0, 0.0), "S3": (2.5, 0.0)})
    return tumour_dir


@pytest.fixture
def single_clone_single_sample_dir(tmp_path):
    tumour_dir = _single_clone_dir(tmp_path / "in", ("MOCK_S1",))
    _override_observations(tumour_dir, SEG2, {"MOCK_S1": (2.5, 0.0)})
    return tumour_dir


@pytest.fixture
def single_sample_dir(tmp_path):
    # default mock: 4 clones, one sample; segment 1 is subclonal (A=2.4, B=0.7)
    return make_mock_tumour_dir(tmp_path / "in")


def _run_tumour(run_alpaca, tumour_dir, out_dir, *extra):
    run_alpaca(
        "run", "--input_tumour_directory", str(tumour_dir), "--output_directory", str(out_dir),
        "--debug", *_PYOMO_SCIP, *extra,
    )
    return read_dataframe(out_dir / f"ALPACA_output_{TUMOUR_ID}.parquet")


def _edge_case_report(out_dir):
    return pd.read_csv(out_dir / "logs_and_reports" / "edge_case_report.csv").set_index("segment")


def test_single_clone_tumour_mode(tmp_path, run_alpaca, single_clone_dir):
    out = tmp_path / "out"
    df = _run_tumour(run_alpaca, single_clone_dir, out).set_index("segment")
    assert (df.loc[SEG1, "pred_CN_A"], df.loc[SEG1, "pred_CN_B"]) == SINGLE_CLONE_TRUTH[SEG1]["clone1"]
    assert (df.loc[SEG2, "pred_CN_A"], df.loc[SEG2, "pred_CN_B"]) == (2, 0)

    warning = (out / "warning_single_clone.txt").read_text()
    assert "clones: 1, samples: 3" in warning
    assert "Segments flagged: 1" in warning
    report = _edge_case_report(out)
    assert set(report["solver_path"]) == {"single_clone_direct"}
    assert pd.isna(report.loc[SEG1, "flags"])
    assert "no_common_integer_A" in report.loc[SEG2, "flags"]
    # per-segment reports are folded into the combined one
    assert not list(out.rglob(f"{TUMOUR_ID}_*_edge_case_report.*"))


def test_single_clone_single_sample_tumour_mode(tmp_path, run_alpaca, single_clone_single_sample_dir):
    out = tmp_path / "out"
    df = _run_tumour(run_alpaca, single_clone_single_sample_dir, out).set_index("segment")
    assert (df.loc[SEG1, "pred_CN_A"], df.loc[SEG1, "pred_CN_B"]) == SINGLE_CLONE_TRUTH[SEG1]["clone1"]
    warning = (out / "warning_single_clone_and_single_sample.txt").read_text()
    assert "clones: 1, samples: 1" in warning
    flags = _edge_case_report(out).loc[SEG2, "flags"]
    assert "ci_contains_no_integer_A" in flags and "rounding_tie_A" in flags


def test_single_sample_tumour_mode(tmp_path, run_alpaca, single_sample_dir):
    out = tmp_path / "out"
    _run_tumour(run_alpaca, single_sample_dir, out)
    warning = (out / "warning_single_sample.txt").read_text()
    assert "clones: 4, samples: 1" in warning
    report = _edge_case_report(out)
    assert set(report["solver_path"]) == {"ilp"}
    assert report["restrict_heterogeneity"].all()
    assert report.loc["1_1000000_5000000", "flags"] == "ambiguous"
    assert "equally good alternative" in report.loc["1_1000000_5000000", "details"]


def test_single_sample_restriction_can_be_disabled(tmp_path, run_alpaca, single_sample_dir):
    out = tmp_path / "out"
    _run_tumour(run_alpaca, single_sample_dir, out, "--single_sample_restrict_heterogeneity", "0")
    assert not _edge_case_report(out)["restrict_heterogeneity"].any()


def test_standard_tumour_has_no_warning(tmp_path, run_alpaca):
    tumour_dir = make_mock_tumour_dir(tmp_path / "in", segments=[SEG1], samples=("S1", "S2"))
    out = tmp_path / "out"
    _run_tumour(run_alpaca, tumour_dir, out)
    assert not list(out.glob("warning_*.txt"))
    assert not (out / "logs_and_reports" / "edge_case_report.csv").exists()


@pytest.mark.parametrize(
    "fixture,category",
    [
        ("single_clone_dir", "single_clone"),
        ("single_sample_dir", "single_sample"),
        ("single_clone_single_sample_dir", "single_clone_and_single_sample"),
    ],
)
def test_segment_mode_writes_per_segment_reports(tmp_path, run_alpaca, request, fixture, category):
    tumour_dir = request.getfixturevalue(fixture)
    segment_files = sorted(p.name for p in (tumour_dir / "segments").glob("*.csv"))
    out = tmp_path / "out"
    run_alpaca(
        "run", "--mode", "segment",
        "--input_tumour_directory", str(tumour_dir),
        "--input_data_directory", str(tumour_dir / "segments"),
        "--input_files", *segment_files,
        "--output_directory", str(out),
        "--debug", *_PYOMO_SCIP,
    )
    for name in segment_files:
        segment = name.replace(f"ALPACA_input_table_{TUMOUR_ID}_", "").replace(".csv", "")
        assert (out / f"optimal_{TUMOUR_ID}_{segment}.parquet").exists()
        report = read_dataframe(out / f"{TUMOUR_ID}_{segment}_edge_case_report.parquet")
        assert report["category"].tolist() == [category]
    # tumour-level warning files are only written in tumour mode (Nextflow writes them after merging)
    assert not list(out.glob("warning_*.txt"))


@pytest.mark.parametrize("solver", [_PYOMO_SCIP, pytest.param(("--solver", "gurobi"), marks=pytest.mark.requires_gurobi)])
def test_single_sample_clonal_segment_is_feasible(tmp_path, run_alpaca, solver):
    """Restricting heterogeneity fixes every clone to (2,1), which complexity 0 cannot reach.

    The iterations must start at the lowest feasible complexity rather than dropping the
    segment as infeasible.
    """
    if "gurobi" in solver:
        pytest.importorskip("gurobipy")
    clones = ("clone1", "clone2", "clone3", "clone4")
    tumour_dir = make_mock_tumour_dir(tmp_path / "in", segments=[SEG1], truth={SEG1: {c: (2, 1) for c in clones}})
    out = tmp_path / "out"
    run_alpaca(
        "run", "--input_tumour_directory", str(tumour_dir), "--output_directory", str(out), "--debug", *solver,
    )
    df = read_dataframe(out / f"ALPACA_output_{TUMOUR_ID}.parquet")
    assert set(zip(df["pred_CN_A"], df["pred_CN_B"])) == {(2, 1)}
    infeasible = pd.read_csv(out / "logs_and_reports" / "infeasibility_report.csv")
    assert infeasible.empty
    assert _edge_case_report(out).loc[SEG1, "restrict_heterogeneity"]

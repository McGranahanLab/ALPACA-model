"""Unit tests for single-clone / single-sample detection, solving and reporting."""

import pandas as pd
import pytest

from alpaca.edge_cases import (
    SINGLE_CLONE,
    SINGLE_CLONE_AND_SINGLE_SAMPLE,
    SINGLE_SAMPLE,
    STANDARD,
    classify_tumour,
    count_integers,
    make_segment_report,
    solve_single_clone,
    write_tumour_edge_case_outputs,
)

SMALL_SEG = "1_1000000_5000000"
LARGE_SEG = "1_1000000_90000000"
MODEL_CONFIG = {
    "limit_homozygous_deletions_threshold_flag": True,
    "homozygous_deletion_threshold": 1,
    "homo_del_strong_evidence_threshold": 0.5,
    "homo_del_size_limit": 5 * 10**7,
}


def _cp(clones, samples):
    return pd.DataFrame(1.0 / len(clones), index=pd.Index(clones, name="clone"), columns=samples)


@pytest.mark.parametrize(
    "clones,samples,expected",
    [
        (["c1", "c2"], ["S1", "S2"], STANDARD),
        (["c1"], ["S1", "S2"], SINGLE_CLONE),
        (["c1", "c2"], ["S1"], SINGLE_SAMPLE),
        (["c1"], ["S1"], SINGLE_CLONE_AND_SINGLE_SAMPLE),
    ],
)
def test_classify_tumour(clones, samples, expected):
    assert classify_tumour(_cp(clones, samples)) == expected


@pytest.mark.parametrize(
    "lower,upper,expected",
    [(0.0, 0.35, 1), (2.1, 2.4, 0), (1.8, 2.2, 1), (1.0, 2.0, 2), (0.6, 2.4, 2)],
)
def test_count_integers(lower, upper, expected):
    assert count_integers(lower, upper) == expected


def _tables(obs, segment=SMALL_SEG, half_width=0.2):
    """obs: {sample: (cpnA, cpnB)} -> (input_table, ci_table)."""
    it = pd.DataFrame(
        [{"tumour_id": "T", "sample": s, "segment": segment, "cpnA": a, "cpnB": b} for s, (a, b) in obs.items()]
    )
    ci = pd.DataFrame(
        [
            {
                "sample": s,
                "segment": segment,
                "lower_CI_A": max(a - half_width, 0),
                "upper_CI_A": a + half_width,
                "lower_CI_B": max(b - half_width, 0),
                "upper_CI_B": b + half_width,
            }
            for s, (a, b) in obs.items()
        ]
    )
    return it, ci


def test_single_clone_consistent_samples():
    it, ci = _tables({"S1": (2.1, 1.0), "S2": (1.85, 0.95), "S3": (2.2, 1.1)})
    (a, b), scores, flags, _ = solve_single_clone(it, ci, SMALL_SEG, MODEL_CONFIG)
    assert (a, b) == (2, 1)
    assert flags == []
    assert scores["CI_score"] == 0
    assert scores["complexity"] == 1  # one gain from diploid


def test_single_clone_contradictory_samples_flags_no_common_integer_and_tie():
    it, ci = _tables({"S1": (2.0, 1.0), "S2": (3.0, 1.0), "S3": (2.5, 1.0)})
    (a, b), scores, flags, details = solve_single_clone(it, ci, SMALL_SEG, MODEL_CONFIG)
    # A=2 and A=3 both miss two CIs at summed distance 1.5; the state closer to diploid wins
    assert (a, b) == (2, 1)
    assert "no_common_integer_A" in flags
    assert "rounding_tie_A" in flags
    assert not any(f.endswith("_B") for f in flags)
    assert scores["CI_score"] == 2
    assert "chose 2 (tied states [2, 3])" in details


def test_single_clone_single_sample_non_integer():
    it, ci = _tables({"S1": (2.5, 1.0)})
    (a, b), _, flags, _ = solve_single_clone(it, ci, SMALL_SEG, MODEL_CONFIG)
    assert (a, b) == (2, 1)
    assert set(flags) == {"ci_contains_no_integer_A", "rounding_tie_A"}


def test_single_clone_homozygous_deletion_allowed_with_strong_evidence():
    it, ci = _tables({"S1": (0.2, 0.1)}, segment=LARGE_SEG)
    (a, b), _, _, _ = solve_single_clone(it, ci, LARGE_SEG, MODEL_CONFIG)
    assert (a, b) == (0, 0)


def test_single_clone_homozygous_deletion_blocked_on_large_segment():
    config = {**MODEL_CONFIG, "homo_del_strong_evidence_threshold": 0.1}
    it, ci = _tables({"S1": (0.3, 0.3)}, segment=LARGE_SEG)
    (a, b), _, _, _ = solve_single_clone(it, ci, LARGE_SEG, config)
    assert a + b >= 1


def test_write_tumour_edge_case_outputs(tmp_path):
    tumour_dir = tmp_path / "T1"
    tumour_dir.mkdir()
    _cp(["c1"], ["S1", "S2"]).to_csv(tumour_dir / "cp_table.csv")
    reports = tmp_path / "reports"
    reports.mkdir()
    for seg, flags in (("1_10_100", []), ("2_10_100", ["no_common_integer_A"])):
        make_segment_report("T1", seg, SINGLE_CLONE, 1, 2, "single_clone_direct", False, flags, "A: ...").to_parquet(
            reports / f"T1_{seg}_edge_case_report.parquet"
        )
    # a report from another tumour sharing the directory must be ignored
    make_segment_report("T10", "1_10_100", SINGLE_CLONE, 1, 2, "single_clone_direct", False, []).to_csv(
        reports / "T10_1_10_100_edge_case_report.csv", index=False
    )
    out = tmp_path / "out"
    category = write_tumour_edge_case_outputs(str(tumour_dir), str(reports), str(out), cohort_status="excluded")
    assert category == SINGLE_CLONE
    text = (out / "warning_single_clone.txt").read_text()
    assert "Segments solved: 2" in text
    assert "Segments flagged: 1" in text
    assert "2_10_100\tno_common_integer_A" in text
    assert "Cohort status: excluded" in text
    combined = pd.read_csv(out / "edge_case_report.csv")
    assert combined["tumour_id"].unique().tolist() == ["T1"]


def test_write_tumour_edge_case_outputs_standard_tumour_writes_nothing(tmp_path):
    tumour_dir = tmp_path / "T1"
    tumour_dir.mkdir()
    _cp(["c1", "c2"], ["S1", "S2"]).to_csv(tumour_dir / "cp_table.csv")
    out = tmp_path / "out"
    assert write_tumour_edge_case_outputs(str(tumour_dir), str(tmp_path), str(out)) is None
    assert not out.exists()


# ---------------------------------------------------------------------------
# heterogeneity restriction and ambiguity search in the ILP backends
# ---------------------------------------------------------------------------

from alpaca.solvers import ModelInputs, create_solver_backend  # noqa: E402
from alpaca.utils import homogeneous_copy_number  # noqa: E402


@pytest.mark.parametrize(
    "lowers,uppers,expected",
    [
        ([1.8], [2.2], 2),
        ([0.0], [0.35], 0),  # truncating the bounds used to report "no integer" here
        ([2.1], [2.4], None),  # no integer: heterogeneity allowed
        ([1.8, 2.1], [2.2, 2.4], None),  # any sample without an integer allows heterogeneity
        ([0.6], [2.4], None),  # two integers: not restricted
        ([], [], None),
    ],
)
def test_homogeneous_copy_number(lowers, uppers, expected):
    assert homogeneous_copy_number(lowers, uppers) == expected


def _single_sample_inputs(cpn_a, cpn_b, half_width=0.2):
    it, ci = _tables({"S1": (cpn_a, cpn_b)}, half_width=half_width)
    proportions = {"clone1": 0.4, "clone2": 0.2, "clone3": 0.1, "clone4": 0.3}
    return ModelInputs(
        segment=SMALL_SEG,
        ci_table=ci,
        fractional_copy_number_table=it,
        tree=[["clone1", "clone2", "clone3"], ["clone1", "clone4"]],
        clone_proportions=pd.DataFrame({"S1": proportions}).rename_axis("clone"),
    )


BACKENDS = [
    pytest.param("pyomo", {"pyomo_solver": "scip"}, id="pyomo"),
    pytest.param("gurobi", {"time_limit": 30}, id="gurobi", marks=pytest.mark.requires_gurobi),
]


def _backend(name, config, inputs):
    if name == "gurobi":
        pytest.importorskip("gurobipy")
    return create_solver_backend(name, inputs, config)


@pytest.mark.parametrize("name,config", BACKENDS)
def test_restrict_heterogeneity_makes_single_sample_solution_clonal(name, config):
    # A=2.1 with CI [1.9, 2.3] holds only the integer 2; B=0.1 with CI [0, 0.3] holds only 0
    inputs = _single_sample_inputs(2.1, 0.1)
    solution = _backend(name, {**config, "restrict_heterogeneity_flag": True}, inputs).solve(10).solution
    solution = solution[solution.clone != "diploid"]
    assert set(solution["pred_CN_A"]) == {2}
    assert set(solution["pred_CN_B"]) == {0}


@pytest.mark.parametrize("name,config", BACKENDS)
def test_find_alternative(name, config):
    # 0.4*2 + 0.2*3 + 0.1*4 + 0.3*2 = 2.4: several assignments explain the observation
    inputs = _single_sample_inputs(2.4, 0.7)
    result = _backend(name, config, inputs).solve(4).solution
    result = result[result.clone != "diploid"].set_index("clone")
    reference = {(a, c): int(result.loc[c, f"pred_CN_{a}"]) for c in result.index for a in ("A", "B")}
    max_ci = float(result["CI_score"].iloc[0])
    max_d = float(result["D_score"].iloc[0]) + 1e-3
    alternative = _backend(name, config, inputs).find_alternative(4, reference, max_ci, max_d)
    assert alternative is not None
    alternative = alternative[alternative.clone != "diploid"].set_index("clone")
    assert any(
        int(alternative.loc[c, f"pred_CN_{a}"]) != v for (a, c), v in reference.items()
    )
    # a diploid segment at complexity 0 has a unique solution
    inputs = _single_sample_inputs(1.0, 1.0)
    reference = {(a, c): 1 for c in ("clone1", "clone2", "clone3", "clone4") for a in ("A", "B")}
    assert _backend(name, config, inputs).find_alternative(0, reference, 0, 1e-3) is None


def _elbow_segment(category):
    from alpaca.ALPACA_segment_solution_class import SegmentSolution

    seg = SegmentSolution.__new__(SegmentSolution)
    seg.tumour_category = category
    seg.tumour_id, seg.segment = "T1", SMALL_SEG
    seg.maximum_complexity = 20
    seg.s_type = "s_strictly_decreasing"
    seg.config = {"model_config": {"d_zero": False}}
    seg.metrics = {
        "solutions": [
            pd.DataFrame({"complexity": [c], "D_score": [5.0 - c], "CI_score": [0.0], "allowed_complexity": [c]})
            for c in range(5)
        ]
    }
    return seg


def test_single_sample_elbow_failure_selects_most_complex(monkeypatch):
    import alpaca.ALPACA_segment_solution_class as module

    monkeypatch.setattr(module, "find_s_values", lambda *args, **kwargs: (None, "default"))
    seg = _elbow_segment(SINGLE_SAMPLE)
    seg.find_elbow()
    assert seg.optimal_solution_index == 4
    assert seg.elbow["raw_code"] == "elbow_failure"
    assert seg.elbow["dec_code"] == "elbow_failure"

    with pytest.raises(AssertionError, match="S_raw not found"):
        _elbow_segment(STANDARD).find_elbow()

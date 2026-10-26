import os

import pandas as pd
import pytest

from alpaca.solvers import ModelInputs, create_solver_backend
from alpaca.solvers.base import SolverFactoryError
from alpaca.solvers.pyomo_backend import PyomoBackend


def _make_minimal_inputs() -> ModelInputs:
    segment = "1_0_1"
    ci_table = pd.DataFrame(
        {
            "segment": [segment],
            "sample": ["S1"],
            "lower_CI_A": [0.5],
            "upper_CI_A": [1.5],
            "lower_CI_B": [0.5],
            "upper_CI_B": [1.5],
            "tumour_id": ["T1"],
            "ci_value": [0.5],
        }
    )
    fractional = pd.DataFrame(
        {
            "segment": [segment],
            "sample": ["S1"],
            "tumour_id": ["T1"],
            "cpnA": [1.0],
            "cpnB": [1.0],
        }
    )
    tree = [["clone1", "clone2"]]
    cp_table = pd.DataFrame([[1.0]], columns=["S1"], index=["clone2"])
    cp_table.loc["clone1"] = 0.0
    cp_table = cp_table.sort_index()
    return ModelInputs(
        segment=segment,
        ci_table=ci_table,
        fractional_copy_number_table=fractional,
        tree=tree,
        clone_proportions=cp_table,
    )


def test_create_solver_backend_pyomo_returns_backend():
    inputs = _make_minimal_inputs()
    backend = create_solver_backend(
        "pyomo",
        inputs,
        {
            "solver": "pyomo",
            "pyomo_solver": "cbc",
            "minimise_events_to_diploid": False,
        },
    )
    assert isinstance(backend, PyomoBackend)
    # verify backend recorded solver choice without running optimization
    assert backend.solver_name == "cbc"


def test_create_solver_backend_unknown():
    inputs = _make_minimal_inputs()
    with pytest.raises(SolverFactoryError):
        create_solver_backend("unknown", inputs, {})


def test_pyomo_backend_raises_scip_probe_timeout(monkeypatch):
    inputs = _make_minimal_inputs()
    monkeypatch.delenv("PYOMO_SOLVER_EXEC_TIMEOUT", raising=False)

    backend = PyomoBackend(
        inputs,
        {
            "solver": "pyomo",
            "pyomo_solver": "scip",
            "minimise_events_to_diploid": False,
        },
    )

    assert backend.solver_name == "scip"
    assert os.environ["PYOMO_SOLVER_EXEC_TIMEOUT"] == "30"


def test_pyomo_backend_keeps_existing_longer_scip_timeout(monkeypatch):
    inputs = _make_minimal_inputs()
    monkeypatch.setenv("PYOMO_SOLVER_EXEC_TIMEOUT", "45")

    PyomoBackend(
        inputs,
        {
            "solver": "pyomo",
            "pyomo_solver": "scipampl",
            "minimise_events_to_diploid": False,
        },
    )

    assert os.environ["PYOMO_SOLVER_EXEC_TIMEOUT"] == "45"


def _single_clone_inputs(cpn_a=2.0, cpn_b=1.0) -> ModelInputs:
    segment = "1_1000000_5000000"
    return ModelInputs(
        segment=segment,
        ci_table=pd.DataFrame(
            {
                "segment": [segment],
                "sample": ["S1"],
                "lower_CI_A": [cpn_a - 0.2],
                "upper_CI_A": [cpn_a + 0.2],
                "lower_CI_B": [cpn_b - 0.2],
                "upper_CI_B": [cpn_b + 0.2],
            }
        ),
        fractional_copy_number_table=pd.DataFrame(
            {"segment": [segment], "sample": ["S1"], "tumour_id": ["T1"], "cpnA": [cpn_a], "cpnB": [cpn_b]}
        ),
        tree=[["clone1"]],
        clone_proportions=pd.DataFrame({"S1": [1.0]}, index=pd.Index(["clone1"], name="clone")),
    )


def test_pyomo_diploid_pseudo_clone_is_fixed_to_one():
    """The trunk's change from diploid must count towards complexity (as in the Gurobi model)."""
    from pyomo import environ as pyo

    backend = create_solver_backend("pyomo", _single_clone_inputs(), {"pyomo_solver": "scip"})
    result = backend.solve(allowed_complexity=0)
    m = backend.model
    assert pyo.value(m.X["A", "diploid"]) == 1
    assert pyo.value(m.X["B", "diploid"]) == 1
    # with no events allowed the clone must stay diploid
    assert result.solution[["pred_CN_A", "pred_CN_B"]].iloc[0].tolist() == [1, 1]

    result = create_solver_backend("pyomo", _single_clone_inputs(), {"pyomo_solver": "scip"}).solve(
        allowed_complexity=1
    )
    row = result.solution.iloc[0]
    assert (row["pred_CN_A"], row["pred_CN_B"], row["complexity"]) == (2, 1, 1)


@pytest.mark.requires_gurobi
@pytest.mark.parametrize("allowed_complexity", [0, 1, 2, 4])
def test_pyomo_and_gurobi_agree_on_scores(allowed_complexity):
    pytest.importorskip("gurobipy")
    from mock_inputs import CLONE_PROPORTIONS, SEGMENT_TRUTH, TREE_PATHS, observed_copy_numbers

    segment = "1_1000000_5000000"
    cpn_a, cpn_b = observed_copy_numbers(SEGMENT_TRUTH[segment], CLONE_PROPORTIONS)
    inputs = _single_clone_inputs(cpn_a, cpn_b)
    inputs.tree = TREE_PATHS
    inputs.clone_proportions = pd.DataFrame(
        {"S1": CLONE_PROPORTIONS}
    ).rename_axis("clone")
    scores = {}
    for name, config in (("pyomo", {"pyomo_solver": "scip"}), ("gurobi", {"time_limit": 30})):
        solution = create_solver_backend(name, inputs, config).solve(allowed_complexity).solution
        solution = solution[solution.clone != "diploid"]
        scores[name] = solution[["complexity", "CI_score", "D_score"]].iloc[0].tolist()
    assert scores["pyomo"] == scores["gurobi"]

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

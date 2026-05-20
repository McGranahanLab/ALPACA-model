"""Factory helpers for solver backends."""

from __future__ import annotations

from typing import Dict, Type

from .base import ModelInputs, SolverBackend, SolverFactoryError
from .gurobi_backend import GurobiBackend
from .pyomo_backend import PyomoBackend

AVAILABLE_SOLVERS: Dict[str, Type[SolverBackend]] = {
    "gurobi": GurobiBackend,
    "pyomo": PyomoBackend,
}


def create_solver_backend(name: str, inputs: ModelInputs, config: dict) -> SolverBackend:
    key = (name or "pyomo").lower()
    if key not in AVAILABLE_SOLVERS:
        raise SolverFactoryError(f"Unknown solver backend '{name}'. Available: {sorted(AVAILABLE_SOLVERS.keys())}")
    backend_cls = AVAILABLE_SOLVERS[key]
    return backend_cls(inputs, config)

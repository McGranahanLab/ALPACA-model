"""Common dataclasses and interfaces for ALPACA solver backends."""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple

import pandas as pd


@dataclass
class ModelInputs:
    """Container with immutable inputs required to instantiate a solver backend."""

    segment: str
    ci_table: pd.DataFrame
    fractional_copy_number_table: pd.DataFrame
    tree: List[List[str]]
    clone_proportions: pd.DataFrame


@dataclass
class SolverResult:
    """Normalized output returned by each solver backend."""

    solution: pd.DataFrame
    runtime: Optional[float]
    backend_name: str
    gap_status: Dict[str, Any] = field(default_factory=dict)
    raw_model: Any = None
    wrapper: Any = None
    objective_values: Dict[str, float] = field(default_factory=dict)
    solver_specific: Dict[str, Any] = field(default_factory=dict)


class SolverBackend(ABC):
    """Abstract base class implemented by all solver backends."""

    name: str

    def __init__(self, inputs: ModelInputs, config: Dict[str, Any]):
        self.inputs = inputs
        self.config = config

    @abstractmethod
    def solve(self, allowed_complexity: int) -> SolverResult:
        """Build, solve, and return a normalized ALPACA solution."""

    @abstractmethod
    def find_alternative(
        self,
        allowed_complexity: int,
        reference: Dict[Tuple[str, str], int],
        max_ci: Optional[float],
        max_d: Optional[float],
    ) -> Optional[pd.DataFrame]:
        """Look for a solution as good as a reference one but with different copy numbers.

        ``reference`` maps (allele, clone) to the reference integer copy number; at
        least one of these must differ. The CI and D objectives are capped at
        ``max_ci`` / ``max_d`` (None leaves an objective uncapped). Returns the
        alternative solution, or None if none exists.
        """

    def supports_feature(self, feature_name: str) -> bool:
        """Override to expose optional features (e.g., IIS reports)."""

        return False


class SolverFactoryError(RuntimeError):
    """Raised when an unknown backend name is requested."""

    pass

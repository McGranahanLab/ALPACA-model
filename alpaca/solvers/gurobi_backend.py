"""Gurobi-backed implementation of the ALPACA solver interface."""

from __future__ import annotations

from typing import Dict, Optional

# get Gurobi backend from the original implementation
from alpaca.ALPACA_model_class import Model as GurobiModel

from .base import ModelInputs, SolverBackend, SolverResult


class GurobiBackend(SolverBackend):
    name = "gurobi"

    def solve(self, allowed_complexity: int) -> SolverResult:
        model_kwargs = {**self.config, "allowed_tree_complexity": allowed_complexity}
        model = GurobiModel(
            segment=self.inputs.segment,
            ci_table=self.inputs.ci_table,
            fractional_copy_number_table=self.inputs.fractional_copy_number_table,
            tree=self.inputs.tree,
            clone_proportions=self.inputs.clone_proportions,
            **model_kwargs,
        )
        model.model.optimize()
        model.get_output()
        runtime = getattr(model.model, "Runtime", None)
        objective_values: Dict[str, float] = {}
        if getattr(model, "add_CI_objective", False):
            objective_values["CI"] = float(model.Z.getValue())
        if getattr(model, "add_D_objective", False):
            objective_values["D"] = float(model.D.getValue())
        solver_specific = {
            "gurobi_model": model.model,
            "gurobi_solution": model.solution,
        }
        return SolverResult(
            solution=model.solution.copy(),
            runtime=runtime,
            backend_name=self.name,
            gap_status=getattr(model, "gap_status", {}),
            raw_model=model.model,
            wrapper=model,
            objective_values=objective_values,
            solver_specific=solver_specific,
        )

    def find_alternative(self, allowed_complexity, reference, max_ci, max_d):
        from gurobipy import GRB, quicksum

        model_kwargs = {**self.config, "allowed_tree_complexity": allowed_complexity}
        model = GurobiModel(
            segment=self.inputs.segment,
            ci_table=self.inputs.ci_table,
            fractional_copy_number_table=self.inputs.fractional_copy_number_table,
            tree=self.inputs.tree,
            clone_proportions=self.inputs.clone_proportions,
            **model_kwargs,
        )
        g = model.model
        # no-good cut on integers: X - ref = up - down with up, down not both positive,
        # and at least one unit of difference over all referenced copy numbers
        differences = []
        for (allele, clone), value in reference.items():
            up = g.addVar(vtype=GRB.INTEGER, lb=0, name=f"alt_up_{allele}_{clone}")
            down = g.addVar(vtype=GRB.INTEGER, lb=0, name=f"alt_down_{allele}_{clone}")
            is_up = g.addVar(vtype=GRB.BINARY, name=f"alt_is_up_{allele}_{clone}")
            g.addConstr(model.X[allele][clone] - value == up - down)
            g.addConstr(up <= 1000 * is_up)
            g.addConstr(down <= 1000 * (1 - is_up))
            differences += [up, down]
        g.addConstr(quicksum(differences) >= 1, name="alt_differs")
        if model.add_CI_objective and max_ci is not None:
            g.addConstr(model.Z <= max_ci, name="alt_max_ci")
        if model.add_D_objective and max_d is not None:
            g.addConstr(model.D <= max_d, name="alt_max_d")
        g.optimize()
        if g.Status in (GRB.INFEASIBLE, GRB.INF_OR_UNBD):
            return None
        if g.Status != GRB.OPTIMAL:
            raise RuntimeError(f"alternative search ended with Gurobi status {g.Status}")
        model.get_output()
        return model.solution.copy()

    def supports_feature(self, feature_name: str) -> bool:
        if feature_name in {"iis", "objective_gap"}:
            return True
        return super().supports_feature(feature_name)

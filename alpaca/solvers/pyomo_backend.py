"""Pyomo-backed implementation of the ALPACA solver interface."""

from __future__ import annotations

from typing import Dict, List, Tuple

import pandas as pd
from pyomo import environ as pyo

from alpaca.utils import (
    find_path_edges,
    flat_list,
    get_length_from_name,
    get_tree_edges,
    SegmentInfeasibleError,
)

from .base import ModelInputs, SolverBackend, SolverResult


class PyomoBackend(SolverBackend):
    """Pyomo-based solver backend for ALPACA."""

    name = "pyomo"

    def __init__(self, inputs: ModelInputs, config: Dict[str, any]):
        super().__init__(inputs, config)
        self.inputs = inputs
        self.solver_name = config.get("pyomo_solver", "cbc")
        self.strict_gap = bool(config.get("strict_gap", True))
        self.time_limit = int(config.get("time_limit", 60))
        self.cpus = int(config.get("cpus", 2))
        self.custom_solver_options = dict(config.get("pyomo_solver_options", {}) or {})
        self.objectives = str(config.get("objectives", "DCI")).upper()
        self.add_D_objective = "D" in self.objectives
        self.add_CI_objective = "CI" in self.objectives
        if not (self.add_CI_objective or self.add_D_objective):
            raise ValueError(
                "At least one objective (D, CI) must be enabled for Pyomo backend"
            )
        self.minimise_events_to_diploid = bool(
            config.get("minimise_events_to_diploid", True)
        )
        self.exclusive_amp_del = bool(config.get("exclusive_amp_del", True))
        self.variability_penalty = int(config.get("variability_penalty", 0))
        self.enforce_tree_complexity = bool(config.get("enforce_tree_complexity", True))
        self.allowed_tree_complexity_default = int(
            config.get("allowed_tree_complexity", 1000)
        )
        self.prevent_increase_from_zero_flag = bool(
            config.get("prevent_increase_from_zero_flag", True)
        )
        self.add_event_count_constraints_flag = bool(
            config.get("add_event_count_constraints_flag", True)
        )
        self.add_state_change_count_constraints_flag = bool(
            config.get("add_state_change_count_constraints_flag", False)
        )
        self.add_allow_only_one_non_directional_event_flag = bool(
            config.get("add_allow_only_one_non_directional_event_flag", True)
        )
        self.limit_homozygous_deletions_threshold_flag = bool(
            config.get("limit_homozygous_deletions_threshold_flag", True)
        )
        self.restrict_heterogeneity_flag = bool(
            config.get("restrict_heterogeneity_flag", False)
        )
        self.restrict_to_clonal_only_flag = bool(
            config.get("restrict_to_clonal_only_flag", False)
        )
        self.homozygous_deletion_threshold = float(
            config.get("homozygous_deletion_threshold", 1)
        )
        self.homo_del_size_limit = int(config.get("homo_del_size_limit", 5 * 10**7))
        self.add_path_variability_penalty_constraints_flag = bool(
            config.get("add_path_variability_penalty_constraints_flag", False)
        )
        if self.variability_penalty not in (0, 0.0):
            raise NotImplementedError(
                "Pyomo backend currently requires variability_penalty == 0 to remain linear."
            )
        self.solver_log_path = config.get("solver_logs", "")
        self.output_directory = str(config.get("output_directory", ""))
        self.simulate_infeasibility = str(config.get("simulate_infeasibility", ""))
        self.big_m = 1000
        self.model: pyo.ConcreteModel | None = None
        self.allowed_tree_complexity = self.allowed_tree_complexity_default
        self.gap_status: Dict[str, float] = {}
        self.last_runtime: float | None = None
        self._init_internal_structures()

    def _init_internal_structures(self) -> None:
        self.alleles = ["A", "B"]
        self.clone_names = list(sorted(set(flat_list(self.inputs.tree))))
        self.tree_edges = list(sorted(get_tree_edges(self.inputs.tree)))
        self.tree_paths = list(self.inputs.tree)
        self.mrca = self.inputs.tree[0][0]
        if self.minimise_events_to_diploid:
            self.tree_edges.append(("diploid", self.mrca))
            if "diploid" not in self.clone_names:
                self.clone_names.append("diploid")
        self.sample_names = (
            self.inputs.fractional_copy_number_table["sample"].unique().tolist()
        )
        self.clone_proportions = self.inputs.clone_proportions.copy()
        self.clone_prop_lookup = self.clone_proportions.stack().to_dict()
        self.segment = self.inputs.segment
        self.tumour_id = self.inputs.fractional_copy_number_table["tumour_id"].iloc[0]
        self.Y = {
            allele: {
                row["sample"]: row[f"cpn{allele}"]
                for _, row in self.inputs.fractional_copy_number_table.iterrows()
            }
            for allele in self.alleles
        }
        self.ci_lookup: Dict[Tuple[str, str, str], float] = {}
        ci_df = self.inputs.ci_table[self.inputs.ci_table.segment == self.segment]
        for _, row in ci_df.iterrows():
            sample = row["sample"]
            for allele in self.alleles:
                self.ci_lookup[(allele, sample, "upper")] = float(
                    row[f"upper_CI_{allele}"]
                )
                self.ci_lookup[(allele, sample, "lower")] = float(
                    row[f"lower_CI_{allele}"]
                )
        self.path_edges_lookup = {
            idx: list(find_path_edges(path, self.tree_edges))
            for idx, path in enumerate(self.tree_paths)
        }

    def solve(self, allowed_complexity: int) -> SolverResult:
        self.allowed_tree_complexity = allowed_complexity
        self._build_model()
        if self.add_CI_objective and self.add_D_objective:
            self._activate_objective("CI")
            primary = self._run_solver(stage="CI")
            self._check_feasibility(primary)
            self._fix_ci_objective()
            secondary = self._run_solver(stage="D")
            results = secondary
        else:
            target = "CI" if self.add_CI_objective else "D"
            self._activate_objective(target)
            results = self._run_solver(stage=target)
            self._check_feasibility(results)
        solution_df = self._build_solution_frame()
        runtime = self._extract_runtime(results)
        if runtime is None:
            runtime = self.last_runtime
        objective_values = self._collect_objectives()
        return SolverResult(
            solution=solution_df,
            runtime=runtime,
            backend_name=self.name,
            gap_status=self.gap_status,
            raw_model=self.model,
            wrapper=self,
            objective_values=objective_values,
            solver_specific={"pyomo_results": results},
        )

    # ------------------------------------------------------------------
    # Model construction helpers
    # ------------------------------------------------------------------
    def _build_model(self) -> None:
        m = pyo.ConcreteModel()
        m.alleles = pyo.Set(initialize=self.alleles)
        m.clones = pyo.Set(initialize=self.clone_names)
        m.samples = pyo.Set(initialize=self.sample_names)
        m.tree_edges = pyo.Set(initialize=self.tree_edges, dimen=2)
        m.paths = pyo.Set(initialize=list(self.path_edges_lookup.keys()))
        self.model = m
        self._declare_variables()
        self._add_yhat_constraints()
        self._add_ci_constraints()
        self._add_absolute_distance_constraints()
        self._add_event_count_constraints()
        self._add_state_change_constraints()
        self._add_path_variability_penalty()
        if self.add_allow_only_one_non_directional_event_flag:
            self._add_allow_only_one_non_directional_event()
        if self.limit_homozygous_deletions_threshold_flag:
            self._limit_homozygous_deletions()
        if self.restrict_heterogeneity_flag:
            self._restrict_heterogeneity()
        if self.prevent_increase_from_zero_flag:
            self._prevent_increase_from_zero()
        if self.restrict_to_clonal_only_flag:
            self._restrict_to_clonal_only()
        self._set_complexity_constraints()
        self._build_objectives()
        if self.simulate_infeasibility and self.simulate_infeasibility == self.segment:
            print('Simulating infeasibility...')
            # Create a dummy variable fixed to 0, then require it to be >= 1 (infeasible)
            m._inf_dummy = pyo.Var(domain=pyo.NonNegativeIntegers, bounds=(0, 0))
            m.simulate_infeasibility_constr = pyo.Constraint(expr=m._inf_dummy >= 1)
            print(
                f"simulate_infeasibility: injected infeasible constraint for segment '{self.segment}'"
            )

    def _declare_variables(self) -> None:
        assert self.model is not None
        m = self.model
        m.X = pyo.Var(m.alleles, m.clones, domain=pyo.NonNegativeIntegers)
        m.Yhat = pyo.Var(m.alleles, m.samples, domain=pyo.NonNegativeReals)
        m.d = pyo.Var(m.alleles, m.samples, domain=pyo.NonNegativeReals)
        m.yhat_above_upper = pyo.Var(m.alleles, m.samples, domain=pyo.Binary)
        m.yhat_below_lower = pyo.Var(m.alleles, m.samples, domain=pyo.Binary)
        m.CI_overlap = pyo.Var(m.alleles, m.samples, domain=pyo.Binary)
        m.CN_diff_amp = pyo.Var(m.alleles, m.tree_edges, domain=pyo.NonNegativeIntegers)
        m.CN_diff_del = pyo.Var(m.alleles, m.tree_edges, domain=pyo.NonNegativeIntegers)
        m.cpn_change_up = pyo.Var(m.alleles, m.tree_edges, domain=pyo.Binary)
        m.cpn_change_down = pyo.Var(m.alleles, m.tree_edges, domain=pyo.Binary)
        m.n = pyo.Var(m.alleles, domain=pyo.NonNegativeIntegers)
        m.total_events_count = pyo.Var(domain=pyo.NonNegativeIntegers)
        m.total_events = pyo.Var(domain=pyo.NonNegativeIntegers, initialize=0)
        m.total_edge_changes = pyo.Var(domain=pyo.NonNegativeIntegers, initialize=0)
        m.total_edge_changes_count = pyo.Var(domain=pyo.NonNegativeIntegers)
        m.total_path_variability_penalty = pyo.Var(
            domain=pyo.NonNegativeIntegers, initialize=0
        )
        m.total_path_variability_penalty_count = pyo.Var(domain=pyo.NonNegativeIntegers)
        m.total_tree_complexity = pyo.Var(domain=pyo.NonNegativeReals)
        m.amps_count_on_path = pyo.Var(
            m.alleles, m.paths, domain=pyo.NonNegativeIntegers
        )
        m.dels_count_on_path = pyo.Var(
            m.alleles, m.paths, domain=pyo.NonNegativeIntegers
        )
        m.more_than_1_amp_change = pyo.Var(m.alleles, m.paths, domain=pyo.Binary)
        m.more_than_1_del_change = pyo.Var(m.alleles, m.paths, domain=pyo.Binary)
        m.path_variability_penalty = pyo.Var(
            m.alleles, m.paths, domain=pyo.NonNegativeIntegers
        )

    def _clone_prop(self, sample: str, clone: str) -> float:
        return float(self.clone_prop_lookup.get((clone, sample), 0.0))

    # ------------------------------------------------------------------
    # Constraint builders
    # ------------------------------------------------------------------
    def _add_yhat_constraints(self) -> None:
        assert self.model is not None
        m = self.model

        def yhat_rule(model, allele, sample):
            return model.Yhat[allele, sample] == sum(
                self._clone_prop(sample, clone) * model.X[allele, clone]
                for clone in self.clone_names
            )

        m.Yhat_constraint = pyo.Constraint(m.alleles, m.samples, rule=yhat_rule)

    def _add_ci_constraints(self) -> None:
        assert self.model is not None
        m = self.model
        m.ci_constraints = pyo.ConstraintList()
        for allele in self.alleles:
            for sample in self.sample_names:
                upper = self.ci_lookup.get((allele, sample, "upper"), 0.0)
                lower = self.ci_lookup.get((allele, sample, "lower"), 0.0)
                m.ci_constraints.add(
                    m.Yhat[allele, sample]
                    >= upper
                    - self.big_m
                    + self.big_m * m.yhat_above_upper[allele, sample]
                )
                m.ci_constraints.add(
                    m.Yhat[allele, sample]
                    <= upper + self.big_m * m.yhat_above_upper[allele, sample]
                )
                m.ci_constraints.add(
                    m.Yhat[allele, sample]
                    <= lower + self.big_m * (1 - m.yhat_below_lower[allele, sample])
                )
                m.ci_constraints.add(
                    m.Yhat[allele, sample]
                    >= lower - self.big_m * m.yhat_below_lower[allele, sample]
                )
                m.ci_constraints.add(
                    m.CI_overlap[allele, sample] >= m.yhat_above_upper[allele, sample]
                )
                m.ci_constraints.add(
                    m.CI_overlap[allele, sample] >= m.yhat_below_lower[allele, sample]
                )
                m.ci_constraints.add(
                    m.CI_overlap[allele, sample]
                    <= m.yhat_above_upper[allele, sample]
                    + m.yhat_below_lower[allele, sample]
                )

    def _add_absolute_distance_constraints(self) -> None:
        assert self.model is not None
        m = self.model
        m.abs_distance_constraints = pyo.ConstraintList()
        for allele in self.alleles:
            for sample in self.sample_names:
                obs = self.Y[allele][sample]
                m.abs_distance_constraints.add(
                    m.d[allele, sample] >= m.Yhat[allele, sample] - obs
                )
                m.abs_distance_constraints.add(
                    -m.d[allele, sample] <= m.Yhat[allele, sample] - obs
                )

    def _add_event_count_constraints(self) -> None:
        assert self.model is not None
        m = self.model

        def edge_balance_rule(model, allele, parent, child):
            return (
                model.X[allele, parent] + model.CN_diff_amp[allele, parent, child]
                == model.X[allele, child] + model.CN_diff_del[allele, parent, child]
            )

        m.edge_balance = pyo.Constraint(m.alleles, m.tree_edges, rule=edge_balance_rule)

        if self.exclusive_amp_del:
            m.exclusive_constraints = pyo.ConstraintList()
            for allele in self.alleles:
                for parent, child in self.tree_edges:
                    m.exclusive_constraints.add(
                        m.CN_diff_amp[allele, parent, child]
                        <= self.big_m * m.cpn_change_up[allele, parent, child]
                    )
                    m.exclusive_constraints.add(
                        m.CN_diff_amp[allele, parent, child]
                        >= m.cpn_change_up[allele, parent, child]
                    )
                    m.exclusive_constraints.add(
                        m.CN_diff_del[allele, parent, child]
                        <= self.big_m * m.cpn_change_down[allele, parent, child]
                    )
                    m.exclusive_constraints.add(
                        m.CN_diff_del[allele, parent, child]
                        >= m.cpn_change_down[allele, parent, child]
                    )
                    m.exclusive_constraints.add(
                        m.cpn_change_up[allele, parent, child]
                        + m.cpn_change_down[allele, parent, child]
                        <= 1
                    )

        def n_rule(model, allele):
            return model.n[allele] == sum(
                model.CN_diff_amp[allele, parent, child]
                + model.CN_diff_del[allele, parent, child]
                for (parent, child) in self.tree_edges
            )

        m.num_event_constraints = pyo.Constraint(m.alleles, rule=n_rule)
        m.total_events_count_constraint = pyo.Constraint(
            expr=m.total_events_count == sum(m.n[allele] for allele in self.alleles)
        )
        if self.add_event_count_constraints_flag:
            m.total_events_def = pyo.Constraint(
                expr=m.total_events == sum(m.n[allele] for allele in self.alleles)
            )
        else:
            m.total_events.fix(0)

    def _add_state_change_constraints(self) -> None:
        assert self.model is not None
        m = self.model
        m.state_change_constraints = pyo.ConstraintList()
        for allele in self.alleles:
            for parent, child in self.tree_edges:
                m.state_change_constraints.add(
                    m.CN_diff_amp[allele, parent, child]
                    <= self.big_m * m.cpn_change_up[allele, parent, child]
                )
                m.state_change_constraints.add(
                    m.CN_diff_amp[allele, parent, child]
                    >= m.cpn_change_up[allele, parent, child]
                )
                m.state_change_constraints.add(
                    m.CN_diff_del[allele, parent, child]
                    <= self.big_m * m.cpn_change_down[allele, parent, child]
                )
                m.state_change_constraints.add(
                    m.CN_diff_del[allele, parent, child]
                    >= m.cpn_change_down[allele, parent, child]
                )

        m.total_edge_changes_count_constraint = pyo.Constraint(
            expr=m.total_edge_changes_count
            == sum(
                m.cpn_change_up[allele, parent, child]
                + m.cpn_change_down[allele, parent, child]
                for allele in self.alleles
                for (parent, child) in self.tree_edges
            )
        )

        if self.add_state_change_count_constraints_flag:
            m.total_edge_changes_def = pyo.Constraint(
                expr=m.total_edge_changes
                == sum(
                    m.cpn_change_up[allele, parent, child]
                    + m.cpn_change_down[allele, parent, child]
                    for allele in self.alleles
                    for (parent, child) in self.tree_edges
                )
            )
        else:
            m.total_edge_changes.fix(0)

    def _add_path_variability_penalty(self) -> None:
        assert self.model is not None
        m = self.model
        # Since variability_penalty == 0, keep penalties at zero while retaining reporting variables
        for allele in self.alleles:
            for path_idx in m.paths:
                m.path_variability_penalty[allele, path_idx].fix(0)
        m.total_path_variability_penalty_count_constraint = pyo.Constraint(
            expr=m.total_path_variability_penalty_count
            == sum(
                m.path_variability_penalty[allele, path_idx]
                for allele in self.alleles
                for path_idx in m.paths
            )
        )
        if self.add_path_variability_penalty_constraints_flag:
            m.total_path_variability_penalty_def = pyo.Constraint(
                expr=m.total_path_variability_penalty
                == sum(
                    m.path_variability_penalty[allele, path_idx]
                    for allele in self.alleles
                    for path_idx in m.paths
                )
            )
        else:
            m.total_path_variability_penalty.fix(0)

    def _add_allow_only_one_non_directional_event(self) -> None:
        assert self.model is not None
        m = self.model
        U = self.big_m
        m.path_direction_constraints = pyo.ConstraintList()
        for allele in self.alleles:
            for path_idx in m.paths:
                path_edges = self.path_edges_lookup[path_idx]
                m.path_direction_constraints.add(
                    m.amps_count_on_path[allele, path_idx]
                    == sum(
                        m.cpn_change_up[allele, parent, child]
                        for parent, child in path_edges
                    )
                )
                m.path_direction_constraints.add(
                    m.dels_count_on_path[allele, path_idx]
                    == sum(
                        m.cpn_change_down[allele, parent, child]
                        for parent, child in path_edges
                    )
                )
                m.path_direction_constraints.add(
                    m.amps_count_on_path[allele, path_idx]
                    >= 2 - U + U * m.more_than_1_amp_change[allele, path_idx]
                )
                m.path_direction_constraints.add(
                    m.amps_count_on_path[allele, path_idx]
                    <= 1 + U * m.more_than_1_amp_change[allele, path_idx]
                )
                m.path_direction_constraints.add(
                    m.dels_count_on_path[allele, path_idx]
                    >= 2 - U + U * m.more_than_1_del_change[allele, path_idx]
                )
                m.path_direction_constraints.add(
                    m.dels_count_on_path[allele, path_idx]
                    <= 1 + U * m.more_than_1_del_change[allele, path_idx]
                )
                m.path_direction_constraints.add(
                    m.more_than_1_amp_change[allele, path_idx]
                    + m.more_than_1_del_change[allele, path_idx]
                    <= 1
                )

    def _prevent_increase_from_zero(self) -> None:
        assert self.model is not None
        m = self.model
        m.prevent_zero_constraints = pyo.ConstraintList()
        for allele in self.alleles:
            for parent, child in self.tree_edges:
                m.prevent_zero_constraints.add(
                    m.X[allele, parent] >= m.cpn_change_up[allele, parent, child]
                )

    def _limit_homozygous_deletions(self) -> None:
        assert self.model is not None
        m = self.model
        threshold = float(self.homozygous_deletion_threshold)
        seg_len = get_length_from_name(self.segment)
        if seg_len >= self.homo_del_size_limit:
            for clone in self.clone_names:
                if clone == "diploid":
                    continue
                m.ci_constraints.add(m.X["A", clone] + m.X["B", clone] >= 1)
            return

        samples_with_low_A = [
            sample for sample, value in self.Y["A"].items() if value < threshold
        ]
        samples_with_low_B = [
            sample for sample, value in self.Y["B"].items() if value < threshold
        ]
        permitted_samples = sorted(set(samples_with_low_A) & set(samples_with_low_B))
        if not permitted_samples:
            for clone in self.clone_names:
                if clone == "diploid":
                    continue
                m.ci_constraints.add(m.X["A", clone] + m.X["B", clone] >= 1)
            return

        cp_table = self.clone_proportions
        valid_samples = [
            sample for sample in permitted_samples if sample in cp_table.columns
        ]
        if not valid_samples:
            valid_samples = []
        clones_in_samples = cp_table.index[
            cp_table[valid_samples].sum(axis=1) > 0
        ].tolist()
        clones_not_present = [
            clone for clone in self.clone_names if clone not in clones_in_samples
        ]
        for clone in clones_not_present:
            if clone == "diploid":
                continue
            m.ci_constraints.add(m.X["A", clone] + m.X["B", clone] >= 1)

    def _restrict_heterogeneity(self) -> None:
        assert self.model is not None
        m = self.model
        ci_table = self.inputs.ci_table[
            self.inputs.ci_table.segment == self.segment
        ].copy()
        for clone in self.clone_names:
            if clone == "diploid":
                continue
            if clone not in self.clone_proportions.index:
                continue
            samples_present = (self.clone_proportions.loc[clone] > 0).astype(bool)
            sample_names = samples_present[samples_present].index.tolist()
            clone_ci = ci_table[ci_table["sample"].isin(sample_names)]
            clone_ci["span_A"] = clone_ci[f"upper_CI_A"].astype(int) - clone_ci[
                f"lower_CI_A"
            ].astype(int)
            clone_ci["span_B"] = clone_ci[f"upper_CI_B"].astype(int) - clone_ci[
                f"lower_CI_B"
            ].astype(int)
            if (clone_ci["span_A"] == 1).any():
                cn_val = int(clone_ci[f"upper_CI_A"][clone_ci["span_A"] == 1].median())
                m.ci_constraints.add(m.X["A", clone] == cn_val)
            if (clone_ci["span_B"] == 1).any():
                cn_val = int(clone_ci[f"upper_CI_B"][clone_ci["span_B"] == 1].median())
                m.ci_constraints.add(m.X["B", clone] == cn_val)

    def _restrict_to_clonal_only(self) -> None:
        assert self.model is not None
        m = self.model
        clonal_clones = [
            clone
            for clone in self.clone_names
            if clone != "diploid" and clone in self.clone_proportions.index
        ]
        num_clones = len(clonal_clones)
        if num_clones == 0:
            return
        for allele in self.alleles:
            total = sum(m.X[allele, clone] for clone in clonal_clones)
            for clone in clonal_clones:
                m.ci_constraints.add(num_clones * m.X[allele, clone] == total)

    def _set_complexity_constraints(self) -> None:
        assert self.model is not None
        m = self.model
        m.total_complexity_components = pyo.Constraint(
            expr=m.total_tree_complexity
            == m.total_path_variability_penalty + m.total_events + m.total_edge_changes
        )
        if self.enforce_tree_complexity:
            m.tree_complexity_limit = pyo.Constraint(
                expr=m.total_tree_complexity <= self.allowed_tree_complexity
            )

    def _build_objectives(self) -> None:
        assert self.model is not None
        m = self.model
        m.Z_expr = pyo.Expression(
            expr=sum(
                m.CI_overlap[allele, sample]
                for allele in self.alleles
                for sample in self.sample_names
            )
        )
        m.D_expr = pyo.Expression(
            expr=sum(
                m.d[allele, sample]
                for allele in self.alleles
                for sample in self.sample_names
            )
        )
        if self.add_CI_objective:
            m.obj_CI = pyo.Objective(expr=m.Z_expr, sense=pyo.minimize)
        if self.add_D_objective:
            m.obj_D = pyo.Objective(expr=m.D_expr, sense=pyo.minimize)
            if self.add_CI_objective:
                m.obj_D.deactivate()

    def _activate_objective(self, name: str) -> None:
        assert self.model is not None
        if name == "CI" and getattr(self.model, "obj_CI", None) is not None:
            self.model.obj_CI.activate()
            if getattr(self.model, "obj_D", None) is not None:
                self.model.obj_D.deactivate()
        elif name == "D" and getattr(self.model, "obj_D", None) is not None:
            self.model.obj_D.activate()
            if getattr(self.model, "obj_CI", None) is not None:
                self.model.obj_CI.deactivate()

    # ------------------------------------------------------------------
    # Solver orchestration
    # ------------------------------------------------------------------
    def _run_solver(self, stage: str):
        assert self.model is not None
        # On macOS, SCIP can be slow to respond to version queries during solver detection.
        # Set environment variable to increase timeout for solver factory initialization.
        solver_name_lower = (self.solver_name or "").lower()
        if solver_name_lower in {"scip", "scipampl"}:
            import os

            # Set environment variable for Pyomo's subprocess timeout handling
            os.environ.setdefault("PYOMO_SOLVER_EXEC_TIMEOUT", "30")

        solver = pyo.SolverFactory(self.solver_name)
        options = self._solver_options()
        for key, value in options.items():
            try:
                solver.options[key] = value
            except Exception:
                pass
        results = solver.solve(self.model, tee=False)
        runtime = self._extract_runtime(results)
        self.last_runtime = runtime
        self._update_gap_status(results, stage)
        return results

    def _extract_runtime(self, results):
        runtime = getattr(results.solver, "time", None)
        if runtime in (None, ""):
            runtime = getattr(results.solver, "wallclock_time", None)
        if runtime not in (None, ""):
            try:
                runtime = float(runtime)
            except (TypeError, ValueError):
                runtime = None
        return runtime

    def _solver_options(self) -> Dict[str, float]:
        options: Dict[str, float] = {}
        solver_name = (self.solver_name or "").lower()
        if self.time_limit:
            if solver_name in {"cbc"}:
                options["seconds"] = self.time_limit
            elif solver_name in {"glpk"}:
                options["tmlim"] = self.time_limit
            else:
                options["TimeLimit"] = self.time_limit
        if self.strict_gap:
            if solver_name in {"cbc"}:
                options["ratioGap"] = 0
            elif solver_name in {"gurobi"}:
                options["MIPGap"] = 0.0
                options["MIPGapAbs"] = 0.0
        if self.cpus:
            if solver_name in {"gurobi"}:
                options["Threads"] = self.cpus
        if solver_name in {"scip", "scipampl"}:
            # on MacOS scipampl can throw an error with shorter timeout
            options.setdefault("ampl_command_timeout", 10)
        if self.custom_solver_options:
            options.update(self.custom_solver_options)
        return options

    def _check_feasibility(self, results) -> None:
        """Raise SegmentInfeasibleError if the solver returned an infeasible status."""
        termination = str(getattr(results.solver, "termination_condition", "")).lower()
        solver_status = str(getattr(results.solver, "status", "unknown"))
        if "infeasible" in termination:
            # Write a JSON summary report so it can be combined with other segments later.
            import json as _json
            import os as _os
            from datetime import datetime as _dt

            report_dir = (
                getattr(self, "output_directory", "")
                or getattr(self, "solver_log_path", "")
                or _os.getcwd()
            )
            _os.makedirs(report_dir, exist_ok=True)
            json_path = _os.path.join(
                report_dir,
                f"{self.tumour_id}_{self.segment}_infeasibility_report.json",
            )
            report_data = {
                "tumour_id": self.tumour_id,
                "segment": self.segment,
                "solver_status": solver_status,
                "termination_condition": termination,
                "timestamp": _dt.now().isoformat(),
                "solver": self.solver_name,
            }
            try:
                with open(json_path, "w") as fh:
                    _json.dump(report_data, fh, indent=4)
                print(f"Wrote infeasibility JSON report: {json_path}")
            except Exception as e:
                print(f"Failed to write infeasibility JSON report: {e}")
            print(
                f"Solver returned infeasible (termination={termination}). "
                f"Skipping segment '{self.segment}'."
            )
            raise SegmentInfeasibleError(self.tumour_id, self.segment, solver_status)

    def _update_gap_status(self, results, stage: str) -> None:
        gap = getattr(results.solver, "relative_gap", None)
        runtime = getattr(results.solver, "time", None)
        termination = str(getattr(results.solver, "termination_condition", "unknown"))
        status = str(getattr(results.solver, "status", "unknown"))
        self.gap_status = {
            "max_gap": gap if gap is not None else -1,
            "gap_reason": termination,
            "runtime": runtime if runtime is not None else self.last_runtime,
            "status": status,
            "stage": stage,
        }

    def _fix_ci_objective(self) -> None:
        assert self.model is not None and hasattr(self.model, "Z_expr")
        best_ci = pyo.value(self.model.Z_expr)
        tolerance = 1e-6
        self.model.ci_objective_fix = pyo.Constraint(
            expr=self.model.Z_expr <= best_ci + tolerance
        )
        self._activate_objective("D")

    # ------------------------------------------------------------------
    # Output helpers
    # ------------------------------------------------------------------
    def _collect_objectives(self) -> Dict[str, float]:
        assert self.model is not None
        values: Dict[str, float] = {}
        if self.add_CI_objective:
            values["CI"] = float(pyo.value(self.model.Z_expr))
        if self.add_D_objective:
            values["D"] = float(pyo.value(self.model.D_expr))
        return values

    def _build_solution_frame(self) -> pd.DataFrame:
        assert self.model is not None
        records: List[Dict[str, float]] = []
        for clone in self.clone_names:
            if clone == "diploid":
                continue
            record = {
                "clone": clone,
                "pred_CN_A": int(round(pyo.value(self.model.X["A", clone]))),
                "pred_CN_B": int(round(pyo.value(self.model.X["B", clone]))),
            }
            records.append(record)
        solution = pd.DataFrame.from_records(records)
        ci_score = (
            int(round(pyo.value(self.model.Z_expr))) if self.add_CI_objective else 0
        )
        d_score = float(pyo.value(self.model.D_expr)) if self.add_D_objective else 0.0
        solution["complexity"] = int(round(pyo.value(self.model.total_tree_complexity)))
        solution["CI_score"] = ci_score
        solution["D_score"] = round(d_score, 3)
        solution["variability_penalty_count"] = int(
            round(pyo.value(self.model.total_path_variability_penalty_count))
        )
        solution["state_change_count"] = int(
            round(pyo.value(self.model.total_edge_changes_count))
        )
        solution["event_count"] = int(round(pyo.value(self.model.total_events_count)))
        solution["allowed_complexity"] = self.allowed_tree_complexity
        solution["gurobi_time_CI"] = -1
        solution["gurobi_gap_CI"] = -1
        solution["gurobi_time_D"] = -1
        solution["gurobi_gap_D"] = -1
        solution = solution.sort_values("clone").reset_index(drop=True)
        solution.index.name = "index"
        solution.reset_index(drop=True, inplace=True)
        return solution

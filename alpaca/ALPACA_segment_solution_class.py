import os
from datetime import datetime
import pandas as pd
import json
import numpy as np
import math
import kneed
from scipy.stats import norm
import typing
from typing import Optional, Dict, Any
import time
from alpaca.solvers import ModelInputs, SolverResult, create_solver_backend
from alpaca.utils import read_tree_json, save_dataframe, SegmentInfeasibleError
from alpaca.edge_cases import (
    AMBIGUITY_CHECK_FAILED,
    AMBIGUOUS,
    SINGLE_CLONE,
    SINGLE_CLONE_AND_SINGLE_SAMPLE,
    SINGLE_SAMPLE,
    STANDARD,
    classify_tumour,
    make_segment_report,
    segment_report_name,
    solve_single_clone,
)
import logging
try:
    import matplotlib.pyplot as plt  # type: ignore
except Exception:
    plt = None


def ensure_elbow_strictly_decreasing(df):
    for col in ["D_score"]:
        df[col] = df[col].round(3)
        current_minimum = 1000
        new_values = []
        for v in df[col]:
            if v < current_minimum:
                new_values.append(v)
                current_minimum = v
            else:
                new_v = current_minimum - 0.001
                current_minimum = new_v
                new_values.append(new_v)
        df[col] = new_values
    return df


def find_s_values(elbow_search_df, max_iterations, x="allowed_complexity", y="D_score"):
    elbow_table = (elbow_search_df[[x, y]]).dropna()
    m = 1
    s = kneed.KneeLocator(
        elbow_table[x],
        elbow_table[y],
        S=m,
        curve="convex",
        direction="decreasing",
        interp_method="interp1d",
        online=True,
    ).knee
    s_code = "default"
    if not s:
        # try with higher sensitivity:
        s_candidates = list(
            set(
                [
                    s
                    for s in [
                        kneed.KneeLocator(
                            elbow_table[x],
                            elbow_table[y],
                            S=m,
                            curve="convex",
                            direction="decreasing",
                            interp_method="interp1d",
                            online=True,
                        ).knee
                        for m in range(0, 200)
                    ]
                    if s is not None
                ]
            )
        )
        if len(s_candidates) > 0:
            s = s_candidates[0]
            s_code = "high_sensitivity"
    return s, s_code


def missing_clones_inherit_from_children(optimal_solution, tree, cp_table):
    """
    Only applies when a missing clone has exactly one child, Otherwise parsimony should be enough to find the correct solution.
    """
    missing_clones = list(cp_table.loc[cp_table.sum(1) == 0].index)
    # exclude leaves as they cannot have any children:
    leaves = [x[-1] for x in tree]
    missing_clones = [x for x in missing_clones if x not in leaves]
    if len(missing_clones) > 0:
        for clone in missing_clones:
            tree_branch = [branch for branch in tree if clone in branch][0]
            missing_in_branch = [c for c in tree_branch if c in missing_clones]
            missing_in_branch.reverse()
            for missing_clone in missing_in_branch:
                non_missing_children = []
                for branch in tree:
                    if missing_clone in branch:
                        child_index = branch.index(missing_clone) + 1
                        child = branch[child_index]
                        if child not in non_missing_children:
                            non_missing_children.append(child)
                if len(non_missing_children) == 1:
                    non_missing_child = non_missing_children[0]
                    optimal_solution.loc[
                        optimal_solution.clone == missing_clone, "pred_CN_A"
                    ] = optimal_solution.loc[
                        optimal_solution.clone == non_missing_child, "pred_CN_A"
                    ].values[
                        0
                    ]
                    optimal_solution.loc[
                        optimal_solution.clone == missing_clone, "pred_CN_B"
                    ] = optimal_solution.loc[
                        optimal_solution.clone == non_missing_child, "pred_CN_B"
                    ].values[
                        0
                    ]
    return optimal_solution


def coerce_bool_like_series(series: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(series):
        return series.astype("boolean")
    if pd.api.types.is_numeric_dtype(series):
        mapped = series.map({1: True, 0: False})
        mapped[series.isna()] = pd.NA
        return mapped.astype("boolean")
    normalized = series.astype(str).str.strip().str.lower()
    mapped = normalized.map(
        {
            "true": True,
            "false": False,
            "1": True,
            "0": False,
            "yes": True,
            "no": False,
            "y": True,
            "n": False,
            "t": True,
            "f": False,
        }
    )
    mapped[series.isna()] = pd.NA
    return mapped.astype("boolean")


def remove_small_clones(cp_table, tree):
    tree_levels = max([len(branch) for branch in tree])
    MRCA = tree[0][0]
    for level in reversed(range(0, tree_levels)):
        for branch in tree:
            try:
                clone = branch[level]
                for region in cp_table.columns:
                    clone_cp = cp_table.loc[clone, region]
                    if (clone_cp > 0) & (clone_cp < 0.1) & (clone != MRCA):
                        parent = branch[branch.index(clone) - 1]
                        cp_table.loc[parent, region] = (
                            cp_table.loc[parent, region] + clone_cp
                        )
                        cp_table.loc[clone, region] = 0
            except IndexError:
                pass
    return cp_table


def split_input_file_name(input_file_name: str, df:pd.DataFrame):
    stripped_name = input_file_name.split("ALPACA_input_table_")[1]
    if stripped_name.count("_") != 3:
        # the segment file name is probably malformed. Try to read the values from the file:
        try:
            tumour_id_all = df.tumour_id.unique()
            segment_all = df.segment.unique()
            assert len(tumour_id_all) == 1, "Multiple tumour ids found in the input file, cannot determine tumour id for this segment"
            assert len(segment_all) == 1, "Multiple segments found in the input file, cannot determine segment name for this segment"
            tumour_id = tumour_id_all[0]
            segment = segment_all[0]
            return tumour_id, segment
        except Exception as e:
            # return original error:
            print(e)
            raise ValueError(
                f"Input file name {input_file_name} is not in the expected format and values could not be read from the file. Expected format: ALPACA_input_table_<tumour_id>_<segment>.<ext>"
            )
    stripped_name = stripped_name.replace(".parquet", "").replace(".csv", "")
    t_id = stripped_name.split("_")[0]
    s_name = "_".join(stripped_name.split("_")[1:])
    return t_id, s_name


def calculate_CI(df_seg_reg, CI):
    for allele in ["A", "B"]:
        data = np.array(df_seg_reg[f"ph_cpn{allele}_vec"])
        # drop nans:
        data = data[~np.isnan(data)]
        m = np.mean(data)
        std = np.std(data, ddof=1)
        lower_CI = norm.ppf(CI / 2, loc=m, scale=std)
        upper_CI = norm.ppf(1 - CI / 2, loc=m, scale=std)
        df_seg_reg[f"lower_CI_{allele}"] = lower_CI
        df_seg_reg[f"upper_CI_{allele}"] = upper_CI
    return df_seg_reg


def get_ci_table(input_table, tumour_dir, segment, ci_table_name="", CI=0.5, min_ci=None):
    # if confidence interval table name is provided, read it:
    if ci_table_name != "":
        ci_table = pd.read_csv(f"{tumour_dir}/{ci_table_name}")
        ci_table = ci_table[ci_table.segment == segment]
        # Enforce minimum CI span if caller provided min_ci via preprocessing_config
        affected_samples = set()
        affected_alleles = set()
        affected_segment = segment
        if min_ci is not None:
            min_ci = float(min_ci)
            required_cols = ["lower_CI_A", "upper_CI_A", "lower_CI_B", "upper_CI_B"]
            if all([c in ci_table.columns for c in required_cols]):
                for allele in ["A", "B"]:
                    lower_col = f"lower_CI_{allele}"
                    upper_col = f"upper_CI_{allele}"
                    span_col = f"span_CI_{allele}"
                    mid_col = f"mid_CI_{allele}"
                    # compute original span
                    ci_table[span_col] = ci_table[upper_col] - ci_table[lower_col]
                    ci_table[mid_col] = (ci_table[upper_col] + ci_table[lower_col]) / 2
                    # if span smaller than min_ci, expand symmetrically around midpoint
                    too_small = ci_table[span_col] < min_ci
                    if too_small.any():
                        affected_rows = ci_table.loc[too_small]
                        for s in affected_rows["sample"].unique():
                            affected_samples.add(s)
                        affected_alleles.add(allele)
                        ci_table.loc[too_small, upper_col] = (
                            ci_table.loc[too_small, mid_col] + min_ci / 2
                        )
                        ci_table.loc[too_small, lower_col] = (
                            ci_table.loc[too_small, mid_col] - min_ci / 2
                        )
                    # ensure lower bound non-negative
                    ci_table[lower_col] = ci_table[lower_col].apply(lambda x: max(x, 0))
                
                if affected_samples:
                    report = {
                        "tumour_id": os.path.basename(tumour_dir),
                        "segment": affected_segment,
                        "affected_samples": sorted(list(affected_samples)),
                        "affected_alleles": sorted(list(affected_alleles)),
                        "min_ci": float(min_ci),
                        "timestamp": datetime.now().isoformat(),
                    }
                else:
                    report = {}
                                
    # if table is not provided, but SNP table exists, calculate CI from SNP table:
    else:
        report = {}
        # try to read SNP table, if not present, create artificial CI table:
        try:
            asas_table = pd.read_csv(f"{tumour_dir}/asas_table.csv")
            asas_table = asas_table[asas_table.segment == segment]
            ci_table = asas_table.groupby(["sample", "segment"]).apply(
                lambda df_seg_reg: calculate_CI(df_seg_reg, CI)
            )
            ci_table = ci_table[
                [
                    "sample",
                    "segment",
                    "lower_CI_A",
                    "upper_CI_A",
                    "lower_CI_B",
                    "upper_CI_B",
                ]
            ].drop_duplicates()
        except FileNotFoundError:
            # create dummy CIs:
            print("No SNP table found, creating artificial CI table")
            ci_table = input_table[["sample", "segment"]].drop_duplicates().copy()
            for x in ["lower_CI_A", "upper_CI_A", "lower_CI_B", "upper_CI_B"]:
                ci_table[x] = float(0)  # to ensure float data type
            for s in ci_table["sample"].unique():
                A = input_table[input_table["sample"] == s].cpnA.median()
                B = input_table[input_table["sample"] == s].cpnB.median()
                ci_table.loc[ci_table["sample"] == s, "lower_CI_A"] = A - 0.5
                ci_table.loc[ci_table["sample"] == s, "lower_CI_B"] = B - 0.5
                ci_table.loc[ci_table["sample"] == s, "upper_CI_A"] = A + 0.5
                ci_table.loc[ci_table["sample"] == s, "upper_CI_B"] = B + 0.5
        for allele in ["A", "B"]:
            ci_table[f"lower_CI_{allele}"] = ci_table[f"lower_CI_{allele}"].apply(
                lambda x: max(x, 0)
            )
            ci_table[f"upper_CI_{allele}"] = ci_table[f"upper_CI_{allele}"].apply(
                lambda x: max(x, 0.01)
            )
        ci_table["ci_value"] = CI
        ci_table = ci_table.reset_index(drop=True).sort_values("sample")
    return ci_table, report


def validate_inputs(
    it: pd.DataFrame, cpt: pd.DataFrame, cit: pd.DataFrame, t: typing.List[typing.List]
):
    # check if tree is a list of lists of strings:
    # e.g. tree=[['cloneA','cloneB'],['cloneA','cloneD','cloneE']]
    if not all([isinstance(x, list) for x in t]):
        raise ValueError("Tree is not a list of lists of strings (clone names)")
    # check if all clones are present:
    cpt_clones = set(cpt.index.unique())
    tree_clones = set([c for branch in t for c in branch])
    if cpt_clones != tree_clones:
        raise ValueError("Clones in cp_table and tree_paths file (json or nwk) do not match")
    # check if all samples are present:
    it_samples = set(it["sample"].unique())
    cpt_samples = set(cpt.columns)
    cit_samples = set(cit["sample"].unique())
    if (
        (it_samples != cpt_samples)
        or (cpt_samples != cit_samples)
        or (it_samples != cit_samples)
    ):
        error_msg = f"""
        Sample names in input table: {sorted(list(it_samples))}
        Sample names in clone proportions table: {sorted(list(cpt_samples))}
        Sample names in confidence intervals table: {sorted(list(cit_samples))}"""
        raise ValueError(
            f"Sample names in input table, cp_table and ci_table do not match\n{error_msg}"
        )
    # check if segment is present in the ci_table:
    it_segments = set(it["segment"].unique())
    cit_segments = set(cit["segment"].unique())
    if not it_segments.issubset(cit_segments):
        raise ValueError("Segments in input table and ci_table do not match")
    # check if all columns are present in the input table:
    expected_columns = ["sample", "cpnA", "cpnB", "segment", "tumour_id"]
    if not set(expected_columns).issubset(set(it.columns)):
        raise ValueError(
            f"""Input table does not contain all expected columns.
            Columns in input table:
            {sorted(list(it.columns))}
            Expected columns:
            {sorted(expected_columns)}"""
        )
    # check if clone proportions sum to 1 for each sample
    proportions_expressed_as_percents = (cpt.sum() > 10).any()
    if proportions_expressed_as_percents:
        raise ValueError(
            "Clone proportions are probably expressed as percents, not fractions (e.g. 80 instead of 0.8)"
        )
    proportions_dont_sum_to_1 = (cpt.sum() != 1).any()
    if proportions_dont_sum_to_1:
        print("------WARNING------")
        print("Clone proportions do not sum to 1 in some samples")
        sum_df = (
            cpt.sum()
            .reset_index()
            .rename(columns={"index": "sample", 0: "proportions"})
        )
        print(sum_df)
        if (abs(cpt.sum() - 1) < 0.05).all():
            print(
                "Clones proportions are close to 1, calibrating them to sum to 1 (likely rounding errors)"
            )
            cpt = calibrate_clone_proportions(cpt)
        else:
            print("Clones proportions are not close to 1, exiting")
            proportions_below_1_in_any_sample = (cpt.sum() < 1).any()
            if proportions_below_1_in_any_sample:
                raise ValueError("Clone proportions sum to less than 1 in some samples")
            proportions_above_1_in_any_sample = (cpt.sum() > 1).any()
            if proportions_above_1_in_any_sample:
                raise ValueError("Clone proportions sum to more than 1 in some samples")
            
            
def calibrate_clone_proportions(cp: pd.DataFrame):
    for r in cp.select_dtypes(include=float).columns:
        cp[r] = cp[r] / cp[r].sum()
    return cp


def rescale_elbow_points(complexities, elbow):
    def rescale_below(arr, elbow):
        min_val = min(arr)
        max_val = elbow
        return [round((x - min_val) / (max_val - min_val) - 1, 2) for x in arr]

    def rescale_above(arr, elbow):
        min_val = elbow
        max_val = max(arr)
        return [round((x - min_val) / (max_val - min_val), 2) for x in arr]

    below_elbow = [x for x in complexities if x < elbow]
    above_elbow = [x for x in complexities if x > elbow]
    rescaled_below = rescale_below(below_elbow, elbow) if len(below_elbow) > 0 else []
    rescaled_above = rescale_above(above_elbow, elbow) if len(above_elbow) > 0 else []
    rescaled = rescaled_below + [0] + rescaled_above
    rescaled_dict = {k: v for k, v in zip(complexities, rescaled)}
    return rescaled_dict


class SegmentSolution:
    def __init__(
        self,
        input_file_name: str,
        config: Optional[Dict[str, Any]] = None,
        logger: Optional[logging.Logger] = None,
    ):
        self.logger = logger
        # get start time:
        self.start_time = time.time()
        if config is None:
            config = {"preprocessing_config": {}, "model_config": {}}
        # Define default values:
        self.ci: float = 0.5
        self.rsc: bool = False
        self.input_data_directory: Optional[str] = None
        self.ccp: bool = True
        self.s_type: str = "s_strictly_decreasing"
        self.missing_clones_inherit_from_children_flag: bool = True
        self.d_zero: int = 0
        self.optimal_solution: Optional[Any] = None
        self.optimal_solution_index: Optional[int] = None
        self.elbow: Dict[str, Any] = {}
        self.maximum_complexity: Optional[int] = None
        self.solutions_combined: pd.DataFrame = pd.DataFrame()
        self.config: Dict[str, Any] = config
        self.metrics: Dict[str, list] = {
            name: []
            for name in [
                "D_scores",
                "solutions",
                "run_time",
                "models",
                "complexity",
                "gap_status",
                "time_limits",
            ]
        }
        self.no_change_in_complexity: bool = False
        self.no_change_in_D_score: bool = False
        self.no_improvement_in_D_score: bool = False
        self.diploid_solution_found: bool = False
        self.compare_with_true_solution: bool = False
        self.ci_table_name: str = ""
        self.output_all_solutions: bool = False
        self.write_elbow_plots: bool = False
        self.output_model_selection_table: bool = False
        self.output_format: str = "parquet"
        self.debug: bool = False
        self.extra_columns: list = []
        self.used_clonal_shortcut: bool = False
        self.elbow_increase_report: pd.DataFrame = pd.DataFrame()
        # load config
        # default values present in the config object will overwrite the default values defined above
        for key, value in self.config["preprocessing_config"].items():
            setattr(self, key, value)
        for key, value in self.config["model_config"].items():
            setattr(self, key, value)
        self.input_file_name = input_file_name
        # define tumour input directory depending on the run environment:
        self.set_directories()
        # load fractional copy numbers:
        self.input_table = pd.read_csv(
            f"{self.segments_dir}/{input_file_name}"
        ).sort_values("sample")
        self.tumour_id, self.segment = split_input_file_name(self.input_file_name, self.input_table)
        # load tree:
        self.tree = read_tree_json(f"{self.tumour_dir}/tree_paths.json")
        # load clone proportions:
        self.cp_table = pd.read_csv(
            f"{self.tumour_dir}/cp_table.csv", index_col="clone"
        )
        self.cp_table = (
            calibrate_clone_proportions(self.cp_table) if self.ccp else self.cp_table
        )
        self.cp_table = (
            remove_small_clones(self.cp_table, self.tree) if self.rsc else self.cp_table
        )
        # get confidence intervals for copy number values:
        self.ci_table, self.ci_report = get_ci_table(
            self.input_table,
            self.tumour_dir,
            self.segment,
            ci_table_name=self.ci_table_name,
            CI=self.ci,
            min_ci=getattr(self, 'min_ci', None),
        )
        # check if inputs are in the expected format and contain all required columns:
        validate_inputs(
            it=self.input_table, cpt=self.cp_table, cit=self.ci_table, t=self.tree
        )
        # single clone / single sample tumours get dedicated handling (see alpaca.edge_cases):
        self.tumour_category = classify_tumour(self.cp_table)
        self.solver_path = "ilp"
        self.edge_case_flags: list = []
        self.edge_case_details: str = ""
        if self.tumour_category == SINGLE_SAMPLE and self.config["model_config"].get(
            "single_sample_restrict_heterogeneity", True
        ):
            # copy so the change cannot leak into other segments sharing this config
            self.config = {
                **self.config,
                "model_config": {**self.config["model_config"], "restrict_heterogeneity_flag": True},
            }
        self.solver_name = self.config["model_config"].get("solver", "pyomo")
        self.model_inputs = ModelInputs(
            segment=self.segment,
            ci_table=self.ci_table,
            fractional_copy_number_table=self.input_table,
            tree=self.tree,
            clone_proportions=self.cp_table,
        )
        self.get_monoclonal_samples_report()
        #
        print(datetime.now())
        print(f"Running: {input_file_name}")
        print(f"Tumour id: {self.tumour_id}")
        print(f"Segment name: {self.segment}")
        print("input table:")
        print(self.input_table)

    def get_model_metrics(self, solver_result):
        solution_df = solver_result.solution
        self.metrics["D_scores"].append(solution_df.D_score.iloc[0])
        self.metrics["solutions"].append(solution_df)
        self.metrics["run_time"].append(solver_result.runtime)
        self.metrics["models"].append(solver_result.raw_model)
        self.metrics["complexity"].append(solution_df.complexity.iloc[0])
        self.metrics["time_limits"].append(self.config["model_config"].get("time_limit", 60))
        gap_status = solver_result.gap_status or {
            "max_gap": -1,
            "gap_reason": "unknown",
            "runtime": solver_result.runtime or -1,
            "status": -1,
        }
        self.metrics["gap_status"].append(gap_status)

    def run_model(self, allowed_complexity, report_infeasibility=True):
        backend_config = {
            **self.config["model_config"],
            "allowed_tree_complexity": allowed_complexity,
            "report_infeasibility": report_infeasibility,
        }
        # Forward output_directory so infeasibility reports are written to the correct location.
        if "output_directory" not in backend_config:
            backend_config["output_directory"] = self.config["preprocessing_config"].get("output_directory", ".")
        solver_name = backend_config.get("solver", "pyomo")
        solver_logs = backend_config.get("solver_logs")
        if solver_name == "gurobi" and solver_logs and not backend_config.get("gurobi_logs"):
            backend_config["gurobi_logs"] = solver_logs
        backend = create_solver_backend(solver_name, self.model_inputs, backend_config)
        solver_result = backend.solve(allowed_complexity=allowed_complexity)
        if self.missing_clones_inherit_from_children_flag:
            solver_result.solution = missing_clones_inherit_from_children(
                solver_result.solution, self.tree, self.cp_table
            )
        self.get_model_metrics(solver_result)

    def _events_from_diploid(self, pred_cn_a: int, pred_cn_b: int) -> int:
        """Complexity of a solution in which every clone has the same copy numbers."""
        if not self.config["model_config"].get("minimise_events_to_diploid", True):
            return 0
        return abs(pred_cn_a - 1) + abs(pred_cn_b - 1)

    def _shortcut_result(
        self, pred_cn_a: int, pred_cn_b: int, ci_score: int, d_score: float, shortcut: str
    ) -> SolverResult:
        """SolverResult for a solution found without the ILP (every clone gets the same state)."""
        complexity = self._events_from_diploid(pred_cn_a, pred_cn_b)
        solution = pd.DataFrame({"clone": self.cp_table.index.to_list()})
        solution["pred_CN_A"] = pred_cn_a
        solution["pred_CN_B"] = pred_cn_b
        solution["complexity"] = complexity
        solution["CI_score"] = ci_score
        solution["D_score"] = d_score
        solution["variability_penalty_count"] = 0
        solution["state_change_count"] = 0
        solution["event_count"] = complexity
        solution["allowed_complexity"] = 0
        solution["gurobi_time_CI"] = -1
        solution["gurobi_gap_CI"] = -1
        solution["gurobi_time_D"] = -1
        solution["gurobi_gap_D"] = -1
        solution = solution.sort_values("clone").reset_index(drop=True)
        return SolverResult(
            solution=solution,
            runtime=0.0,
            backend_name=shortcut,
            gap_status={
                "max_gap": 0,
                "gap_reason": shortcut,
                "runtime": 0.0,
                "status": 0,
            },
            objective_values={"D_score": d_score, "CI_score": ci_score},
            solver_specific={"shortcut": shortcut},
        )

    def _build_clonal_shortcut_result(self) -> Optional[SolverResult]:
        if "is_clonal" not in self.input_table.columns:
            return None
        is_clonal = coerce_bool_like_series(self.input_table["is_clonal"])
        if is_clonal.isna().any() or not bool(is_clonal.all()):
            return None
        rounded_cns = self.input_table[["cpnA", "cpnB"]].round().astype(int)
        unique_states = rounded_cns.drop_duplicates().reset_index(drop=True)
        if len(unique_states) != 1:
            return None
        return self._shortcut_result(
            int(unique_states.at[0, "cpnA"]),
            int(unique_states.at[0, "cpnB"]),
            ci_score=0,
            d_score=0.0,
            shortcut="clonal_shortcut",
        )

    def _build_single_clone_result(self) -> SolverResult:
        (pred_cn_a, pred_cn_b), scores, flags, details = solve_single_clone(
            self.input_table, self.ci_table, self.segment, self.config["model_config"]
        )
        self.edge_case_flags += flags
        self.edge_case_details = details
        return self._shortcut_result(
            pred_cn_a,
            pred_cn_b,
            ci_score=scores["CI_score"],
            d_score=scores["D_score"],
            shortcut="single_clone_direct",
        )

    def _use_shortcut_result(self, result: SolverResult) -> None:
        shortcut = result.backend_name
        self.used_clonal_shortcut = True
        self.solver_path = shortcut
        self.maximum_complexity = 0
        self.get_model_metrics(result)
        self.solutions_combined = result.solution.copy()
        self.elbow_search_df = self.solutions_combined[
            ["complexity", "D_score", "CI_score", "allowed_complexity"]
        ].drop_duplicates(subset="allowed_complexity", keep="first")
        self.elbow_search_df_strictly_decreasing = self.elbow_search_df.copy()
        self.elbow = {
            "s_min": 0,
            "s_raw": 0,
            "s_strictly_decreasing": 0,
            "raw_code": shortcut,
            "dec_code": shortcut,
        }

    def _run_lowest_feasible_model(self) -> int:
        """Run the model at complexity 0, or at the lowest complexity that is feasible.

        Restricting heterogeneity fixes some clones' copy numbers, which can take more
        events from diploid than a low complexity cap allows. Those caps are skipped
        without writing infeasibility reports. Returns the complexity that was run.
        """
        if not self.config["model_config"].get("restrict_heterogeneity_flag", False):
            self.run_model(allowed_complexity=0)
            return 0
        for c in range(0, self.maximum_complexity):
            try:
                self.run_model(allowed_complexity=c, report_infeasibility=False)
                return c
            except SegmentInfeasibleError:
                continue
        # infeasible at every complexity: run once more so the infeasibility is reported
        self.run_model(allowed_complexity=self.maximum_complexity)
        return self.maximum_complexity

    def stop_conditions_check(self, oft):
        optimization_time = self.metrics["run_time"][-1]
        time_limit = self.metrics["time_limits"][-1] if self.metrics["time_limits"] else None
        slow_iteration = (
            optimization_time is not None
            and time_limit is not None
            and optimization_time >= time_limit
        )
        no_improvement_in_D_score = (
            sum(abs(np.diff(self.metrics["D_scores"])) < oft) > 3
        )
        at_least_3_complexities = len(set(self.metrics["complexity"])) >= 3
        return no_improvement_in_D_score & at_least_3_complexities & slow_iteration

    def run_iterations(self):
        if self.tumour_category in (SINGLE_CLONE, SINGLE_CLONE_AND_SINGLE_SAMPLE):
            print(f"Skipping solver for {self.segment}: the tumour has a single clone.")
            self._use_shortcut_result(self._build_single_clone_result())
            return
        shortcut_result = self._build_clonal_shortcut_result()
        if shortcut_result is not None:
            print(
                f"Skipping solver for {self.segment}: segment is clonal and integer "
                "copy numbers are identical across samples."
            )
            self._use_shortcut_result(shortcut_result)
            return
        # use heuristics to determine maximum complexity:
        self.maximum_complexity = max(
            20,
            len(self.input_table["sample"].unique())
            * math.ceil(self.input_table[["cpnA", "cpnB"]].max().max()),
        )
        
        objective_function_threshold = 0.1  # iterations will stop if D score does not improve by more than this value in 3 consecutive iterations
        # If user provided both a debug_solution_file and an explicit complexity,
        # run a single iteration with that complexity and emit raw Gurobi logs
        override_complexity = None
        try:
            override_complexity = self.config["model_config"].get("complexity")
        except Exception:
            override_complexity = None
        if self.config["model_config"].get("debug_solution_file") and override_complexity is not None:
            print(f"Running single iteration with user-specified complexity: {override_complexity}")
            # ensure gurobi logs are written to a file in the output directory
            out_dir = self.config["preprocessing_config"].get("output_directory", ".")
            out_dir = os.path.expanduser(out_dir)
            os.makedirs(out_dir, exist_ok=True)
            gurobi_log_path = os.path.join(out_dir, f"gurobi_log_{self.tumour_id}_{self.segment}.txt")
            # temporarily set model_config gurobi_logs so Model will write logs
            self.config["model_config"]["gurobi_logs"] = gurobi_log_path
            # run diploid and then the specified complexity only
            print(f"**Single iteration with complexity: {override_complexity}")
            self.run_model(allowed_complexity=override_complexity)
        else:
            # run diploid model (or the lowest complexity the restrictions allow):
            first_complexity = self._run_lowest_feasible_model()
            self.elbow_increase_report = pd.DataFrame()
            # don't iterate if solution is likely to be diploid:
            if self.metrics["D_scores"][0] > objective_function_threshold:
                complexity_range = range(first_complexity + 1, self.maximum_complexity)
                for c in complexity_range:
                    print(f"**Iterating with complexity: {c}")
                    self.run_model(allowed_complexity=c)
                    stop_conditions = self.stop_conditions_check(
                        objective_function_threshold
                    )
                    if stop_conditions:
                        # check if elbow can be found
                        self.find_elbow()
                        elbow_findable = (
                            self.elbow["s_min"] < 1000
                            and self.elbow["raw_code"] != "elbow_failure"
                        )
                        if elbow_findable:
                            print(
                                f"** Stopping iterations at complexity {c} due to lack of improvement in D score"
                            )
                            break
        self.solutions_combined = pd.concat(self.metrics["solutions"])

    def find_elbow(self):
        assert self.metrics is not None, "Metrics not found, run iterations first"
        solutions_combined = pd.concat(self.metrics["solutions"])
        self.elbow_search_df = (
            solutions_combined[
                ["complexity", "D_score", "CI_score", "allowed_complexity"]
            ]
            .drop_duplicates(subset="allowed_complexity", keep="first")
            .sort_values("allowed_complexity")
            .reset_index(drop=True)
        )
        # check if objective function increases at any point. If yes, it signifies and issue, as error values should never increase.
        # this can happen if confidence intervals are very narrow and if both objectives are used for optimization. The model
        # always trys to minimize CI score, and in some rare cases this can lead to an increase in D score. In such a scenario, it is
        # recommended to use min_ci parameter to enforce a minimum confidence interval width.
        d_score_increases = (np.diff(self.elbow_search_df.D_score) > 0).any()
        if d_score_increases:
            elbow_increase_report_df = self.elbow_search_df.copy()
            elbow_increase_report_df["issue"] = ['']+['D_score_increase' if x > 0 else '' for x in np.diff(self.elbow_search_df.D_score)]
            elbow_increase_report_df['tumour_id'] = self.tumour_id
            elbow_increase_report_df['segment'] = self.segment
            self.elbow_increase_report = elbow_increase_report_df

        self.elbow_search_df_strictly_decreasing = ensure_elbow_strictly_decreasing(
            self.elbow_search_df.copy()
        )
        s_raw, raw_code = find_s_values(
            self.elbow_search_df,
            self.maximum_complexity,
            "allowed_complexity",
            "D_score",
        )
        s_strictly_decreasing, dec_code = find_s_values(
            self.elbow_search_df_strictly_decreasing,
            self.maximum_complexity,
            "allowed_complexity",
            "D_score",
        )
        # in single sample tumours the elbow curve can be a straight line with no knee;
        # fall back to the most complex solution instead of failing:
        elbow_failure = (s_raw is None or s_strictly_decreasing is None) and self.tumour_category == SINGLE_SAMPLE
        if elbow_failure:
            most_complex = int(self.elbow_search_df["allowed_complexity"].max())
            print(f"Elbow not found for {self.segment}: selecting the most complex solution ({most_complex})")
            s_raw = s_strictly_decreasing = most_complex
            raw_code = dec_code = "elbow_failure"
        assert s_raw is not None, "S_raw not found"
        assert s_strictly_decreasing is not None, "S_strictly_decreasing not found"
        s_min = min(s_raw, s_strictly_decreasing)
        s_values = {
            "s_min": s_min,
            "s_raw": s_raw,
            "s_strictly_decreasing": s_strictly_decreasing,
            "raw_code": raw_code,
            "dec_code": dec_code,
        }
        self.elbow = s_values
        self.optimal_solution_index = self.elbow[self.s_type]
        # required for certain simulated scenarios:
        if (
            self.config["model_config"]["d_zero"]
            & (self.elbow_search_df.D_score == 0).any()
        ):
            self.optimal_solution_index = int(
                self.elbow_search_df.query("D_score == 0")["allowed_complexity"].iloc[0]
            )

    def find_optimal_solution(self):
        # check if diploid solution was found:
        assert self.metrics is not None, "Metrics not found, run iterations first"
        diploid_solution_found = len(set(self.metrics["D_scores"])) == 1
        if diploid_solution_found:
            # the lowest complexity that was run (0 unless restrictions made it infeasible)
            self.optimal_solution_index = int(
                pd.concat(self.metrics["solutions"])["allowed_complexity"].min()
            )
        else:
            self.find_elbow()
            # add metadata to the elbow search dataframe:
            self.elbow_search_df_strictly_decreasing["segment"] = self.segment
            self.elbow_search_df_strictly_decreasing["tumour_id"] = self.tumour_id
            self.elbow_search_df_strictly_decreasing["opt_sol_indx"] = (
                self.optimal_solution_index
            )

    def _get_columns_to_keep(self, base_columns, available_columns):
        columns_to_keep = list(base_columns)
        if self.extra_columns:
            expanded_extra = []
            for col in self.extra_columns:
                expanded_extra.append(col)
                if col == "gurobi_time":
                    expanded_extra.extend(["gurobi_time_D", "gurobi_time_CI"])
                if col == "gurobi_gap":
                    expanded_extra.extend(["gurobi_gap_D", "gurobi_gap_CI"])
            columns_to_keep.extend(expanded_extra)
            
        # remove duplicates preserving order
        columns_to_keep = list(dict.fromkeys(columns_to_keep))
        
        # Only keep columns that actually exist in the dataframe
        columns_to_keep = [c for c in columns_to_keep if c in available_columns]
        return columns_to_keep

    def get_solution(self, s=None):
        if s is None:
            s = self.optimal_solution_index
        self.optimal_solution = self.solutions_combined.query(
            f"allowed_complexity == {s}"
        ).copy()
        self.optimal_solution["tumour_id"] = self.tumour_id
        self.optimal_solution["segment"] = self.segment
        if not self.debug:
            columns_to_drop = [
                "allowed_complexity",
                "variability_penalty_count",
                "state_change_count",
                "event_count",
            ]
            if self.extra_columns:
                columns_to_drop = [c for c in columns_to_drop if c not in self.extra_columns]

            self.optimal_solution.drop(
                columns=columns_to_drop,
                inplace=True,
                errors="ignore",
            )
            
            base_columns = ["tumour_id", "segment", "clone", "pred_CN_A", "pred_CN_B", "complexity"]
            columns_to_keep = self._get_columns_to_keep(base_columns, self.optimal_solution.columns)
            
            self.optimal_solution = self.optimal_solution[columns_to_keep]

    def check_ambiguity(self):
        """Flag the segment if another assignment fits a single-sample tumour equally well.

        The ILP is re-solved at the selected complexity with the CI and D scores capped
        at the selected solution's values and the copy numbers of the clones present
        in the sample required to differ. Clones absent from the sample are left out:
        no data constrains them.
        """
        if self.tumour_category != SINGLE_SAMPLE or self.used_clonal_shortcut:
            return
        chosen = self.solutions_combined.query(
            f"allowed_complexity == {self.optimal_solution_index}"
        )
        chosen = chosen[chosen.clone != "diploid"].set_index("clone")
        present = self.cp_table.index[self.cp_table.sum(axis=1) > 0]
        reference = {
            (allele, clone): int(chosen.loc[clone, f"pred_CN_{allele}"])
            for clone in present
            for allele in ("A", "B")
        }
        objectives = str(self.config["model_config"].get("objectives", "DCI")).upper()
        max_ci = float(chosen["CI_score"].iloc[0]) if "CI" in objectives else None
        # D scores are rounded to 3 decimals in the solution frame
        max_d = float(chosen["D_score"].iloc[0]) + 1e-3 if "D" in objectives else None
        backend_config = {
            **self.config["model_config"],
            "allowed_tree_complexity": self.optimal_solution_index,
            "output_directory": self.config["preprocessing_config"].get("output_directory", "."),
        }
        try:
            backend = create_solver_backend(self.solver_name, self.model_inputs, backend_config)
            alternative = backend.find_alternative(
                self.optimal_solution_index, reference, max_ci, max_d
            )
        except Exception as e:
            self.edge_case_flags.append(AMBIGUITY_CHECK_FAILED)
            self.edge_case_details = f"ambiguity check failed: {e}"
            return
        if alternative is None:
            return
        alternative = alternative.set_index("clone")
        differences = [
            f"{clone}: ({chosen.loc[clone, 'pred_CN_A']},{chosen.loc[clone, 'pred_CN_B']}) or "
            f"({alternative.loc[clone, 'pred_CN_A']},{alternative.loc[clone, 'pred_CN_B']})"
            for clone in present
            if (chosen.loc[clone, "pred_CN_A"], chosen.loc[clone, "pred_CN_B"])
            != (alternative.loc[clone, "pred_CN_A"], alternative.loc[clone, "pred_CN_B"])
        ]
        self.edge_case_flags.append(AMBIGUOUS)
        self.edge_case_details = "equally good alternative: " + "; ".join(differences)

    def _save_edge_case_report(self, output_dir: str) -> None:
        if self.tumour_category == STANDARD:
            return
        report = make_segment_report(
            tumour_id=self.tumour_id,
            segment=self.segment,
            category=self.tumour_category,
            n_clones=len(self.cp_table.index),
            n_samples=len(self.cp_table.columns),
            solver_path=self.solver_path,
            restrict_heterogeneity=bool(
                self.config["model_config"].get("restrict_heterogeneity_flag", False)
            ),
            flags=self.edge_case_flags,
            details=self.edge_case_details,
        )
        fmt = getattr(self, "output_format", "parquet").lower().lstrip(".")
        save_dataframe(
            report,
            output_dir,
            segment_report_name(self.tumour_id, self.segment, fmt),
            output_format=fmt,
        )

    def _get_all_solutions_subdir_name(self, all_dir: str) -> str:
        all_dir_seg = os.path.join(all_dir, self.segment)
        os.makedirs(all_dir_seg, exist_ok=True)
        return all_dir_seg

    def _save_all_solutions(self, all_dir: str, all_solutions: pd.DataFrame) -> str:
        all_dir_seg = self._get_all_solutions_subdir_name(all_dir)
        all_solutions_output_name = os.path.basename(self.create_output_path()).replace(
            "optimal", "all"
        )

        all_solutions_output_path = os.path.join(all_dir_seg, all_solutions_output_name)
        fmt = getattr(self, "output_format", "parquet")
        save_dataframe(all_solutions, all_dir_seg, all_solutions_output_name, output_format=fmt)
        return all_solutions_output_path

    def _save_elbow_table(self, all_dir: str) -> typing.Tuple[typing.Optional[str], typing.Optional[pd.DataFrame]]:
        try:
            if hasattr(self, "elbow_search_df") and self.elbow_search_df is not None:
                elbow_df = self.elbow_search_df.copy()
            else:
                elbow_df = (
                    self.solutions_combined[["complexity", "D_score", "CI_score", "allowed_complexity"]]
                    .drop_duplicates(subset="allowed_complexity", keep="first")
                    .sort_values("allowed_complexity")
                    .reset_index(drop=True)
                )
            elbow_meta = {
                "knee_s_min": self.elbow.get("s_min") if isinstance(self.elbow, dict) else None,
                "knee_s_raw": self.elbow.get("s_raw") if isinstance(self.elbow, dict) else None,
                "knee_s_strictly_decreasing": self.elbow.get("s_strictly_decreasing") if isinstance(self.elbow, dict) else None,
                "knee_raw_code": self.elbow.get("raw_code") if isinstance(self.elbow, dict) else None,
                "knee_dec_code": self.elbow.get("dec_code") if isinstance(self.elbow, dict) else None,
                "selected_by_s_type": self.optimal_solution_index,
                "selected_s_type": self.s_type,
            }
            for k, v in elbow_meta.items():
                elbow_df[k] = v
            all_dir_seg = self._get_all_solutions_subdir_name(all_dir)
            fmt = getattr(self, "output_format", "parquet").lower().lstrip(".")
            elbow_output_filename = f"{self.tumour_id}_{self.segment}_elbow_table.{fmt}"
            elbow_output_path = os.path.join(all_dir_seg, elbow_output_filename)
            # add segment and tumour_id
            elbow_df["segment"] = self.segment if hasattr(self, "segment") else None
            elbow_df["tumour_id"] = self.tumour_id if hasattr(self, "tumour_id") else None
            elbow_df["optimal_complexity"] = self.optimal_solution['complexity'].unique()[0]
            save_dataframe(elbow_df, all_dir_seg, elbow_output_filename, output_format=fmt)
            return elbow_output_path, elbow_df
        except Exception:
            return None, None

    def _plot_elbow(self, all_dir: str, elbow_df: pd.DataFrame) -> typing.Optional[str]:
        """Create a simple elbow plot (D_score vs allowed_complexity) and mark the selected complexity. Returns path or None on failure."""
        try:
            if plt is None:
                return None
            fig, ax = plt.subplots()
            x = elbow_df["allowed_complexity"]
            y = elbow_df["D_score"]
            ax.plot(x, y, marker="o", linestyle="-", label="D_score")
            selected = self.optimal_solution_index
            if selected is not None:
                ax.axvline(selected, color="orange", linestyle="--", label=f"selected: {selected}")
                sel_y = None
                try:
                    sel_y = float(elbow_df.query(f"allowed_complexity == {selected}").D_score.iloc[0])
                except Exception:
                    sel_y = None
                if sel_y is not None:
                    ax.plot([selected], [sel_y], marker="*", color="red", markersize=12)
            ax.set_xlabel("allowed_complexity")
            ax.set_ylabel("D_score")
            ax.set_title(f"Elbow curve: {self.tumour_id} {self.segment}")
            ax.legend()
            all_dir_seg = self._get_all_solutions_subdir_name(all_dir)
            plot_path = os.path.join(all_dir_seg, f"{self.tumour_id}_{self.segment}_elbow_plot.png")
            fig.savefig(plot_path, bbox_inches="tight")
            plt.close(fig)
            return plot_path
        except Exception:
            return None

    def _run_and_save_unconstrained(self, all_dir: str) -> typing.Optional[str]:
        """Run the model once without enforcing the tree complexity constraint and save the resulting solution.
        The result is saved as `all_max_<tumour_id>_<segment>.csv` inside all_dir. This run is NOT used for elbow search,
        as it can interfere with how kneed finds the elbow.
        """
        try:
            model_kwargs = {**self.config.get("model_config", {}), "enforce_tree_complexity": False}
            # ensure we do not change allowed_tree_complexity here and let model compute unconstrained optimum
            backend = create_solver_backend(
                model_kwargs.get("solver", "pyomo"),
                self.model_inputs,
                model_kwargs,
            )
            unconstrained_complexity = model_kwargs.get(
                "allowed_tree_complexity",
                self.maximum_complexity or self.config["model_config"].get("allowed_tree_complexity", 0),
            )
            solver_result = backend.solve(allowed_complexity=unconstrained_complexity)
            sol = solver_result.solution
            sol = sol[sol.clone != "diploid"]
            # add tumour and segment columns if needed
            sol["tumour_id"] = self.tumour_id
            sol["segment"] = self.segment
            fmt = getattr(self, "output_format", "parquet").lower().lstrip(".")
            out_name = f"all_max_{self.tumour_id}_{self.segment}.{fmt}"
            all_dir_seg = self._get_all_solutions_subdir_name(all_dir)
            out_path = os.path.join(all_dir_seg, out_name)
            save_dataframe(sol, all_dir_seg, out_name, output_format=self.output_format)
            return out_path
        except Exception:
            return None

    def get_monoclonal_samples_report(self):
        clones_per_sample = (self.cp_table > 0).sum()
        monoclonal_samples = clones_per_sample[clones_per_sample == 1].index.values
        monoclonal_samples_copynumbers = self.input_table[self.input_table['sample'].isin(monoclonal_samples)].copy()
        monoclonal_samples_copynumbers['distance_to_integer_A'] = abs(monoclonal_samples_copynumbers['cpnA'] - monoclonal_samples_copynumbers['cpnA'].round())
        monoclonal_samples_copynumbers['distance_to_integer_B'] = abs(monoclonal_samples_copynumbers['cpnB'] - monoclonal_samples_copynumbers['cpnB'].round())
        self.monoclonal_samples_report = monoclonal_samples_copynumbers

    def _save_run_summary_report(self, output_dir: str):
        """Create a run_summary report showing segments with non-zero gap and reasons."""
        try:
            if not hasattr(self, "metrics") or "gap_status" not in self.metrics:
                return
            
            # Get gap status for the optimal solution
            optimal_idx = self.optimal_solution_index
            if optimal_idx is None:
                return
            run_complexities = [int(sol["allowed_complexity"].iloc[0]) for sol in self.metrics["solutions"]]
            if optimal_idx not in run_complexities:
                return
            gap_info = self.metrics["gap_status"][run_complexities.index(optimal_idx)]
            max_gap = gap_info.get("max_gap", -1)
            
            # Only report if gap > 0 (non-optimal solution)
            if max_gap <= 0:
                return
            
            # Create report entry
            report_entry = {
                "tumour_id": self.tumour_id,
                "segment": self.segment,
                "max_gap": max_gap,
                "gap_reason": gap_info.get("gap_reason", "unknown"),
                "runtime": gap_info.get("runtime", -1),
                "optimal_complexity": optimal_idx,
                "strict_gap_enabled": getattr(self, "strict_gap", True),
            }
            
            # Save to segment-specific run_gap_summary file
            fmt = getattr(self, "output_format", "parquet").lower().lstrip(".")
            run_summary_name = f"{self.tumour_id}_{self.segment}_run_gap_summary.{fmt}"
            run_summary_path = os.path.join(output_dir, run_summary_name)
            df = pd.DataFrame([report_entry])
            save_dataframe(df, output_dir, run_summary_name, output_format=self.output_format)
        except Exception:
            # Best effort; do not fail if run_summary cannot be created
            pass

    def set_directories(self):
        is_running_in_nextflow = self.config["preprocessing_config"].get('mode') == 'segment'
        self.tumour_dir = self.config["preprocessing_config"].get('input_tumour_directory')
        if is_running_in_nextflow:
            self.segments_dir = self.input_data_directory
        else:
            self.segments_dir = f"{self.tumour_dir}/segments"

    def output_exists(self):
        """
        Checks if the output file already exists.
        """
        output_path = self.create_output_path()
        return os.path.exists(output_path)

    def create_output_path(self):
        """
        Creates the output path for the solution based on file name and options.
        """
        output_format = getattr(self, "output_format", "parquet")
        fmt = output_format.lower().lstrip(".")
        base = self.input_file_name.split("ALPACA_input_table_")[1]
        base = base.replace(".parquet", "").replace(".csv", "")
        output_name = f"optimal_{base}.{fmt}"
        output_dir = self.config["preprocessing_config"]["output_directory"]
        output_path = os.path.join(output_dir, output_name)
        return output_path

    def save_output(self):
        logger = self.logger
        end_time = time.time()
        output_path = self.create_output_path()
        output_dir = os.path.dirname(output_path)
        os.makedirs(output_dir, exist_ok=True)
        # discard diploid clone:
        assert self.optimal_solution is not None
        self.optimal_solution = self.optimal_solution[
            self.optimal_solution.clone != "diploid"
        ]
        if self.output_all_solutions:
            try:
                # reconstruct optimal complexity from the combined solutions
                opt_rows = self.solutions_combined.query(
                    f"allowed_complexity == {self.optimal_solution_index}"
                ).copy()
                elbow = opt_rows["complexity"].iloc[0]
            except Exception:
                # fallback: if complexity still present on optimal_solution use it,
                # otherwise set elbow to None and skip elbow_offset computation
                elbow = (
                    self.optimal_solution["complexity"].iloc[0]
                    if "complexity" in self.optimal_solution.columns
                    else None
                )

            base_columns = ["clone", "pred_CN_A", "pred_CN_B", "complexity", "allowed_complexity"]
            columns_to_keep = self._get_columns_to_keep(base_columns, self.solutions_combined.columns)

            all_solutions = self.solutions_combined[columns_to_keep].copy()
            # remove diploid clone from the all_solutions table as well
            all_solutions = all_solutions[all_solutions.clone != "diploid"]
            all_solutions["tumour_id"] = self.tumour_id
            all_solutions["segment"] = self.segment
            if elbow is not None:
                complexities = all_solutions.complexity.unique()
                rescaled = rescale_elbow_points(complexities, elbow)
                all_solutions["elbow_offset"] = all_solutions.complexity.map(rescaled)
            else:
                all_solutions["elbow_offset"] = 0

            # write all solutions and elbow artifacts into a dedicated subdirectory for clarity
            all_dir = os.path.join(output_dir, "all_solutions")
            os.makedirs(all_dir, exist_ok=True)
            try:
                self._save_all_solutions(all_dir, all_solutions)
                _, elbow_df = self._save_elbow_table(all_dir)
                if self.write_elbow_plots and elbow_df is not None:
                    self._plot_elbow(all_dir, elbow_df)
                # run an unconstrained (no tree_complexity_constr) max-complexity solution
                # and save it separately so it does not interfere with elbow search
                if not self.used_clonal_shortcut:
                    try:
                        self._run_and_save_unconstrained(all_dir)
                    except Exception:
                        pass
            except Exception:
                # best-effort; do not fail the run if elbow saving fails
                pass
        if self.output_model_selection_table:
            output_model_selection_table = self.elbow_search_df_strictly_decreasing
            fmt = getattr(self, "output_format", "parquet").lower().lstrip(".")
            model_sel_name = f"{self.tumour_id}_{self.segment}_model_selection_table.{fmt}"
            save_dataframe(
                output_model_selection_table,
                output_dir,
                model_sel_name,
                output_format=self.output_format,
            )
        # add runtime metadata when debugging
        total_run_time = round(end_time - self.start_time)
        if self.debug:
            self.optimal_solution["run_time_seconds"] = total_run_time

        # Before saving the final optimal solution, drop internal columns including
        # 'complexity' to keep the output minimal for non-debug runs.
        if not self.debug:
            columns_to_drop = [
                "allowed_complexity",
                "complexity",
                "variability_penalty_count",
                "state_change_count",
                "event_count",
                "gurobi_time",
                "gurobi_gap",
            ]
            # do not drop columns requested by user
            if self.extra_columns:
                columns_to_drop = [c for c in columns_to_drop if c not in self.extra_columns]

            self.optimal_solution.drop(
                columns=columns_to_drop,
                inplace=True,
                errors="ignore",
            )

        # always write the optimal solution file
        save_dataframe(
            self.optimal_solution,
            output_dir,
            os.path.basename(output_path),
            output_format=self.output_format,
        )
        # save reports if available:
        # ci_report:
        if self.ci_report:
            ci_report_path = os.path.join(
                output_dir,
                f"{self.tumour_id}_{self.segment}_ci_report.json",
            )
            with open(ci_report_path, "w") as f:
                json.dump(self.ci_report, f, indent=4)
        if not self.monoclonal_samples_report.empty:
            fmt = getattr(self, "output_format", "parquet").lower().lstrip(".")
            mono_report_name = f"{self.tumour_id}_{self.segment}_monoclonal_samples_report.{fmt}"
            save_dataframe(
                self.monoclonal_samples_report,
                output_dir,
                mono_report_name,
                output_format=self.output_format,
            )
        if not self.elbow_increase_report.empty:
            fmt = getattr(self, "output_format", "parquet").lower().lstrip(".")
            elbow_report_name = f"{self.tumour_id}_{self.segment}_elbow_increase_report.{fmt}"
            save_dataframe(
                self.elbow_increase_report,
                output_dir,
                elbow_report_name,
                output_format=self.output_format,
            )
        
        # Save run gap summary report with gap status
        self._save_run_summary_report(output_dir)
        self._save_edge_case_report(output_dir)
        
        if os.path.exists(output_path):
            logger.info("Segment output created")
        else:
            logger.error(f"Output not saved to {output_path}")

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
from alpaca.utils import read_tree_json
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


def split_input_file_name(input_file_name: str):
    stripped_name = input_file_name.split("ALPACA_input_table_")[1]
    assert stripped_name.count("_") == 3, "Input name has to many underscores"
    stripped_name = stripped_name.replace(".csv", "")
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
        self.output_model_selection_table: bool = False
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
        self.tumour_id, self.segment = split_input_file_name(self.input_file_name)
        # define tumour input directory depending on the run environment:
        self.set_directories()
        # load fractional copy numbers:
        self.input_table = pd.read_csv(
            f"{self.segments_dir}/{input_file_name}"
        ).sort_values("sample")
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
        self.solver_name = self.config["model_config"].get("solver", "gurobi")
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

    def run_model(self, allowed_complexity):
        backend_config = {**self.config["model_config"], "allowed_tree_complexity": allowed_complexity}
        solver_name = backend_config.get("solver", "gurobi")
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

        pred_cn_a = int(unique_states.at[0, "cpnA"])
        pred_cn_b = int(unique_states.at[0, "cpnB"])
        solution = pd.DataFrame({"clone": self.cp_table.index.to_list()})
        solution["pred_CN_A"] = pred_cn_a
        solution["pred_CN_B"] = pred_cn_b
        solution["complexity"] = 0
        solution["CI_score"] = 0
        solution["D_score"] = 0.0
        solution["variability_penalty_count"] = 0
        solution["state_change_count"] = 0
        solution["event_count"] = 0
        solution["allowed_complexity"] = 0
        solution["gurobi_time_CI"] = -1
        solution["gurobi_gap_CI"] = -1
        solution["gurobi_time_D"] = -1
        solution["gurobi_gap_D"] = -1
        solution = solution.sort_values("clone").reset_index(drop=True)
        return SolverResult(
            solution=solution,
            runtime=0.0,
            backend_name="clonal_shortcut",
            gap_status={
                "max_gap": 0,
                "gap_reason": "clonal_shortcut",
                "runtime": 0.0,
                "status": 0,
            },
            objective_values={"D_score": 0.0, "CI_score": 0},
            solver_specific={"shortcut": "is_clonal_identical_integer_cns"},
        )

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
        shortcut_result = self._build_clonal_shortcut_result()
        if shortcut_result is not None:
            print(
                f"Skipping solver for {self.segment}: segment is clonal and integer "
                "copy numbers are identical across samples."
            )
            self.used_clonal_shortcut = True
            self.maximum_complexity = 0
            self.get_model_metrics(shortcut_result)
            self.solutions_combined = shortcut_result.solution.copy()
            self.elbow_search_df = self.solutions_combined[
                ["complexity", "D_score", "CI_score", "allowed_complexity"]
            ].drop_duplicates(subset="allowed_complexity", keep="first")
            self.elbow_search_df_strictly_decreasing = self.elbow_search_df.copy()
            self.elbow = {
                "s_min": 0,
                "s_raw": 0,
                "s_strictly_decreasing": 0,
                "raw_code": "clonal_shortcut",
                "dec_code": "clonal_shortcut",
            }
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
            # run diploid model:
            self.run_model(allowed_complexity=0)
            self.elbow_increase_report = pd.DataFrame() 
            # don't iterate if solution is likely to be diploid:
            if self.metrics["D_scores"][0] > objective_function_threshold:
                complexity_range = range(1, self.maximum_complexity)
                for c in complexity_range:
                    print(f"**Iterating with complexity: {c}")
                    self.run_model(allowed_complexity=c)
                    stop_conditions = self.stop_conditions_check(
                        objective_function_threshold
                    )
                    if stop_conditions:
                        # check if elbow can be found
                        self.find_elbow()
                        elbow_findable = self.elbow["s_min"] < 1000
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
            self.optimal_solution_index = self.elbow_search_df.query(
                "D_score == 0"
            ).index[0]

    def find_optimal_solution(self):
        # check if diploid solution was found:
        assert self.metrics is not None, "Metrics not found, run iterations first"
        diploid_solution_found = len(set(self.metrics["D_scores"])) == 1
        if diploid_solution_found:
            self.optimal_solution_index = 0
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
        all_solutions.to_csv(all_solutions_output_path, index=False)
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
            elbow_output_path = os.path.join(all_dir_seg, f"{self.tumour_id}_{self.segment}_elbow_table.csv")
            elbow_df.to_csv(elbow_output_path, index=False)
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
                model_kwargs.get("solver", "gurobi"),
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
            out_name = f"all_max_{self.tumour_id}_{self.segment}.csv"
            all_dir_seg = self._get_all_solutions_subdir_name(all_dir)
            out_path = os.path.join(all_dir_seg, out_name)
            sol.to_csv(out_path, index=False)
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
            if optimal_idx is None or optimal_idx >= len(self.metrics["gap_status"]):
                return
            
            gap_info = self.metrics["gap_status"][optimal_idx]
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
            run_summary_path = os.path.join(
                output_dir,
                f"{self.tumour_id}_{self.segment}_run_gap_summary.csv",
            )
            df = pd.DataFrame([report_entry])
            df.to_csv(run_summary_path, index=False)
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
        output_name = "optimal_" + self.input_file_name.split("ALPACA_input_table_")[1]
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
                elbow_path, elbow_df = self._save_elbow_table(all_dir)
                if elbow_df is not None:
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
            output_model_selection_table.to_csv(
                f"{output_dir}/{self.tumour_id}_{self.segment}_model_selection_table.csv",
                index=False,
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
        self.optimal_solution.to_csv(output_path, index=False)
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
            mono_report_path = os.path.join(
                output_dir,
                f"{self.tumour_id}_{self.segment}_monoclonal_samples_report.csv",
            )
            self.monoclonal_samples_report.to_csv(mono_report_path, index=False)            
        if not self.elbow_increase_report.empty:
            elbow_report_path = os.path.join(
                output_dir,
                f"{self.tumour_id}_{self.segment}_elbow_increase_report.csv",
            )
            self.elbow_increase_report.to_csv(elbow_report_path, index=False)
        
        # Save run gap summary report with gap status
        self._save_run_summary_report(output_dir)
        
        if os.path.exists(output_path):
            logger.info("Segment output created")
        else:
            logger.error(f"Output not saved to {output_path}")

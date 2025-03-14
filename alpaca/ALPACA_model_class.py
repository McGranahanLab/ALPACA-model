import gurobipy as gp
import pandas as pd
import os
from gurobipy import GRB
from utils import find_path_edges, get_tree_edges, flat_list, get_length_from_name


class Model:
    """
    Main class for ALPACA model
    """

    def __init__(
        self,
        segment,
        ci_table,
        fractional_copy_number_table,
        tree,
        clone_proportions,
        **kwargs,
    ):
        # default parameters:
        self.homozygous_deletion_threshold = 1
        self.homo_del_size_limit = 5 * 10**7
        self.limit_homozygous_deletions_threshold_flag = bool(
            self.homozygous_deletion_threshold
        )
        self.add_event_count_constraints_flag = True
        self.prevent_increase_from_zero_flag = True
        self.add_state_change_count_constraints_flag = False
        self.add_path_variability_penalty_constraints_flag = False
        self.add_allow_only_one_non_directional_event_flag = True
        self.variability_penalty = 0
        self.allowed_tree_complexity = 1000
        self.minimise_events_to_diploid = True
        self.exclusive_amp_del = True
        self.two_objectives = True
        self.restrict_heterogeneity_flag = False
        self.restrict_to_clonal_only_flag = False
        self.time_limit = 60
        self.cpus = 2
        self.BestObjStop = None
        self.license = "local"

        # override defaults:
        self.__dict__.update(kwargs)
        # ::::: ILP variables:

        # copy number variables:
        self.X = (
            {}
        )  # Dictionary containing predicted integer copy number values per clone
        self.Y = (
            {}
        )  # Dictionary containing observed fractional copy number values for each sample
        self.d = (
            {}
        )  # Dictionary containing distance between predicted fractional and observed fractional copy number values
        self.Yhat = {}  # Dictionary containing predicted fractional copy number values

        # confidence intervals variables:
        self.CI_upper = {}
        self.CI_lower = {}
        self.yhat_above_upper_CI = {}
        self.yhat_below_lower_CI = {}
        self.CI_overlap = {}

        # events variables:
        self.CN_diff_edges_amp = {}
        self.CN_diff_edges_del = {}
        self.n = {}
        self.total_events = None
        self.cpn_change_up = {}
        self.cpn_change_down = {}
        self.total_edge_changes = None
        self.path_variability_penalty = {}
        self.total_path_variability_penalty = None
        self.amps_count_on_path = {}
        self.more_than_1_amp_change = {}
        self.dels_count_on_path = {}
        self.more_than_1_del_change = {}

        # :::::  inputs:
        self.fractional_copy_number_table = fractional_copy_number_table
        self.ci_table = ci_table
        self.tree = tree
        self.clone_proportions = clone_proportions
        self.segment = segment
        self.tumour_id = fractional_copy_number_table["tumour_id"].unique()[0]
        self.sample_names = fractional_copy_number_table["sample"].unique()
        self.clone_names = list(set(flat_list(tree)))
        self.mrca = self.tree[0][0]
        self.tree_edges = get_tree_edges(self.tree)

        # :::::  outputs:
        self.complexity = None
        self.total_score = None
        self.B = None
        self.A = None
        self.solution = None

        # ::::: initialise model
        print("Initialise model::")
        print(
            f"===={self.tumour_id}_{segment}_allowed_complexity_{self.allowed_tree_complexity}===="
        )
        # activate gurobi license:
        if self.license == "remote":
            options = {
                "WLSACCESSID": os.getenv(
                    "WLSACCESSID", "306b7e47-cdee-4e94-ab85-de6fca0dd709"
                ),
                "WLSSECRET": os.getenv(
                    "WLSSECRET", "c7349a51-4893-47fa-81d4-714b0445a9e8"
                ),
                "LICENSEID": int(os.getenv("LICENSEID", 2507321)),
            }
            env = gp.Env(params=options)
            self.model = gp.Model("ALPACA", env=env)
            print(f'Using remote license: {options["LICENSEID"]}')
        else:
            print(f"Using local license")
            self.model = gp.Model("ALPACA")
        self.model.params.TimeLimit = self.time_limit
        self.model.params.Threads = self.cpus * 2
        if self.BestObjStop:
            self.model.params.BestObjStop = self.BestObjStop
        for allele in ["A", "B"]:
            self.X[allele] = self.model.addVars(
                self.clone_names, name=f"X{allele}", lb=0, vtype=GRB.INTEGER
            )
            self.Y[allele] = {
                row["sample"]: row[f"cpn{allele}"]
                for _, row in self.fractional_copy_number_table.iterrows()
            }
            self.d[allele] = self.model.addVars(
                self.sample_names, name=f"d{allele}", lb=0
            )
            self.Yhat[allele] = self.model.addVars(
                self.sample_names, name=f"Yhat{allele}", lb=0
            )
        # Introduce diploid pseudo-clone:
        if self.minimise_events_to_diploid:
            self.add_diploid_pseudo_clone()

        # ::::: mandatory constraints

        self.add_Yhat_constraints()  # constraints defining predicted fractional copy number values
        self.add_CI_constraints()  # constraints indicating if predicted fractional copy number values are withing CIs
        self.add_absolute_distance_constraint()  # constraint distance to be absolute
        self.add_event_count_variables()  # create variables for events on each edge

        if self.add_event_count_constraints_flag:
            self.add_event_count_constraints()  # count events on each edge (sum of magnitudes of changes)
        else:
            self.total_events = 0
        # ::::: facultative constraints
        self.add_state_change_count_variables()
        self.add_path_variability_penalty_variables()

        if self.add_state_change_count_constraints_flag:
            self.add_state_change_count_constraints()  # count state changes on each edge (sum of number of changes, binary for each edge)
        else:
            self.total_edge_changes = 0

        if self.add_path_variability_penalty_constraints_flag:
            self.add_path_variability_penalty_constraints()  # penalise variability of copy number states (positive vs negative changes) on each path
        else:
            self.total_path_variability_penalty = 0

        if self.add_allow_only_one_non_directional_event_flag:
            self.add_allow_only_one_non_directional_event()

        if self.limit_homozygous_deletions_threshold_flag:
            self.limit_homozygous_deletions_threshold()  # constraint for no homozygous deletions, allowed by default (threshold==0)

        if self.restrict_heterogeneity_flag:
            self.restrict_heterogeneity()  # allow subclonal solutions only if in at least one sample CIs don't overlap with any integer

        if self.prevent_increase_from_zero_flag:
            self.prevent_increase_from_zero()  # prevent increase from zero between parent and child

        if self.restrict_to_clonal_only_flag:
            self.restrict_to_clonal_only()  # restrict solutions to clonal only

        # ::::: define complexity
        self.total_tree_complexity = self.model.addVar(
            name="total_tree_complexity", lb=0, vtype=GRB.CONTINUOUS
        )

        self.model.addConstr(
            self.total_tree_complexity
            == self.total_path_variability_penalty
            + self.total_events
            + self.total_edge_changes,
            name="total_complexity_components",
        )
        # total complexity cannot exceed allowed complexity
        self.model.addConstr(
            self.total_tree_complexity <= self.allowed_tree_complexity,
            name="tree_complexity_constr",
        )

        # ::::: set model sense and define objectives:
        self.model.ModelSense = GRB.MINIMIZE

        self.Z = gp.quicksum(
            [
                self.CI_overlap[allele][sample]
                for allele in ["A", "B"]
                for sample in self.sample_names
            ]
        )

        self.D = gp.quicksum((gp.quicksum(self.d[allele]) for allele in ["A", "B"]))

        # noinspection PyArgumentList
        self.model.setObjectiveN(
            self.Z, index=0, priority=1
        )  # larger values indicate higher priorities
        if self.two_objectives:
            # noinspection PyArgumentList
            self.model.setObjectiveN(self.D, index=1, priority=0)

    def add_Yhat_constraints(self):
        for allele in ["A", "B"]:
            self.model.addConstrs(
                (
                    self.Yhat[allele][sample]
                    == gp.quicksum(
                        (
                            self.clone_proportions[sample].loc[clone]
                            * self.X[allele][clone]
                            for clone in self.clone_names
                        )
                    )
                    for sample in self.sample_names
                ),
                name=f"Yhat{allele}_constr",
            )

    def add_CI_constraints(self):
        """
        Introduce indicator variable to check if Yhat (predicted fractional copy-number) is above the upper CI or below the lower CI:
        Variables:
        Yhat (continuous, predicted fractional copy number)
        L = lower CI (continuous)
        U = upper CI (continuous)
        z = binary indicator variable (1 if Yhat is above upper CI or below lower CI, 0 otherwise)
        M = bigM (large positive number)

        Constraints:
        Yhat >= U - M + M * zu: if Yhat is above upper CI, zu = 1, otherwise zu = 0
        Yhat <= U + M * zu: if Yhat is above upper CI, zu = 1, otherwise zu = 0

        Yhat <= L + M * (1-zl): if Yhat is below lower CI, zl = 1, otherwise zl = 0
        Yhat >= L - M * zl: if Yhat is below lower CI, zl = 1, otherwise zl = 0

        Implementing the OR constraint:

        implementation using gurobi or_() function:
        self.model.addConstrs((self.CI_overlap[allele][sample] == gp.or_([self.yhat_above_upper_CI[allele][sample], self.yhat_below_lower_CI[allele][sample]]) for sample in self.sample_names), name='CI_overlap')

        implementation using binary variables:
        self.CI_overlap[allele][sample] is a binary variable, CIO (is predicted value above upper CI or below lower CI?)
        self.yhat_above_upper_CI[allele][sample] is a binary variable, U (is predicted value above upper CI?)
        self.yhat_below_lower_CI[allele][sample] is a binary variable, L (is predicted value below lower CI?)

        we need to linearize the OR constraint:
        CIO >= z_u
        CIO >= z_l
        CIO <= z_u + z_l

        """
        for allele in ["A", "B"]:
            self.CI_upper[allele] = self.model.addVars(
                self.sample_names,
                name=f"{allele}_CI_upper",
                vtype=GRB.CONTINUOUS,
                lb=float("-inf"),
            )
            self.model.addConstrs(
                (
                    self.CI_upper[allele][sample]
                    == (
                        self.ci_table[
                            (self.ci_table["segment"] == self.segment)
                            & (self.ci_table["sample"] == sample)
                        ][f"upper_CI_{allele}"].iloc[0]
                    )
                    for sample in self.sample_names
                ),
                name="upper_CI",
            )

            self.CI_lower[allele] = self.model.addVars(
                self.sample_names,
                name=f"{allele}_CI_lower",
                vtype=GRB.CONTINUOUS,
                lb=float("-inf"),
            )
            self.model.addConstrs(
                (
                    self.CI_lower[allele][sample]
                    == (
                        self.ci_table[
                            (self.ci_table["segment"] == self.segment)
                            & (self.ci_table["sample"] == sample)
                        ][f"lower_CI_{allele}"].iloc[0]
                    )
                    for sample in self.sample_names
                ),
                name="lower_CI",
            )

            self.yhat_above_upper_CI[allele] = self.model.addVars(
                self.sample_names, name=f"Yhat{allele}_CI_above_upper", vtype=GRB.BINARY
            )
            self.yhat_below_lower_CI[allele] = self.model.addVars(
                self.sample_names, name=f"Yhat{allele}_CI_below_lower", vtype=GRB.BINARY
            )
            self.CI_overlap[allele] = self.model.addVars(
                self.sample_names, name=f"Yhat{allele}_CI_overlap", vtype=GRB.BINARY
            )

            M = 1000
            # Yhat is above upper CI
            self.model.addConstrs(
                (
                    self.Yhat[allele][sample]
                    >= self.CI_upper[allele][sample]
                    - M
                    + M * (self.yhat_above_upper_CI[allele][sample])
                    for sample in self.sample_names
                ),
                name=f"bigM_constr1L_CI_{allele}",
            )

            self.model.addConstrs(
                (
                    self.Yhat[allele][sample]
                    <= self.CI_upper[allele][sample]
                    + (M * self.yhat_above_upper_CI[allele][sample])
                    for sample in self.sample_names
                ),
                name=f"bigM_constr2L_CI_{allele}",
            )

            # Yhat is below lower CI
            self.model.addConstrs(
                (
                    self.Yhat[allele][sample]
                    <= self.CI_lower[allele][sample]
                    + M * (1 - self.yhat_below_lower_CI[allele][sample])
                    for sample in self.sample_names
                ),
                name=f"bigM_constr1U_CI_{allele}",
            )

            self.model.addConstrs(
                (
                    self.Yhat[allele][sample]
                    >= self.CI_lower[allele][sample]
                    - (M * self.yhat_below_lower_CI[allele][sample])
                    for sample in self.sample_names
                ),
                name=f"bigM_constr1U_CI_{allele}",
            )

            self.model.addConstrs(
                (
                    self.CI_overlap[allele][sample]
                    >= self.yhat_above_upper_CI[allele][sample]
                    for sample in self.sample_names
                ),
                name="CI_overlap_U",
            )
            self.model.addConstrs(
                (
                    self.CI_overlap[allele][sample]
                    >= self.yhat_below_lower_CI[allele][sample]
                    for sample in self.sample_names
                ),
                name="CI_overlap_L",
            )
            self.model.addConstrs(
                (
                    self.CI_overlap[allele][sample]
                    <= (
                        self.yhat_above_upper_CI[allele][sample]
                        + self.yhat_below_lower_CI[allele][sample]
                    )
                    for sample in self.sample_names
                ),
                name="CI_overlap_UL",
            )

    def add_absolute_distance_constraint(self):
        for allele in ["A", "B"]:
            self.model.addConstrs(
                (
                    self.d[allele][sample]
                    >= self.Yhat[allele][sample] - self.Y[allele][sample]
                    for sample in self.sample_names
                ),
                name=f"d{allele}_constraint_abs_1",
            )
            self.model.addConstrs(
                (
                    -self.d[allele][sample]
                    <= self.Yhat[allele][sample] - self.Y[allele][sample]
                    for sample in self.sample_names
                ),
                name=f"d{allele}_constraint_abs_2",
            )

    def add_diploid_pseudo_clone(self):
        self.tree_edges.add(("diploid", self.mrca))
        for allele in ["A", "B"]:
            self.X[allele]["diploid"] = self.model.addVar(name=f"X{allele}_diploid")
            self.model.addConstr(self.X[allele]["diploid"] == 1)

    def add_event_count_variables(self):
        if not self.minimise_events_to_diploid:
            # remove events between diploid and mrca - required for correct event count in both scenarios
            self.tree_edges = [edge for edge in self.tree_edges if edge[0] != "diploid"]
        # create variables for events on each edge:

        for allele in ["A", "B"]:
            self.CN_diff_edges_amp[allele] = self.model.addVars(
                self.tree_edges, name=f"amp_{allele}", vtype=GRB.INTEGER, lb=0
            )
            self.CN_diff_edges_del[allele] = self.model.addVars(
                self.tree_edges, name=f"del_{allele}", vtype=GRB.INTEGER, lb=0
            )
            self.cpn_change_up[allele] = self.model.addVars(
                self.tree_edges, name=f"{allele}_cpn_change_up", vtype=GRB.BINARY
            )
            self.cpn_change_down[allele] = self.model.addVars(
                self.tree_edges, name=f"{allele}_cpn_change_down", vtype=GRB.BINARY
            )
            # count events on edges:
            for edge in self.tree_edges:
                parent, child = edge
                self.model.addConstr(
                    self.X[allele][parent] + self.CN_diff_edges_amp[allele][edge]
                    == self.X[allele][child] + self.CN_diff_edges_del[allele][edge],
                    name=f"events_{allele}_{edge}",
                )
                # Constraint below ensures that deletion and amplification cannot be present on the same edge #
                if self.exclusive_amp_del:
                    U = 1000
                    # each edge and each allele can have one change up or down. To reflect this we take the number of events (i.e. the magnitude of change) and binarize it to 0 and 1 states.
                    self.model.addConstr(
                        self.CN_diff_edges_amp[allele][edge]
                        <= U * (self.cpn_change_up[allele][edge]),
                        name=f"Ueps_constr1_{edge}_{allele}_up",
                    )
                    self.model.addConstr(
                        self.CN_diff_edges_amp[allele][edge]
                        >= self.cpn_change_up[allele][edge],
                        name=f"Ueps_constr2_{edge}_{allele}_up",
                    )
                    self.model.addConstr(
                        self.CN_diff_edges_del[allele][edge]
                        <= U * (self.cpn_change_down[allele][edge]),
                        name=f"Ueps_constr1_{edge}_{allele}_down",
                    )
                    self.model.addConstr(
                        self.CN_diff_edges_del[allele][edge]
                        >= self.cpn_change_down[allele][edge],
                        name=f"Ueps_constr2_{edge}_{allele}_down",
                    )
                    self.model.addConstr(
                        self.cpn_change_up[allele][edge]
                        + self.cpn_change_down[allele][edge]
                        <= 1,
                        name=f"exc_amp_del_{allele}_{edge}",
                    )
            # constraint total magnitude of cn change per segment:
            self.n[allele] = self.model.addVar(
                name=f"num_events_{allele}", lb=0, vtype=GRB.INTEGER
            )
            self.model.addConstr(
                self.n[allele]
                == gp.quicksum(
                    (
                        self.CN_diff_edges_amp[allele][edge]
                        + self.CN_diff_edges_del[allele][edge]
                        for edge in self.tree_edges
                    )
                ),
                name=f"num_events_{allele}",
            )
        self.total_events_count = self.model.addVar(
            vtype=GRB.INTEGER, lb=0, name="total_events"
        )
        self.model.addConstr(
            self.total_events_count == gp.quicksum([self.n["A"], self.n["B"]]),
            name="num_SCNA_events_count",
        )

    def add_event_count_constraints(self):
        """
        Event = magnitude of change on edge,
        e.g. +4 = 4 events, -1 = 1 event
        """
        self.total_events = self.model.addVar(
            vtype=GRB.INTEGER, lb=0, name="total_events"
        )
        self.model.addConstr(
            self.total_events == gp.quicksum([self.n["A"], self.n["B"]]),
            name="num_SCNA_events_constr",
        )

    def add_state_change_count_variables(self):
        """
        Given indicator variable z, constants U and eps, and variable x (copy number change on edge):
        Define constraints
        constraint_1: x ≤ U * z
        constraint_2: x ≥ z
        z can only be zero when x is zero
        z can only be one when x is greater than zero
        """
        U = 1000  # U must be higher than any anticipated copy number state
        if not self.minimise_events_to_diploid:
            # remove events between diploid and mrca - required for correct event count in both scenarios
            self.tree_edges = [edge for edge in self.tree_edges if edge[0] != "diploid"]
        for allele in ["A", "B"]:
            self.cpn_change_up[allele] = self.model.addVars(
                self.tree_edges, name=f"{allele}_cpn_change_up", vtype=GRB.BINARY
            )
            self.cpn_change_down[allele] = self.model.addVars(
                self.tree_edges, name=f"{allele}_cpn_change_down", vtype=GRB.BINARY
            )
            for edge in self.tree_edges:
                # each edge and each allele can have one change up or down. To reflect this we take the number of events (i.e. the magnitude of change) and binarize it to 0 and 1 states.
                self.model.addConstr(
                    self.CN_diff_edges_amp[allele][edge]
                    <= U * (self.cpn_change_up[allele][edge]),
                    name=f"Ueps_constr1_{edge}_{allele}_up",
                )
                self.model.addConstr(
                    self.CN_diff_edges_amp[allele][edge]
                    >= self.cpn_change_up[allele][edge],
                    name=f"Ueps_constr2_{edge}_{allele}_up",
                )
                self.model.addConstr(
                    self.CN_diff_edges_del[allele][edge]
                    <= U * (self.cpn_change_down[allele][edge]),
                    name=f"Ueps_constr1_{edge}_{allele}_down",
                )
                self.model.addConstr(
                    self.CN_diff_edges_del[allele][edge]
                    >= self.cpn_change_down[allele][edge],
                    name=f"Ueps_constr2_{edge}_{allele}_down",
                )

        self.total_edge_changes_count = self.model.addVar(
            name="num_SCNA_events_count", lb=0, vtype=GRB.INTEGER
        )  # this component is just to report the number
        self.model.addConstr(
            self.total_edge_changes_count
            == gp.quicksum(
                (
                    gp.quicksum(self.cpn_change_up[allele].values())
                    + gp.quicksum(self.cpn_change_down[allele].values())
                    for allele in ["A", "B"]
                )
            ),
            name="edge_change_count_constr",
        )

    def add_state_change_count_constraints(self):
        self.total_edge_changes = self.model.addVar(
            name="num_SCNA_events", lb=0, vtype=GRB.INTEGER
        )  # this component is used in final complexity calculation
        self.model.addConstr(
            self.total_edge_changes
            == gp.quicksum(
                (
                    gp.quicksum(self.cpn_change_up[allele].values())
                    + gp.quicksum(self.cpn_change_down[allele].values())
                    for allele in ["A", "B"]
                )
            ),
            name="edge_change_count_constr",
        )

    def add_allow_only_one_non_directional_event(self):
        """
        Each path from MRCA to a leaf can only have one event going in the opposite direction to the rest of the path. For example, if there are two gains/amplifications, then we can have only one loss on the same path.

        Binarization of the event count:
        Let x be an integer representing the number of edges with amplification on a path (0, 1, 2, 3, ...).
        Let y be an indicator variable which is 1 if x > 1, and 0 otherwise.
        1st constraint ensures that this is the case when x is 0 or 1:
        x >= 2 - U + U*y
        if x is 0 or 1, then y must be 0, otherwise (if y is 1) -U+U*1 = 0 and neither 0 or 1 is greater/equal to 2
        if x is above 1, then y can be either 0 or 1, and the constraint is satisfied
        2nd constraint ensures that this is the case when x is 2 or more:
        x <= 1 + U*y
        if x is 0 or 1, y can have any value because both 0 and 1 are less or equal to 1
        if x is 2 or more, then y must be 1, because right-hand side must be at least 2, so 1 + non-zero number
        """
        U = 1000
        for allele in ["A", "B"]:
            self.amps_count_on_path[allele] = {}
            self.dels_count_on_path[allele] = {}
            self.more_than_1_amp_change[allele] = {}
            self.more_than_1_del_change[allele] = {}
            for path_index, path in enumerate(self.tree):
                path_edges = find_path_edges(path, self.tree_edges)

                # amps
                self.amps_count_on_path[allele][path_index] = self.model.addVar(
                    name=f"{allele}_{path_index}_amps_count_on_path",
                    lb=0,
                    vtype=GRB.INTEGER,
                )
                self.model.addConstr(
                    self.amps_count_on_path[allele][path_index]
                    == gp.quicksum(
                        [self.cpn_change_up[allele][edge] for edge in path_edges]
                    )
                )

                self.more_than_1_amp_change[allele][path_index] = self.model.addVar(
                    name=f"{allele}_{path_index}_more_than_1_amp_change",
                    vtype=GRB.BINARY,
                )

                # self.model.addGenConstrIndicator(self.more_than_1_amp_change[allele][path_index], True, self.amps_count_on_path[allele][path_index] >= 2, name=f'{allele}_{path_index}_more_than_1_amp_change_ctr')
                self.model.addConstr(
                    self.amps_count_on_path[allele][path_index]
                    >= 2 - U + U * self.more_than_1_amp_change[allele][path_index],
                    f"{allele}_{path_index}_more_than_1_amp_change_ctr",
                )
                self.model.addConstr(
                    self.amps_count_on_path[allele][path_index]
                    <= 1 + U * self.more_than_1_amp_change[allele][path_index],
                    f"{allele}_{path_index}_more_than_1_amp_change_ctr",
                )

                # dels
                self.dels_count_on_path[allele][path_index] = self.model.addVar(
                    name=f"{allele}_{path_index}_dels_count_on_path",
                    lb=0,
                    vtype=GRB.INTEGER,
                )
                self.model.addConstr(
                    self.dels_count_on_path[allele][path_index]
                    == gp.quicksum(
                        [self.cpn_change_down[allele][edge] for edge in path_edges]
                    )
                )

                self.more_than_1_del_change[allele][path_index] = self.model.addVar(
                    name=f"{allele}_{path_index}_more_than_1_del_change",
                    vtype=GRB.BINARY,
                )

                # self.model.addGenConstrIndicator(self.more_than_1_del_change[allele][path_index], True, self.dels_count_on_path[allele][path_index] >= 2, name=f'{allele}_{path_index}_more_than_1_amp_change_ctr')
                self.model.addConstr(
                    self.dels_count_on_path[allele][path_index]
                    >= 2 - U + U * self.more_than_1_del_change[allele][path_index],
                    f"{allele}_{path_index}_more_than_1_del_change_ctr",
                )
                self.model.addConstr(
                    self.dels_count_on_path[allele][path_index]
                    <= 1 + U * self.more_than_1_del_change[allele][path_index],
                    f"{allele}_{path_index}_more_than_1_del_change_ctr",
                )

                # amps and dels cannot have more than one change on the same path:
                self.model.addConstr(
                    self.more_than_1_amp_change[allele][path_index]
                    + self.more_than_1_del_change[allele][path_index]
                    <= 1,
                    f"{allele}_{path_index}_more_than_1_change_ctr",
                )

    def prevent_increase_from_zero(self):
        """
        Each edge has both its parent cn value, and a binary variable indicating whether the change on the edge was positive.
        To prevent increase from zero, we introduce a constraint that parent cn value must be greater or equal to the binary variable representing positive change on the edge.
        Therefore, if there is any gain on the edge (binary variable == 1), parent must also be at least 1. If there is no gain (binary variable == 0), parent cn value can be zero.
        """
        for allele in ["A", "B"]:
            for edge in self.tree_edges:
                parent, _ = edge
                self.model.addConstr(
                    self.X[allele][parent] >= self.cpn_change_up[allele][edge],
                    name=f"prevent_increase_from_zero_{edge}_{allele}",
                )

    def add_path_variability_penalty_variables(self):
        """
        Introduce variables to represent each path and each allele
        For each path and each allele, total additional cost equals: variability_penalty * number_of_edges_with_positive_cpn_change * number_of_edges_with_negative_cpn_change
        This approach allows for different states on different paths, but penalizes paths with more state changes
        """
        for allele in ["A", "B"]:
            self.path_variability_penalty[allele] = self.model.addVars(
                range(0, len(self.tree)),
                name=f"path_variability_penalty_{allele}",
                vtype=GRB.INTEGER,
            )
            for path_index, path in enumerate(self.tree):
                path_edges = find_path_edges(path, self.tree_edges)
                self.model.addConstr(
                    self.path_variability_penalty[allele][path_index]
                    == (
                        self.variability_penalty
                        * gp.quicksum(
                            [self.cpn_change_up[allele][edge] for edge in path_edges]
                        )
                        * gp.quicksum(
                            [self.cpn_change_down[allele][edge] for edge in path_edges]
                        )
                    )
                )
        self.total_path_variability_penalty_count = self.model.addVar(
            name="total_path_variability_penalty_count", lb=0, vtype=GRB.INTEGER
        )
        self.model.addConstr(
            self.total_path_variability_penalty_count
            == (
                gp.quicksum(self.path_variability_penalty["A"].values())
                + gp.quicksum(self.path_variability_penalty["B"].values())
            ),
            name="path_variability_constr",
        )

    def add_path_variability_penalty_constraints(self):
        self.total_path_variability_penalty = self.model.addVar(
            name="total_path_variability_penalty", lb=0, vtype=GRB.INTEGER
        )
        self.model.addConstr(
            self.total_path_variability_penalty
            == (
                gp.quicksum(self.path_variability_penalty["A"].values())
                + gp.quicksum(self.path_variability_penalty["B"].values())
            ),
            name="path_variability_constr",
        )

    def limit_homozygous_deletions_threshold(self):
        # clones can have homozygous deletion if they or their descendants are present in a sample with fractional copy number below threshold and is segment is below size limit
        threshold = float(self.homozygous_deletion_threshold)
        thr_size = int(self.homo_del_size_limit)
        seg_len = get_length_from_name(self.segment)
        if seg_len < thr_size:
            samples_with_low_fractional_A = [
                k[0] for k in self.Y["A"].items() if k[1] < threshold
            ]
            samples_with_low_fractional_B = [
                k[0] for k in self.Y["B"].items() if k[1] < threshold
            ]
            samples_where_homozygous_deletion_is_permitted = list(
                set(
                    [
                        x
                        for x in samples_with_low_fractional_A
                        + samples_with_low_fractional_B
                        if (x in samples_with_low_fractional_A)
                        and (x in samples_with_low_fractional_B)
                    ]
                )
            )
            all_clones = list(self.clone_proportions.index)
            clones_present_in_these_samples = list(
                self.clone_proportions.index[
                    self.clone_proportions[
                        samples_where_homozygous_deletion_is_permitted
                    ].sum(axis=1)
                    > 0
                ]
            )
            clones_not_present_in_these_samples = list(
                self.clone_proportions.index[
                    self.clone_proportions[
                        samples_where_homozygous_deletion_is_permitted
                    ].sum(axis=1)
                    == 0
                ]
            )
            # exclude absent clones:
            absent_clones = (
                self.clone_proportions.sum(axis=1)[
                    self.clone_proportions.sum(axis=1) == 0
                ]
            ).index
            clones_not_present_in_these_samples = [
                c for c in clones_not_present_in_these_samples if c not in absent_clones
            ]
            for clone in clones_not_present_in_these_samples:
                self.model.addConstr(
                    self.X["A"][clone] + self.X["B"][clone] >= 1,
                    name=f"no_homo_del_{clone}",
                )
        else:
            self.model.addConstrs(
                (
                    self.X["A"][clone] + self.X["B"][clone] >= 1
                    for clone in self.clone_names
                ),
                name="no_homo_del_{clone}",
            )
        return self.model

    def restrict_heterogeneity(self):
        # lenient criterion: if at least one sample containing clone x has no overlap with an integer, clone x can be heterogeneous
        homogenous_clones_A = []
        heterogeneous_clones_A = []
        homogenous_clones_B = []
        heterogeneous_clones_B = []
        for clone in self.clone_names:
            samples_where_clone_is_present = (
                (self.clone_proportions[self.clone_proportions.index == clone] > 0)
                .iloc[0][
                    (
                        self.clone_proportions[self.clone_proportions.index == clone]
                        > 0
                    ).iloc[0]
                ]
                .index.to_list()
            )
            clone_ci_table = self.ci_table[
                self.ci_table["sample"].isin(samples_where_clone_is_present)
            ]
            clone_ci_table = clone_ci_table[clone_ci_table.segment == self.segment]
            clone_ci_table["span_A"] = clone_ci_table["upper_CI_A"].astype(
                int
            ) - clone_ci_table["lower_CI_A"].astype(
                int
            )  # 0 = no overlap with integer, 1 = overlap with one integer, 2 = overlap with two integers etc
            clone_ci_table["span_B"] = clone_ci_table["upper_CI_B"].astype(
                int
            ) - clone_ci_table["lower_CI_B"].astype(int)
            if (
                clone_ci_table["span_A"] == 0
            ).any():  # if at least one equals zero it means that in at least one sample CI don't overlap with any integer
                heterogeneous_clones_A.append(clone)
            elif (clone_ci_table["span_A"] == 1).any():
                homogenous_clones_A.append(clone)
                cn_value = int(
                    (
                        clone_ci_table["upper_CI_A"].astype(int)[
                            clone_ci_table["span_A"] == 1
                        ]
                    ).median()
                )
                # constrain to integer copy number:
                self.model.addConstr(
                    self.X["A"][clone] == cn_value,
                    name=f"homogenous_integer_constr_{clone}_A",
                )
            if (clone_ci_table["span_B"] == 0).any():
                heterogeneous_clones_B.append(clone)
            elif (clone_ci_table["span_B"] == 1).any():
                homogenous_clones_B.append(clone)
                cn_value = int(
                    (
                        clone_ci_table["upper_CI_B"].astype(int)[
                            clone_ci_table["span_B"] == 1
                        ]
                    ).median()
                )
                # constrain to integer copy number:
                self.model.addConstr(
                    self.X["B"][clone] == cn_value,
                    name=f"homogenous_integer_constr_{clone}_B",
                )

    def restrict_to_clonal_only(self):
        if "diploid" in self.clone_names:
            print(
                "Error - remove diploid from clone lists before applying this constraint"
            )
        number_of_clones = len(self.clone_names)
        for clone in self.clone_names:
            for allele in ["A", "B"]:
                self.model.addConstr(
                    (number_of_clones * self.X[allele][clone])
                    == (gp.quicksum(self.X[allele])),
                    name=f"clonal_{allele}",
                )

    def get_output(self):
        A = pd.DataFrame(
            {
                c: [int(round(cn_val.X))]
                for c, cn_val in self.X["A"].items()
                if c[0] != "diploid"
            },
            index=["pred_CN_A"],
        ).T
        B = pd.DataFrame(
            {
                c: [int(round(cn_val.X))]
                for c, cn_val in self.X["B"].items()
                if c[0] != "diploid"
            },
            index=["pred_CN_B"],
        ).T
        solution = pd.merge(A, B, left_index=True, right_index=True)
        self.A = A
        self.B = B
        self.total_score = self.Z.getValue() + self.D.getValue()
        self.complexity = int(self.total_tree_complexity.X)
        solution["complexity"] = self.complexity
        solution["CI_score"] = int(self.Z.getValue())
        solution["D_score"] = round(self.D.getValue(), 3)
        solution["variability_penalty_count"] = int(
            self.total_path_variability_penalty_count.X
        )
        solution["state_change_count"] = int(self.total_edge_changes_count.X)
        solution["event_count"] = int(self.total_events_count.X)
        solution["allowed_complexity"] = self.allowed_tree_complexity
        solution.index.name = "clone"
        solution.reset_index(inplace=True)
        self.solution = solution

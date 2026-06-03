import argparse
import ast
import os
import sys
from alpaca.ALPACA_model_class import Model as ALPACA_Model
from alpaca.utils import SUPPORTED_GENOME_BUILDS

"""
ALPACA can operate in two modes: 'tumour' and 'segment'.
In default 'tumour' mode, expected input is a directory with all the inputs required for a single tumour:
fractional copy numbers, clone proportions, confidence intervals and tree.
In this mode, ALPACA will write a single output file for each tumour.
In 'segment' mode, expected input is an array of files to segment files (can be from different tumours).
In this mode, ALPACA will create separate outputs for each segment.
In this mode, ALPACA also requires the input_data_directory to be specified: it is the parent directory
containing subdirectories for each tumour. ALPACA will automatically find the correct tumour subdirectory
by parsing the name of each of the input files.
"""


def get_parser():
    parser = argparse.ArgumentParser(
        description="Run ALPACA with specified parameters."
    )
    # SEGMENT mode arguments:
    parser.add_argument(
        "--input_data_directory",
        type=str,
        required=False,
        help="Directory where input data is stored. Should contain subdirectories for each tumour",
    )
    parser.add_argument(
        "--input_files",
        nargs="+",
        type=str,
        required=False,
        help="Space-separated list of input tables for one or multiple segments.",
    )
    # TUMOUR mode arguments:
    parser.add_argument(
        "--input_tumour_directory",
        type=str,
        required=False,
        help="Directory with all the inputs required for a single tumour.",
    )
    # Common arguments:
    parser.add_argument(
        "--mode",
        type=str,
        default="tumour",
        help="Mode of operation. If 'tumour', expect single file with all the segments and output a single file.\
            If 'segment' expect array of files to segment files (can be from different tumours) and create separate outputs for each segment.",
    )
    parser.add_argument(
        "--nextflow_config",
        type=str,
        default="",
        help="Path to a Nextflow config file to run segment-level distributed workers (only used when --mode segment).",
    )
    parser.add_argument(
        "--overwrite_output",
        type=int,
        default=1,
        help="If set to 0, ALPACA will check if solution file for each segment already exits, and will skip iteration if it does. \
            If set to 1, ALPACA will overwrite existing solution files.",
    )
    parser.add_argument(
        "--output_directory",
        type=str,
        default="./",
        help="Directory where output data is stored. Defaults to current directory.",
    )
    parser.add_argument(
        "--plot_output_mode",
        type=str,
        choices=["pdf", "notebook", "none"],
        default="notebook",
        help="Control how ALPACA emits visualisations after a run: 'pdf' saves static PDFs, 'notebook' (default) writes an interactive Jupyter notebook, 'none' skips plotting.",
    )
    parser.add_argument(
        "--genome_build",
        type=str,
        choices=SUPPORTED_GENOME_BUILDS,
        default="hg19",
        help="Reference genome build for chromosome lengths when plotting (default: hg19).",
    )
    parser.add_argument(
        "--gurobi_logs",
        type=str,
        default="",
        help="Directory where gurobi logs should be stored. If no value is speficied, logs will not be saved.",
    )
    parser.add_argument(
        "--solver",
        type=str,
        default="pyomo",
        help="Solver backend to use (options: pyomo, gurobi).",
    )
    parser.add_argument(
        "--pyomo_solver",
        type=str,
        default="scip",
        help="Pyomo solver plugin to use when --solver pyomo (e.g. scip, glpk, gurobi).",
    )
    parser.add_argument(
        "--solver_logs",
        type=str,
        default="",
        help="Optional log directory for solver-agnostic backends. Falls back to --gurobi_logs when empty.",
    )
    parser.add_argument(
        "--pyomo_solver_options",
        nargs="+",
        default=[],
        help="Additional Pyomo solver options in key=value form (e.g. ratioGap=1e-5).",
    )
    parser.add_argument(
        "--objectives",
        type=str,
        default="DCI",
        help="Which objectives to include: 'D' (distance only), 'CI' (confidence-interval violations only), or 'DCI' (both, default).",
    )
    parser.add_argument(
        "--minimise_events_to_diploid",
        type=int,
        default=1,
        help="Whether to minimize events to diploid or not. If true, ALPACA will introduce diploid pseudo-clone at the root of the tree.",
    )
    parser.add_argument(
        "--prevent_increase_from_zero_flag", type=int, default=1, help=""
    )
    parser.add_argument(
        "--add_event_count_constraints_flag", type=int, default=1, help=""
    )
    parser.add_argument(
        "--add_allow_only_one_non_directional_event_flag",
        type=int,
        default=1,
        help="If true only 1 type of change (negative change (loss) or positivechange (gain)) can happen multiple times on a single path from MRCA to a leaf.",
    )
    parser.add_argument(
        "--time_limit",
        default=300,
        type=int,
        help="Time limit in seconds for each model run",
    )
    parser.add_argument(
        "--homozygous_deletion_threshold",
        type=float,
        default=1,
        help="Model will be allowed to postulate homozygous deletion only if fractional copy number per segment is below this value.",
    )
    parser.add_argument(
        "--homo_del_size_limit",
        type=float,
        default=50000000,
        help="Model will be allowed to postulate homozygous deletion only if width of segment is below this value.",
    )
    parser.add_argument(
        "--missing_clones_inherit_from_children_flag",
        default=1,
        type=int,
        help="Ensure that missing clones inherit cn from childen (events go up in the tree)",
    )
    parser.add_argument(
        "--d_zero",
        default=0,
        type=int,
        help="Set to true if objective function is expected to reach value of zero (only in certain simulated scenarios)",
    )
    parser.add_argument("--cpus", default=1, type=int, help="number of available cpus")
    parser.add_argument("--rsc", default=0, type=int, help="remove small clones")
    parser.add_argument(
        "--ccp", default=0, type=int, help="calibrate clone proportions"
    )
    parser.add_argument(
        "--ci_table_name",
        default="ci_table.csv",
        type=str,
        help="Name of file containing confidence intervals for SNP copynumbers",
    )
    parser.add_argument(
        "--min_ci",
        type=float,
        default=0.1,
        help="Minimum allowed confidence interval span for each allele (float). CI spans tighter than this will be expanded.",
    )
    parser.add_argument("--debug", default=False, action="store_true")
    parser.add_argument(
        "--debug_solution_file",
        type=str,
        default="",
        help="Path to CSV with debug solution values (columns: tumour_id,segment,clone,pred_CN_A,pred_CN_B)",
    )
    parser.add_argument(
        "--output_all_solutions",
        default=False,
        action="store_true",
        help="If set, write all model solutions (not only the optimal) into a subdirectory of the output directory",
    )

    parser.add_argument(
        "--complexity",
        type=int,
        default=None,
        help="If provided together with --debug_solution_file, run a single iteration with this complexity",
    )
    parser.add_argument(
        "--extra_columns",
        nargs="+",
        type=str,
        default=[],
        help="List of extra columns to include in the output. Options: 'gurobi_gap', 'gurobi_time', 'complexity', 'CI_score', 'D_score'. "
        "Note: requesting 'gurobi_time' or 'gurobi_gap' will automatically include objective-specific metrics (e.g. _D, _CI) if available.",
    )
    parser.add_argument(
        "--strict_gap",
        type=int,
        choices=[0, 1],
        default=1,
        help="If set to 1 (default), enforce zero gap tolerance for Gurobi optimization to improve reproducibility. "
        "When enabled, the solver will only stop when it proves optimality (gap=0) or hits the time limit.",
    )
    parser.add_argument(
        "--simulate_infeasibility",
        type=str,
        default="",
        help="If set to a segment name, add an impossible constraint to that segment's model to force infeasibility. "
        "Used to debug infeasible-solution handling.",
    )

    return parser


def _coerce_solver_option_value(raw: str):
    token = raw.strip()
    lowered = token.lower()
    if lowered == "true":
        return True
    if lowered == "false":
        return False
    if lowered == "none":
        return None
    try:
        return ast.literal_eval(token)
    except (ValueError, SyntaxError):
        return token


def _parse_solver_option_list(values):
    options = {}
    for item in values or []:
        if "=" not in item:
            print(
                f"Error: Invalid --pyomo_solver_options entry '{item}'. Use key=value format.",
                file=sys.stderr,
            )
            sys.exit(1)
        key, value = item.split("=", 1)
        key = key.strip()
        if not key:
            print(
                f"Error: Invalid --pyomo_solver_options entry '{item}'. Option name cannot be empty.",
                file=sys.stderr,
            )
            sys.exit(1)
        options[key] = _coerce_solver_option_value(value)
    return options


def validate_args(args):
    if args.mode == "tumour":
        if not args.input_tumour_directory:
            print(
                "Error: --input_tumour_directory is required when --mode is 'tumour'.",
                file=sys.stderr,
            )
            sys.exit(1)
    else:  # mode is 'segment'
        if not args.input_data_directory or not args.input_files:
            print(
                "Error: --input_data_directory and --input_files are required when --mode is 'segment'.",
                file=sys.stderr,
            )
            sys.exit(1)


def make_config(args_in):
    ENV = os.environ.get("APP_ENV", "prod").lower()
    parser = get_parser()
    # If called with sys.argv[1:] the first element may be a sub-command
    # (e.g., 'run'). Strip a leading command token if present so argparse
    # receives only option/flag arguments.
    if args_in and not args_in[0].startswith("-"):
        args_list = args_in[1:]
    else:
        args_list = args_in
    args, remaining_args = parser.parse_known_args(args_list)
    validate_args(args)
    pyomo_solver_options = _parse_solver_option_list(args.pyomo_solver_options)
    # If there are any unknown args left over, fail fast (unless running in dev
    # mode where remaining args are forwarded to dev.parse_optional_args()).
    # also, first argument is a command, so skip that:
    if remaining_args[1:] and ENV not in ("dev", "development"):
        print(
            "Error: Unknown or unsupported arguments passed to ALPACA:\n  "
            + " ".join(remaining_args),
            file=sys.stderr,
        )
        sys.exit(1)
    # make config dictionary
    # Start from Model defaults to avoid duplicating default values across files
    model_config = ALPACA_Model.default_model_config()

    # Override defaults with CLI args and map CLI arg names to Model property names
    model_config.update(
        {
            "solver": args.solver,
            "pyomo_solver": args.pyomo_solver,
            "objectives": args.objectives,
            "minimise_events_to_diploid": bool(args.minimise_events_to_diploid),
            "prevent_increase_from_zero_flag": bool(
                args.prevent_increase_from_zero_flag
            ),
            "add_event_count_constraints_flag": bool(
                args.add_event_count_constraints_flag
            ),
            "add_allow_only_one_non_directional_event_flag": bool(
                args.add_allow_only_one_non_directional_event_flag
            ),
            "homozygous_deletion_threshold": args.homozygous_deletion_threshold,
            "homo_del_size_limit": args.homo_del_size_limit,
            "time_limit": args.time_limit,
            "cpus": args.cpus,
            "gurobi_logs": args.gurobi_logs,
            "solver_logs": args.solver_logs or args.gurobi_logs,
            "missing_clones_inherit_from_children_flag": args.missing_clones_inherit_from_children_flag,
            "d_zero": args.d_zero,
            "debug": args.debug,
            "debug_solution_file": args.debug_solution_file,
            "complexity": args.complexity,
            "strict_gap": bool(args.strict_gap),
            "pyomo_solver_options": pyomo_solver_options,
            "simulate_infeasibility": args.simulate_infeasibility,
        }
    )
    preprocessing_config = {
        "mode": args.mode,
        "overwrite_output": args.overwrite_output,
        "ci_table_name": args.ci_table_name,
        "debug": args.debug,
        "output_all_solutions": args.output_all_solutions,
        "env": ENV,
        "output_directory": args.output_directory,
        "min_ci": args.min_ci,
        "extra_columns": args.extra_columns,
        "plot_output_mode": args.plot_output_mode,
        "genome_build": args.genome_build,
    }
    if args.mode == "tumour":
        preprocessing_config["input_tumour_directory"] = args.input_tumour_directory
    else:
        preprocessing_config["input_files"] = args.input_files
        preprocessing_config["input_data_directory"] = args.input_data_directory
        preprocessing_config["input_tumour_directory"] = args.input_tumour_directory
    if ENV == "dev":
        print("Starting ALPACA in development mode")
        from dev.parse_optional_args import get_config

        _model_config, _preprocessing_config = get_config(remaining_args)
        model_config.update(_model_config)
        preprocessing_config.update(_preprocessing_config)

    config = {
        "model_config": model_config,
        "preprocessing_config": preprocessing_config,
    }
    return config

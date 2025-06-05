import argparse
import os
import sys

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
        "--gurobi_logs",
        type=str,
        default="",
        help="Directory where gurobi logs should be stored. If no value is speficied, logs will not be saved.",
    )
    parser.add_argument(
        "--use_two_objectives",
        type=int,
        default=False,
        help="Whether to use two objectives or not. First objective minimises number of segments outside CI, second objective minimises error.",
    )
    parser.add_argument(
        "--use_minimise_events_to_diploid",
        type=int,
        default=False,
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
        default=60,
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
    parser.add_argument("--debug", default=False, action="store_true")

    return parser


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
    args, remaining_args = parser.parse_known_args(args_in)
    validate_args(args)
    # make config dictionary
    model_config = {
        "use_two_objectives": args.use_two_objectives,
        "use_minimise_events_to_diploid": args.use_minimise_events_to_diploid,
        "prevent_increase_from_zero_flag": args.prevent_increase_from_zero_flag,
        "add_event_count_constraints_flag": args.add_event_count_constraints_flag,
        "add_allow_only_one_non_directional_event_flag": args.add_allow_only_one_non_directional_event_flag,
        "homozygous_deletion_threshold": args.homozygous_deletion_threshold,
        "homo_del_size_limit": args.homo_del_size_limit,
        "time_limit": args.time_limit,
        "cpus": args.cpus,
        "gurobi_logs": args.gurobi_logs,
    }
    preprocessing_config = {
        "mode": args.mode,
        "overwrite_output": args.overwrite_output,
        "ci_table_name": args.ci_table_name,
        "debug": args.debug,
        "env": ENV,
        "output_directory": args.output_directory,
    }
    if args.mode == "tumour":
        preprocessing_config["input_tumour_directory"] = args.input_tumour_directory
    else:
        preprocessing_config["input_files"] = (args.input_files[0].strip().split(" "),)
        preprocessing_config["input_data_directory"] = args.input_data_directory
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

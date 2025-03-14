import argparse
import os


def get_parser():
    parser = argparse.ArgumentParser(
        description="Run ALPACA with specified parameters."
    )
    parser.add_argument(
        "--input_data_directory",
        type=str,
        required=True,
        help="Directory where input data is stored. Should contain subdirectories for each tumour",
    )

    parser.add_argument(
        "--output_directory",
        type=str,
        default="./",
        help="Directory where output data is stored. Defaults to current directory.",
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
        "--input_files",
        nargs="+",
        type=str,
        required=True,
        help="Space-separated list of input tables for one or multiple segments.",
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


def make_config(args_in):
    ENV = os.environ.get("APP_ENV", "prod").lower()
    parser = get_parser()
    args, remaining_args = parser.parse_known_args(args_in)
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
    }

    preprocessing_config = {
        "ci_table_name": args.ci_table_name,
        "input_data_directory": args.input_data_directory,
        "debug": args.debug,
        "env": ENV,
        "input_files": args.input_files[0].strip().split(" "),
        "output_directory": args.output_directory,
    }

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

import argparse


def get_parser():
    parser = argparse.ArgumentParser(description="Optional / development arguments.")
    parser.add_argument(
        "--use_binary_search",
        type=int,
        default=False,
        help="Whether to use binary search or not. Binary search is faster but may not find optimal solution.",
    )

    parser.add_argument("--exclusive_amp_del", type=int, default=1, help="")

    parser.add_argument(
        "--add_state_change_count_constraints_flag", type=int, default=0, help=""
    )
    parser.add_argument(
        "--add_path_variability_penalty_constraints_flag", type=int, default=0, help=""
    )
    parser.add_argument("--s_type", default="s_strictly_decreasing", type=str)
    parser.add_argument("--compare_with_true_solution", default=0, type=int)
    parser.add_argument("--run_with_qc", default=0, type=int)
    parser.add_argument("--output_all_solutions", default=0, type=int)
    parser.add_argument("--output_model_selection_table", default=0, type=int)
    parser.add_argument(
        "--ci", default=0.9, type=float, help="Confidence interval for SNP copynumbers"
    )


def get_config(remaining_args):
    parser = get_parser()
    args = parser.parse_args(remaining_args)
    # make config dictionary
    model_config = {
        "add_state_change_count_constraints_flag": args.add_state_change_count_constraints_flag,
        "add_path_variability_penalty_constraints_flag": args.add_path_variability_penalty_constraints_flag,
        "exclusive_amp_del": args.exclusive_amp_del,
        "use_binary_search": args.use_binary_search,
    }
    preprocessing_config = {
        "rsc": args.rsc,
        "ccp": args.ccp,
        "ci": args.ci,
        "compare_with_true_solution": args.compare_with_true_solution,
        "s_type": args.s_type,
        "run_with_qc": args.run_with_qc,
        "output_all_solutions": args.output_all_solutions,
    }
    return model_config, preprocessing_config

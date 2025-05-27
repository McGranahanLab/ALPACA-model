import os
import subprocess
import sys
from importlib.resources import files
from alpaca.analysis import get_cn_change_to_ancestor
from alpaca.analysis import calculate_ccd
import argparse
from datetime import datetime
import logging
from alpaca.utils import create_logger, save_dataframe_to_csv


def input_conversion():
    """
    Wrapper function to execute input_conversion.sh from submodule
    """
    print("Running input_conversion - it may take a few minutes")
    try:
        # check if all file in sys.argv exist:
        for i, arg in enumerate(sys.argv[1:], start=1):
            if "=" not in arg and ("/" in arg or "\\" in arg or "." in arg):
                exists = os.path.exists(arg)
                print(
                    f"Argument {i} ({arg}): {'Exists' if exists else 'DOES NOT EXIST'}"
                )
        # Locate the input_conversion.sh script
        script_path = str(
            files("alpaca").joinpath(
                "scripts/submodules/alpaca_input_formatting/input_conversion.sh"
            )
        )
        print(script_path)
        # Locate the submodules directory
        submodules_path = str(files("alpaca").joinpath("scripts/submodules"))
        # Set environment variable to help script locate its dependencies
        env = os.environ.copy()
        env["SUBMODULES_PATH"] = submodules_path
        cmd_args = sys.argv[1:]
        print(f"Command line arguments: {cmd_args}")
        # Execute the shell script with all passed arguments
        result = subprocess.run(
            [script_path] + sys.argv[2:],
            check=True,
            text=True,
            capture_output=True,
            env=env,
        )
        print(f"Return code: {result.returncode}")
        print(f"Output: {result.stdout}")
        print(f"Stderr: {result.stderr}")
        print("")
        print("")
        print("---------------------------------------")
        print("Input conversion completed successfully.")
        return 0
    except subprocess.CalledProcessError as e:
        print(f"Error executing input_conversion: {e.stderr}", file=sys.stderr)
        return e.returncode
    except Exception as e:
        print(f"Unexpected error: {e}", file=sys.stderr)
        return 1


def run_get_cn_change_to_ancestor():
    """CLI wrapper for get_cn_change_to_ancestor"""
    logger = create_logger(name="ancestor-delta", log_dir="logs")
    parser = argparse.ArgumentParser(
        description="Compute copy number changes to ancestor and save to CSV."
    )
    parser.add_argument("command", choices=["ancestor-delta"], help="Command to run")
    parser.add_argument("--tree_path", help="Path to the tree file", required=True)
    parser.add_argument(
        "--tumour_df_path",
        help="Path to the tumour dataframe file (CSV format)",
        required=True,
    )
    parser.add_argument(
        "--output_directory",
        help="Directory to save the output CSV file",
        required=True,
    )

    args = parser.parse_args()
    # Validate input files exist
    if not os.path.isfile(args.tree_path):
        logger.error(f"Tree file not found: {args.tree_path}")
        exit(1)

    if not os.path.isfile(args.tumour_df_path):
        logger.error(f"Tumour dataframe file not found: {args.tumour_df_path}")
        exit(1)

    try:
        logger.info("Starting analysis...")
        cn_change_to_ancestor_df = get_cn_change_to_ancestor(
            args.tree_path, args.tumour_df_path
        )
        save_dataframe_to_csv(
            df=cn_change_to_ancestor_df,
            output_dir=args.output_directory,
            output_filename="cn_change_to_ancestor.csv",
        )
        logger.info(
            f"Analysis completed successfully. Output saved to: {args.output_directory}"
        )

    except Exception as e:
        logger.exception(f"An error occurred during analysis: {e}")
        exit(1)


def run_calculate_ccd():
    """CLI wrapper for calculate_ccd"""
    logger = create_logger(name="ccd_analysis", log_dir="logs")
    parser = argparse.ArgumentParser(
        description="Compute clone copy number diversity and save results."
    )
    parser.add_argument("command", choices=["ccd"], help="Command to run")
    parser.add_argument(
        "--alpaca_output_path",
        help="Path to the results dataframe file (CSV format), either the entire cohort or a single tumour",
        required=True,
    )
    parser.add_argument(
        "--output_directory", help="Path to save the output CSV file", required=True
    )

    args = parser.parse_args()
    # Validate input files exist
    if not os.path.isfile(args.alpaca_output_path):
        logger.error(f"Tumour dataframe file not found: {args.alpaca_output_path}")
        exit(1)
    # check if first row of the file contains columns 'tumour_id' and 'pred_CN_A' and 'pred_CN_B':
    with open(args.alpaca_output_path, "r") as f:
        header = f.readline().strip().split(",")
        required_columns = ["tumour_id", "clone", "segment", "pred_CN_A", "pred_CN_B"]
        missing_columns = [col for col in required_columns if col not in header]
        if missing_columns:
            logger.error(
                f"Tumour dataframe file does not contain required columns: {missing_columns}"
            )
            exit(1)
    try:
        logger.info("Starting CCD analysis...")
        ccd_scores_df = calculate_ccd(args.alpaca_output_path)

        # Ensure output directory exists
        output_dir = args.output_directory
        os.makedirs(output_dir, exist_ok=True)
        output_name = f"{output_dir}/clone_copy_number_diversity_scores.csv"
        ccd_scores_df.to_csv(output_name, index=False)
        logger.info(f"Analysis completed successfully. Output saved to: {output_name}")

    except Exception as e:
        logger.exception(f"An error occurred during analysis: {e}")
        exit(1)

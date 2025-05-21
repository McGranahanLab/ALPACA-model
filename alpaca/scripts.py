import os
import subprocess
import sys
from importlib.resources import files
from alpaca.analysis import get_cn_change_to_ancestor
import argparse

import logging

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)


def input_conversion():
    """
    Wrapper function to execute input_conversion.sh from submodule
    """
    print("Running input_conversion - it may take a few minutes")
    try:
        # check if all file in sys.argv exist:
        for i, arg in enumerate(sys.argv[1:], start=1):
            if '=' not in arg and ('/' in arg or '\\' in arg or '.' in arg):
                exists = os.path.exists(arg)
                print(f"Argument {i} ({arg}): {'Exists' if exists else 'DOES NOT EXIST'}")
        # Locate the input_conversion.sh script
        script_path = str(files("alpaca").joinpath("scripts/submodules/alpaca_input_formatting/input_conversion.sh"))
        print(script_path)
        # Locate the submodules directory
        submodules_path = str(files("alpaca").joinpath("scripts/submodules"))

        # Ensure the script is executable
        # os.chmod(script_path, 0o755)

        # Set environment variable to help script locate its dependencies
        env = os.environ.copy()
        env["SUBMODULES_PATH"] = submodules_path

        # Execute the shell script with all passed arguments
        result = subprocess.run(
            [script_path] + sys.argv[1:],
            check=True,
            text=True,
            capture_output=True,
            env=env,
        )
        print(f"Return code: {result.returncode}")
        print(f"Output: {result.stdout}")
        print(f"Stderr: {result.stderr}")
        return 0
    except subprocess.CalledProcessError as e:
        print(f"Error executing input_conversion: {e.stderr}", file=sys.stderr)
        return e.returncode
    except Exception as e:
        print(f"Unexpected error: {e}", file=sys.stderr)
        return 1


def run_get_cn_change_to_ancestor():
    """CLI wrapper for get_cn_change_to_ancestor"""
    parser = argparse.ArgumentParser(
        description="Compute copy number changes to ancestor and save to CSV."
    )
    parser.add_argument("--tree_path", help="Path to the tree file", required=True)
    parser.add_argument(
        "--tumour_df_path",
        help="Path to the tumour dataframe file (CSV format)",
        required=True,
    )
    parser.add_argument(
        "--output_directory", help="Directory to save the output CSV file", required=True
    )

    args = parser.parse_args()
    # Validate input files exist
    if not os.path.isfile(args.tree_path):
        logging.error(f"Tree file not found: {args.tree_path}")
        exit(1)

    if not os.path.isfile(args.tumour_df_path):
        logging.error(f"Tumour dataframe file not found: {args.tumour_df_path}")
        exit(1)

    try:
        logging.info("Starting analysis...")
        cn_change_to_ancestor_df = get_cn_change_to_ancestor(
            args.tree_path, args.tumour_df_path
        )

        # Ensure output directory exists
        output_dir = args.output_directory
        if output_dir and not os.path.exists(output_dir):
            os.makedirs(output_dir)
        output_name = f"{args.output_directory}/cn_change_to_ancestor.csv"
        cn_change_to_ancestor_df.to_csv(output_name, index=False)
        logging.info(f"Analysis completed successfully. Output saved to: {output_name}")

    except Exception as e:
        logging.exception(f"An error occurred during analysis: {e}")
        exit(1)

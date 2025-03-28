import os
import subprocess
import sys
import pkg_resources
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
        # Locate the input_conversion.sh script
        script_path = pkg_resources.resource_filename(
            "alpaca", "scripts/submodules/alpaca_input_formatting/input_conversion.sh"
        )

        # Locate the submodules directory
        submodules_path = pkg_resources.resource_filename(
            "alpaca", "scripts/submodules"
        )

        # Ensure the script is executable
        os.chmod(script_path, 0o755)

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
        print(result.stdout)
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
        "--output_path", help="Path to save the output CSV file", required=True
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
        output_dir = os.path.dirname(args.output_path)
        if output_dir and not os.path.exists(output_dir):
            os.makedirs(output_dir)

        cn_change_to_ancestor_df.to_csv(args.output_path, index=False)
        logging.info(
            f"Analysis completed successfully. Output saved to: {args.output_path}"
        )

    except Exception as e:
        logging.exception(f"An error occurred during analysis: {e}")
        exit(1)


if __name__ == "__main__":
    if sys.argv[1] == "input_conversion":
        input_conversion()
    elif sys.argv[1] == "get_cn_change_to_ancestor":
        run_get_cn_change_to_ancestor()

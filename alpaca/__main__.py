#!/usr/bin/env python3
import os
import sys
from tqdm import tqdm
from io import StringIO
from alpaca.ALPACA_segment_solution_class import SegmentSolution
from alpaca.utils import (
    show_version,
    show_help,
    print_logo,
    concatenate_output,
    set_run_mode,
    create_logger,
    save_dataframe_to_csv,
)
from alpaca.make_configuration import make_config
from alpaca.analysis import get_cn_change_to_ancestor
import alpaca.scripts as scripts


def main():

    if len(sys.argv) > 1:
        command = sys.argv[1].lower()
    if command in ["version", "--version", "-v"]:
        show_version()
        return
    elif command in ["help", "--help", "-h"]:
        show_help()
        return
    elif command == "input-conversion":
        scripts.input_conversion()
        return
    elif command == "ancestor-delta":
        scripts.run_get_cn_change_to_ancestor()
        return
    elif command == "ccd":
        scripts.run_calculate_ccd()
        return
    elif command == "run":
        run_alpaca()
    else:
        print(f"Unknown command: {command}")
        print("Run 'alpaca help' for available commands.")
        sys.exit(1)


def run_alpaca():
    # Configure logging
    logger = create_logger(name="ALPACA", log_dir="logs")
    logger.info("Starting ALPACA")
    config = make_config(sys.argv[1:])
    debug = config["preprocessing_config"]["debug"]
    # determine running mode:
    # if 'tumour', expect single file with all the segments and output a single file
    # if 'segment' expect array of files to segment files (can be from different tumours) and create separate outputs for each segment
    config, run_mode = set_run_mode(config)
    logger.info("-------------------------------------------------")
    logger.info("Running ALPACA with the following parameters:")
    # print value of each parameter:
    logger.info(config)
    print_logo()
    # initiate progress bar:
    if not debug:
        progress_bar = tqdm(
            total=len(config["preprocessing_config"]["input_files"]),
            desc="Processing files",
            unit="file",
            file=sys.stderr,
        )
        original_stdout = sys.stdout
        if os.name == "nt":  # Windows
            sys.stdout = open(os.devnull, "w")
        else:  # Unix/Linux
            sys.stdout = StringIO()
    try:
        for input_file_name in config["preprocessing_config"]["input_files"]:
            SS = SegmentSolution(input_file_name, config, logger)
            if (
                not config["preprocessing_config"]["overwrite_output"]
                and SS.output_exists()
            ):
                logger.warning(
                    f"Output for {input_file_name} already exists. Use '--overwrite_output 1' option to overwrite existing output. Skipping this segment."
                )
                continue
            SS.run_iterations()
            SS.find_optimal_solution()
            SS.get_solution()
            SS.save_output()
            if not debug:
                progress_bar.update(1)
                progress_bar.set_description(f"Processing {input_file_name}")
            else:
                logger.info(f"Segment {input_file_name} solved.")
        if run_mode == "tumour":
            concatenated_output_path = concatenate_output(
                config["preprocessing_config"]["output_directory"]
            )
            logger.info("Calculating copy number change to ancestor...")
            cn_change_to_ancestor_df = get_cn_change_to_ancestor(
                f"{SS.tumour_dir}/tree_paths.json", concatenated_output_path
            )
            save_dataframe_to_csv(
                df=cn_change_to_ancestor_df,
                output_dir=SS.config["preprocessing_config"]["output_directory"],
                output_filename="cn_change_to_ancestor.csv",
            )
            logger.info(
                f"""Analysis completed successfully. Output saved to: {SS.config["preprocessing_config"]["output_directory"]}"""
            )
        logger.info("Done")
    except Exception as e:
        logger.error(f"An error occurred: {e}")
        raise e
    finally:
        if not debug:
            sys.stdout.close()
            sys.stdout = original_stdout
            progress_bar.close()


if __name__ == "__main__":
    main()

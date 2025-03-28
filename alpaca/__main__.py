#!/usr/bin/env python3
import os
import sys
from tqdm import tqdm
from io import StringIO
from alpaca.ALPACA_segment_solution_class import SegmentSolution
from alpaca.utils import print_logo, concatenate_output, set_run_mode
from alpaca.make_configuration import make_config


def main():
    print("Starting ALPACA")
    config = make_config(sys.argv[1:])
    debug = config["preprocessing_config"]["debug"]
    # determine running mode:
    # if 'tumour', expect single file with all the segments and output a single file
    # if 'segment' expect array of files to segment files (can be from different tumours) and create separate outputs for each segment
    config, run_mode = set_run_mode(config)
    print("-------------------------------------------------")
    print("Running ALPACA with the following parameters:")
    # print value of each parameter:
    print(config)
    print_logo()
    # initiate progress bar:
    if not debug:
        progress_bar = tqdm(
            total=len(config["preprocessing_config"]["input_files"]),
            desc="Processing files",
            unit="file",
        )
        original_stdout = sys.stdout
        if os.name == "nt":  # Windows
            sys.stdout = open(os.devnull, "w")
        else:  # Unix/Linux
            sys.stdout = StringIO()
    try:
        for input_file_name in config["preprocessing_config"]["input_files"]:
            SS = SegmentSolution(input_file_name, config)
            SS.run_iterations()
            SS.find_optimal_solution()
            SS.get_solution()
            SS.save_output()
            if not debug:
                progress_bar.update(1)
                progress_bar.set_description(f"Processing {input_file_name}")
            else:
                print(f"Segment {input_file_name} solved.")
        if run_mode == "tumour":
            concatenate_output(config["preprocessing_config"]["output_directory"])
    finally:
        if not debug:
            sys.stdout.close()
            sys.stdout = original_stdout
            progress_bar.close()
            print("Done")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
import os
import sys

from tqdm import tqdm
from io import StringIO
from ALPACA_segment_solution_class import SegmentSolution
from utils import print_logo
from make_configuration import make_config


def main():
    print("Starting ALPACA")
    config = make_config(sys.argv[1:]) 
    print("-------------------------------------------------")
    print("Running ALPACA with the following parameters:")
    # print value of each parameter:
    print(config)
    debug = config["preprocessing_config"]["debug"]
    # exit if input_files is empty:
    if len(config['preprocessing_config']['input_files']) == 0:
        raise ValueError("No input files provided. Please specify at least one input file.")

    # initiate progress bar:
    if not config["preprocessing_config"]["debug"]:
        progress_bar = tqdm(total=len(config['preprocessing_config']['input_files']), desc="Processing files", unit="file")
        original_stdout = sys.stdout
        if os.name == "nt":  # Windows
            sys.stdout = open(os.devnull, "w")
        else:  # Unix/Linux
            sys.stdout = StringIO()
    try:
        for input_file_name in config['preprocessing_config']['input_files']:
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
    finally:
        if not debug:
            sys.stdout.close()
            sys.stdout = original_stdout
            progress_bar.close()
            print("Done")
            print_logo()


if __name__ == "__main__":
    main()
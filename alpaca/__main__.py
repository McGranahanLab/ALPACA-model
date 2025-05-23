#!/usr/bin/env python3
import os
import sys
from tqdm import tqdm
from io import StringIO
from alpaca.ALPACA_segment_solution_class import SegmentSolution
from alpaca.utils import print_logo, concatenate_output, set_run_mode
from alpaca.make_configuration import make_config
from datetime import datetime
import logging

def main():
    
    # Configure logging
    log_time = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_filename = f"run_log_{log_time}.log"
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        handlers=[
            logging.StreamHandler(),
            logging.FileHandler(log_filename)
        ])
    
    
    logging.info("Starting ALPACA")
    config = make_config(sys.argv[1:])
    debug = config["preprocessing_config"]["debug"]
    # determine running mode:
    # if 'tumour', expect single file with all the segments and output a single file
    # if 'segment' expect array of files to segment files (can be from different tumours) and create separate outputs for each segment
    config, run_mode = set_run_mode(config)
    logging.info("-------------------------------------------------")
    logging.info("Running ALPACA with the following parameters:")
    # print value of each parameter:
    logging.info(config)
    print_logo()
    # initiate progress bar:
    if not debug:
        progress_bar = tqdm(
            total=len(config["preprocessing_config"]["input_files"]),
            desc="Processing files",
            unit="file",
            file=sys.stderr 
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
        logging.info("Done")
    except Exception as e:
        logging.error(f"An error occurred: {e}")
        raise e
    finally:
        if not debug:
            print('-')
            sys.stdout.close()
            sys.stdout = original_stdout
            progress_bar.close()


if __name__ == "__main__":
    main()

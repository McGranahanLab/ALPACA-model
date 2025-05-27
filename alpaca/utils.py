import os
import pandas as pd
from pathlib import Path
import json
import importlib
import logging
from datetime import datetime
from typing import Optional
import logging


def show_version():
    try:
        version = importlib.metadata.version("alpaca")
        print(f"alpaca {version}")
    except importlib.metadata.PackageNotFoundError:
        print("alpaca version unknown (not installed)")


def show_help():
    print(
        "ALPACA = ALlele-specific Phylogenetic Analysis of clone Copy-number Alterations"
    )
    print_logo()
    print("")
    print("Usage:")
    print("  alpaca [command]")
    print("")
    print("Commands:")
    print("  version              Show version")
    print("  help                 Show this help")
    print("  run                  Run ALPACA")
    print("  input-conversion     Run input conversion")
    print("  ccd                  Calculate clone copy number diversity")
    print("")


def create_logger(name: str, log_dir: Optional[str] = None) -> logging.Logger:
    """
    Create a named logger with both console and file handlers.

    Args:
        name: Name for the logger
        log_dir: Optional directory for log files (defaults to current directory)

    Returns:
        Configured logger instance
    """
    # create logger
    logger = logging.getLogger(name)
    # check for active handlers
    if logger.handlers:
        return logger

    logger.setLevel(logging.INFO)
    formatter = logging.Formatter(
        fmt="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    # Console handler
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)

    # File handler
    log_time = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_filename = f"{name}_log_{log_time}.log"

    if log_dir:
        os.makedirs(log_dir, exist_ok=True)
        log_path = os.path.join(log_dir, log_filename)
    else:
        log_path = log_filename

    file_handler = logging.FileHandler(log_path)
    file_handler.setLevel(logging.INFO)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

    logger.info(f"Logger '{name}' initialized. Log file: {log_path}")

    return logger


def split_to_segments(tumour_dir: str) -> list[str]:
    segments_dir_path = f"{tumour_dir}/segments"
    df_path = f"{tumour_dir}/ALPACA_input_table.csv"
    os.makedirs(segments_dir_path, exist_ok=True)
    df = pd.read_csv(df_path)
    tumour_id = df["tumour_id"].iloc[0]
    assert len(
        df["tumour_id"].unique()
    ), "Found multiple tumour ids. In tummour mode only one tumour_id is allowed per input csv"
    segments = []
    for segment, segment_df in df.groupby("segment"):
        segment_df_path = (
            f"{segments_dir_path}/ALPACA_input_table_{tumour_id}_{segment}.csv"
        )
        segment_df.to_csv(segment_df_path, index=False)
        segments.append(segment_df_path)
    return segments


def concatenate_output(output_dir: str) -> str:
    logger = logging.getLogger("ALPACA")
    # keep only segment files in output files list
    output_files = [
        f
        for f in os.listdir(output_dir)
        if f.endswith(".csv") and (("optimal" in f) or ("all" in f))
    ]
    dfs = [pd.read_csv(f"{output_dir}/{f}") for f in output_files]
    concatenated_df = pd.concat(dfs)
    tumour_id = concatenated_df["tumour_id"].iloc[0]
    output_name = f"{output_dir}/ALPACA_output_{tumour_id}.csv"
    concatenated_df.to_csv(output_name, index=False)
    if os.path.exists(output_name):
        logger.info(f"Combined output saved to {output_name}")
        # remove segment files
        for f in output_files:
            if f != os.path.basename(output_name):
                os.remove(f"{output_dir}/{f}")
    else:
        logger.error(f"Failed to save combined output to {output_name}")
        raise FileNotFoundError(f"Output file not found: {output_name}")
    return output_name


def set_run_mode(config: dict) -> tuple[dict, str]:
    run_mode = config["preprocessing_config"]["mode"]
    if run_mode == "tumour":
        print("Running in tumour mode")
        # create segment files:
        config["preprocessing_config"]["input_files"] = [
            Path(x).name
            for x in split_to_segments(
                config["preprocessing_config"]["input_tumour_directory"]
            )
        ]
        config["preprocessing_config"]["input_data_directory"] = Path(
            config["preprocessing_config"]["input_tumour_directory"]
        ).parent
    return config, run_mode


def read_tree_json(json_path: str) -> list[list[str]]:
    with open(json_path, "r") as f:
        tree = json.load(f)
    return tree


def find_path_edges(branch, tree_edges):
    branch_edges = []
    for edge in tree_edges:
        if (edge[0] in branch) and (edge[1] in branch):
            branch_edges.append(edge)
    return set(branch_edges)


def get_tree_edges(tree):
    all_edges = list()
    for branch in tree:
        if len(branch) == 2:
            all_edges.append(tuple(branch))
        else:
            for i in range(len(branch) - 1):
                all_edges.append((branch[i], branch[i + 1]))
    unique_edges = set(all_edges)
    return unique_edges


def find_parent(tree, clone_name):
    for branch in tree:
        if branch[0] == clone_name:
            return "diploid"
        if clone_name in branch:
            clone_index = branch.index(clone_name)
            return branch[clone_index - 1]


def flat_list(target_list):
    if isinstance(target_list[0], list):
        return [item for sublist in target_list for item in sublist]
    else:
        return target_list


def get_length_from_name(segment):
    e = int(segment.split("_")[-1])
    s = int(segment.split("_")[-2])
    return e - s


def print_logo():
    print(
        """
     _____ __    _____ _____ _____ _____
    |  _  |  |  |  _  |  _  |     |  _  |
    |     |  |__|   __|     |   --|     |
    |__|__|_____|__|  |__|__|_____|__|__|
    /\\⌒⌒⌒/\\
    (⦿   ⦿)
    ( 'Y' )
     (   )
     (   )
     (   )
     (~ ~~~~~~~~~~)
     ( ~ ~~   ~~  )
     ( ~  ~ ~  ~  )
     (~  ~~~~~   ~)
     │ │      │ │
     │ │      │ │
    """
    )


def save_dataframe_to_csv(df: pd.DataFrame, output_dir: str, output_filename: str):
    """
    Save a DataFrame to a CSV file.

    Args:
        df: DataFrame to save
        output_dir: Path to directory where the CSV file will be saved
        output_filename: Name of the output CSV file
    """
    output_path = os.path.join(output_dir, output_filename)
    os.makedirs(os.path.dirname(output_dir), exist_ok=True)
    df.to_csv(output_path, index=False)

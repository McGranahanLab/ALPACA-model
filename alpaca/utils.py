import os
import pandas as pd
from pathlib import Path
import json


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


def concatenate_output(output_dir: str):
    output_files = [f for f in os.listdir(output_dir) if f.endswith(".csv")]
    dfs = [pd.read_csv(f"{output_dir}/{f}") for f in output_files]
    concatenated_df = pd.concat(dfs)
    tumour_id = concatenated_df["tumour_id"].iloc[0]
    concatenated_df.to_csv(f"{output_dir}/final_{tumour_id}.csv", index=False)
    # remove segment files
    for f in output_files:
        os.remove(f"{output_dir}/{f}")


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

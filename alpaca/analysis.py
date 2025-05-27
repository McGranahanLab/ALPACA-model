import pandas as pd
from .utils import find_parent, read_tree_json
from scipy.spatial import distance
import logging


### get_cn_change_to_ancestor ###
def get_parent_copynumbers(
    tree: list[list[str]], tumour_df: pd.DataFrame
) -> pd.DataFrame:
    assert (
        tumour_df["tumour_id"].nunique() == 1
    ), "Output should contain only one tumour_id"
    clone_parent_map = (
        tumour_df[["clone"]]
        .drop_duplicates()
        .apply(
            lambda x: pd.Series(
                {"clone": x["clone"], "parent": find_parent(tree, x["clone"])}
            ),
            axis=1,
        )
    )
    output_with_parent_clones = tumour_df.merge(clone_parent_map, on="clone")
    parents_df = tumour_df[["clone", "segment", "pred_CN_A", "pred_CN_B"]].copy()
    parents_df.rename(columns={"pred_CN_A": "parent_pred_cpnA"}, inplace=True)
    parents_df.rename(columns={"pred_CN_B": "parent_pred_cpnB"}, inplace=True)
    parents_df.rename(columns={"clone": "parent"}, inplace=True)
    output_with_parent_clones_copynumbers = output_with_parent_clones.merge(
        parents_df, left_on=["parent", "segment"], right_on=["parent", "segment"]
    )
    output_with_parent_clones_copynumbers["cn_dist_to_parent_A"] = (
        output_with_parent_clones_copynumbers["pred_CN_A"]
        - output_with_parent_clones_copynumbers["parent_pred_cpnA"]
    )
    output_with_parent_clones_copynumbers["cn_dist_to_parent_B"] = (
        output_with_parent_clones_copynumbers["pred_CN_B"]
        - output_with_parent_clones_copynumbers["parent_pred_cpnB"]
    )
    return output_with_parent_clones_copynumbers


def get_cn_change_to_ancestor(tree_path: str, tumour_df_path: str) -> pd.DataFrame:
    tumour_df = pd.read_csv(tumour_df_path)
    tree = read_tree_json(tree_path)
    return get_parent_copynumbers(tree, tumour_df)


### calculate_ccd ###
def calculate_ccd(results_path, metric="euclidean"):
    """
    example data:
        clone  pred_CN_A  pred_CN_B   segment
    0   clone1          1          3  1_10_100
    1  clone10          1          3  1_10_100
    2  clone12          1          2  1_10_100
    """
    logger = logging.getLogger("ccd")
    results_df = pd.read_csv(results_path)
    # validate input:
    # check if columns tumour_id, clone and segment are strings, while pred_CN_A and pred_CN_B are integers:
    if not all(
        results_df[col].dtype == "object" for col in ["tumour_id", "clone", "segment"]
    ):
        raise ValueError("tumour_id, clone and segment should be strings")
    if not all(results_df[col].dtype == "int" for col in ["pred_CN_A", "pred_CN_B"]):
        raise ValueError("pred_CN_A and pred_CN_B should be integers")
    # check for any NaN or empty values in the required columns:
    required_columns = ["tumour_id", "clone", "segment", "pred_CN_A", "pred_CN_B"]
    for col in required_columns:
        if results_df[col].isnull().any() or (results_df[col] == "").any():
            raise ValueError(f"Column {col} contains NaN or empty values")
    # check if the segment column contains the expected format:
    if not all(results_df["segment"].str.match(r"^\d+_\d+_\d+$")):
        raise ValueError(
            "segment column should contain the format 'chromosome_start_end'"
        )
    results_per_tumour = []
    tumour_ids = results_df["tumour_id"].unique()
    logger.info(f"Found {len(tumour_ids)} unique tumours")
    for tumour_id, tumour_df in results_df.groupby("tumour_id"):
        tumour_df["chromosome"] = (
            tumour_df["segment"].str.split("_", expand=True)[0].astype(int)
        )
        tumour_df["start"] = (
            tumour_df["segment"].str.split("_", expand=True)[1].astype(int)
        )
        distance_func = getattr(distance, metric)
        tumour_df = tumour_df.sort_values(["chromosome", "start", "clone"])
        vectors = [
            tumour_df[tumour_df.clone == clone][["pred_CN_A", "pred_CN_B"]]
            .unstack()
            .values
            for clone in tumour_df.clone.unique()
        ]
        vectors_names = [clone for clone in tumour_df.clone.unique()]
        max_difference = 0
        differences = []
        for i in range(len(vectors)):
            for j in range(i + 1, len(vectors)):
                diff = distance_func(vectors[i], vectors[j])
                differences.append(diff)
                if diff > max_difference:
                    max_difference = diff
                    most_different_vectors = (vectors_names[i], vectors_names[j])
        tumour_df["CCD"] = max_difference
        tumour_results_df = tumour_df[["tumour_id", "CCD"]].drop_duplicates().copy()
        results_per_tumour.append(tumour_results_df)
    ccd_df = pd.concat(results_per_tumour, ignore_index=True)
    return ccd_df

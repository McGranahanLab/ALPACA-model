import pandas as pd
from .utils import find_parent, read_tree_json, get_length_from_name, read_dataframe
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
    # add diploid clone so that delta for MRCA can also be calcualted
    diploid_frame = pd.DataFrame({
        "parent": ["diploid"] * tumour_df["segment"].nunique(),
        "segment": tumour_df["segment"].unique(),
        "parent_pred_cpnA": [1] * tumour_df["segment"].nunique(),
        "parent_pred_cpnB": [1] * tumour_df["segment"].nunique(),
    })
    parents_df = pd.concat([parents_df, diploid_frame], ignore_index=True)
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
    tumour_df = read_dataframe(tumour_df_path)
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
    results_df = read_dataframe(results_path)
    # validate input:
    # check if columns tumour_id, clone and segment are strings, while pred_CN_A and pred_CN_B are integers:
    if not all(results_df[col].dtype in ["str", "object"] for col in ["tumour_id", "clone", "segment"]):
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


### calculate_wgd_ratios ###
def calculate_wgd_ratios(results_path: str, tree_path: str) -> pd.DataFrame:
    """
    Calculate weighted copy-number increase ratios for each clone vs its parent.

    Returns a per-clone table with allele-specific and total ratios. For MRCA
    clones the parent is assumed to be diploid (CN=1,1).
    """
    logger = logging.getLogger("wgd")
    results_df = read_dataframe(results_path)
    required_columns = ["tumour_id", "clone", "segment", "pred_CN_A", "pred_CN_B"]
    missing = [col for col in required_columns if col not in results_df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")
    for col in required_columns:
        if results_df[col].isnull().any() or (results_df[col] == "").any():
            raise ValueError(f"Column {col} contains NaN or empty values")
    if not all(results_df["segment"].str.match(r"^\d+_\d+_\d+$")):
        raise ValueError(
            "segment column should contain the format 'chromosome_start_end'"
        )
    if results_df["tumour_id"].nunique() != 1:
        raise ValueError("calculate_wgd_ratios expects a single tumour input file")

    tree = read_tree_json(tree_path)

    df = results_df[required_columns].copy()
    # segment length for weighted averages
    df["segment_length"] = df["segment"].apply(get_length_from_name).astype(float)
    if (df["segment_length"] <= 0).any():
        raise ValueError("Found non-positive segment lengths in segment names")

    clone_parent_map = (
        df[["clone"]]
        .drop_duplicates()
        .apply(
            lambda x: pd.Series(
                {"clone": x["clone"], "parent": find_parent(tree, x["clone"])}
            ),
            axis=1,
        )
    )
    if clone_parent_map["parent"].isnull().any():
        missing_clones = clone_parent_map[clone_parent_map["parent"].isnull()][
            "clone"
        ].tolist()
        raise ValueError(f"Missing parent for clones: {missing_clones}")

    df = df.merge(clone_parent_map, on="clone")

    parents_df = df[["clone", "segment", "pred_CN_A", "pred_CN_B"]].copy()
    parents_df.rename(columns={"pred_CN_A": "parent_pred_cpnA"}, inplace=True)
    parents_df.rename(columns={"pred_CN_B": "parent_pred_cpnB"}, inplace=True)
    parents_df.rename(columns={"clone": "parent"}, inplace=True)

    diploid_frame = pd.DataFrame({
        "parent": ["diploid"] * df["segment"].nunique(),
        "segment": df["segment"].unique(),
        "parent_pred_cpnA": [1] * df["segment"].nunique(),
        "parent_pred_cpnB": [1] * df["segment"].nunique(),
    })
    parents_df = pd.concat([parents_df, diploid_frame], ignore_index=True)

    df = df.merge(parents_df, on=["parent", "segment"], how="left")
    if df[["parent_pred_cpnA", "parent_pred_cpnB"]].isnull().any().any():
        raise ValueError("Missing parent copy numbers for some clone/segment pairs")

    df["parent_pred_total"] = df["parent_pred_cpnA"] + df["parent_pred_cpnB"]
    df["child_pred_total"] = df["pred_CN_A"] + df["pred_CN_B"]

    def _weighted_ratio(group: pd.DataFrame, child_col: str, parent_col: str) -> tuple[float, int, float]:
        valid = group[parent_col] > 0
        if not valid.any():
            return float("nan"), 0, 0.0
        used = group[valid]
        ratios = used[child_col] / used[parent_col]
        weights = used["segment_length"]
        weighted_avg = (ratios * weights).sum() / weights.sum()
        return weighted_avg, int(valid.sum()), float(weights.sum())

    def _calc_group(group: pd.DataFrame) -> pd.Series:
        ratio_a, used_a, weight_a = _weighted_ratio(
            group, "pred_CN_A", "parent_pred_cpnA"
        )
        ratio_b, used_b, weight_b = _weighted_ratio(
            group, "pred_CN_B", "parent_pred_cpnB"
        )
        ratio_total, used_total, weight_total = _weighted_ratio(
            group, "child_pred_total", "parent_pred_total"
        )
        return pd.Series(
            {
                "tumour_id": group.name[0],
                "clone": group.name[1],
                "parent": group.name[2],
                "ratio_A": ratio_a,
                "ratio_B": ratio_b,
                "ratio_total": ratio_total,
                "segments_used_A": used_a,
                "segments_used_B": used_b,
                "segments_used_total": used_total,
                "weight_sum_A": weight_a,
                "weight_sum_B": weight_b,
                "weight_sum_total": weight_total,
            }
        )

    results = (
        df.groupby(["tumour_id", "clone", "parent"], sort=False)
        .apply(_calc_group)
        .reset_index(drop=True)
    )
    logger.info(f"Computed WGD ratios for {results['clone'].nunique()} clones")
    return results

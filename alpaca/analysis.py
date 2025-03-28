import pandas as pd
from .utils import find_parent, read_tree_json


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

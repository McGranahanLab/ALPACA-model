"""Builders for small, fully synthetic ALPACA input directories.

The mock tumour has a known ground truth (clone tree, clone proportions and
clone-level allele-specific copy numbers), so the observed sample-level copy
numbers are computed from it instead of being stored as opaque fixtures.
"""

import json
from pathlib import Path

import pandas as pd

TUMOUR_ID = "MOCK-Tumour1"
SAMPLE = "MOCK_S1"

# Root-to-leaf paths of the mock clone tree:
#   clone1 -> clone2 -> clone3
#   clone1 -> clone4
TREE_PATHS = [["clone1", "clone2", "clone3"], ["clone1", "clone4"]]

# Single-sample cancer cell fractions per clone; sums to 1.
CLONE_PROPORTIONS = {"clone1": 0.4, "clone2": 0.2, "clone3": 0.1, "clone4": 0.3}

# Ground-truth (A, B) copy numbers per clone for each mock segment.
SEGMENT_TRUTH = {
    "1_1000000_5000000": {
        "clone1": (2, 1),
        "clone2": (3, 1),
        "clone3": (4, 1),
        "clone4": (2, 0),
    },
    "2_1000000_5000000": {
        "clone1": (1, 1),
        "clone2": (1, 1),
        "clone3": (1, 0),
        "clone4": (2, 1),
    },
}

# Single-clone (monoclonal) variant: a tree with only the trunk clone.
SINGLE_CLONE_TREE_PATHS = [["clone1"]]
SINGLE_CLONE_PROPORTIONS = {"clone1": 1.0}
SINGLE_CLONE_TRUTH = {
    "1_1000000_5000000": {"clone1": (2, 1)},
    "2_1000000_5000000": {"clone1": (3, 0)},
}

CI_HALF_WIDTH = 0.2


def observed_copy_numbers(truth, proportions):
    """Proportion-weighted sample-level (cpnA, cpnB) implied by clone-level truth."""
    cpn_a = sum(proportions[c] * a for c, (a, _) in truth.items())
    cpn_b = sum(proportions[c] * b for c, (_, b) in truth.items())
    return cpn_a, cpn_b


def make_mock_tumour_dir(
    root,
    segments=None,
    samples=(SAMPLE,),
    tumour_id=TUMOUR_ID,
    tree_paths=None,
    proportions=None,
    truth=None,
):
    """Write a complete ALPACA tumour input directory under ``root``.

    Every sample gets the same clone proportions. Pass ``tree_paths``,
    ``proportions`` and ``truth`` to override the default 4-clone tumour. Returns the tumour directory.
    """
    tree_paths = tree_paths or TREE_PATHS
    proportions = proportions or CLONE_PROPORTIONS
    truth = truth or SEGMENT_TRUTH
    segments = list(segments or truth)
    tumour_dir = Path(root) / tumour_id
    seg_dir = tumour_dir / "segments"
    seg_dir.mkdir(parents=True, exist_ok=True)

    (tumour_dir / "tree_paths.json").write_text(json.dumps(tree_paths))

    cp = pd.DataFrame(
        {s: proportions for s in samples}
    ).rename_axis("clone").reset_index()
    cp.to_csv(tumour_dir / "cp_table.csv", index=False)

    rows, ci_rows = [], []
    for seg in segments:
        cpn_a, cpn_b = observed_copy_numbers(truth[seg], proportions)
        for s in samples:
            rows.append(
                {"tumour_id": tumour_id, "sample": s, "segment": seg, "cpnA": cpn_a, "cpnB": cpn_b}
            )
            ci_rows.append(
                {
                    "segment": seg,
                    "sample": s,
                    "lower_CI_A": max(cpn_a - CI_HALF_WIDTH, 0),
                    "upper_CI_A": cpn_a + CI_HALF_WIDTH,
                    "lower_CI_B": max(cpn_b - CI_HALF_WIDTH, 0),
                    "upper_CI_B": cpn_b + CI_HALF_WIDTH,
                    "tumour_id": tumour_id,
                    "ci_value": 0.5,
                }
            )
    table = pd.DataFrame(rows)
    table.to_csv(tumour_dir / "ALPACA_input_table.csv", index=False)
    pd.DataFrame(ci_rows).to_csv(tumour_dir / "ci_table.csv", index=False)
    for seg, group in table.groupby("segment"):
        group.to_csv(seg_dir / f"ALPACA_input_table_{tumour_id}_{seg}.csv", index=False)
    return tumour_dir

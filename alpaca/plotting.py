from pathlib import Path
import warnings
import argparse
import math
import json

import pandas as pd
import numpy as np
import seaborn as sns
import plotly.graph_objects as go
import networkx as nx
from plotly.subplots import make_subplots

from alpaca.utils import (
    read_tree_json,
    flat_list,
    ensure_chr_table,
    SUPPORTED_GENOME_BUILDS,
)
from alpaca.plotting_helpers import (
    get_chr_table,
    remove_duplicates_preserve_order,
    get_unique_lists,
    get_tree_edges,
    clean_output,
    find_parent,
)


_DEFAULT_GENOME_BUILD = "hg19"
_DEFAULT_HEATMAP_PALETTE = "classic"

_GAIN_PALETTE_MAP = {
    "classic": "Reds",
    "reds": "Reds",
    "reds_r": "Reds_r",
    "orrd": "OrRd",
    "ylorrd": "YlOrRd",
    "rocket": "rocket",
    "rocket_r": "rocket_r",
    "magma": "magma",
    "magma_r": "magma_r",
    "flare": "flare",
    "rdpu": "RdPu",
    "pink": "RdPu",
}
_SUPPORTED_HEATMAP_CHOICES = sorted(_GAIN_PALETTE_MAP.keys())

_LOSS_STATE_COLOUR = tuple(
    int(round(channel * 255)) for channel in sns.color_palette("Blues", 8)[-1]
)
_NEUTRAL_STATE_COLOUR = (255, 255, 255)


def _resolve_gain_palette_name(choice):
    normalized = (choice or _DEFAULT_HEATMAP_PALETTE).lower()
    palette = _GAIN_PALETTE_MAP.get(normalized)
    if palette is None:
        raise ValueError(
            "Unsupported heatmap palette '{choice}'. Choose one of: {options}.".format(
                choice=choice,
                options=", ".join(_SUPPORTED_HEATMAP_CHOICES),
            )
        )
    return palette


def _format_rgb(rgb_tuple):
    return f"rgb{tuple(int(channel) for channel in rgb_tuple)}"


def _normalize_chr_label(value):
    """Normalize chromosome labels to UCSC form used by chr tables.

    Handles ALPACA segment conventions where sex chromosomes may be encoded as
    23 (X) and 24 (Y), and accepts either raw tokens ("23") or UCSC-like
    labels ("chr23", "chrX").
    """

    if pd.isna(value):
        return np.nan

    token = str(value).strip()
    if not token:
        return np.nan

    if token.lower().startswith("chr"):
        token = token[3:]

    upper_token = token.upper()
    if upper_token in {"X", "Y"}:
        return f"chr{upper_token}"

    if token == "23":
        return "chrX"
    if token == "24":
        return "chrY"

    if token.isdigit():
        return f"chr{int(token)}"

    return f"chr{token}"


def build_copy_number_palette(max_state, palette_name=_DEFAULT_HEATMAP_PALETTE):
    """Generate a discrete palette mapping copy-number states to RGB tuples.

    Regardless of the palette selected, copy-number 0 is blue (loss),
    copy-number 1 stays white (diploid), and states >=2 are filled with
    progressively deeper shades drawn from a red-focused seaborn palette.
    """

    if max_state is None or math.isnan(max_state):
        max_state = 0
    max_state = int(math.ceil(max(0, max_state)))

    palette = {0: _LOSS_STATE_COLOUR}
    if max_state >= 1:
        palette[1] = _NEUTRAL_STATE_COLOUR

    red_states = max(0, max_state - 1)
    if not red_states:
        return palette

    gain_palette_name = _resolve_gain_palette_name(palette_name)
    reds = sns.color_palette(gain_palette_name, red_states)
    for idx, color in enumerate(reds, start=2):
        palette[idx] = tuple(int(round(channel * 255)) for channel in color)
    return palette


def load_chr_table(custom_path=None, genome_build=_DEFAULT_GENOME_BUILD):
    """Load chromosome lengths from a user path or the cached genome build."""
    if custom_path:
        return get_chr_table(custom_path)

    table_path = ensure_chr_table(genome_build)
    return get_chr_table(table_path)


def _read_table(table_path):
    sep = "\t" if str(table_path).lower().endswith((".tsv", ".txt")) else ","
    return pd.read_csv(table_path, sep=sep)


def load_mutation_table(input_dir, explicit_path=None):
    """Read the driver mutation table if it is present."""
    if explicit_path:
        path = Path(explicit_path)
        if not path.exists():
            raise FileNotFoundError(f"Mutation table not found: {path}")
        return _read_table(path)

    input_dir = Path(input_dir)
    candidates = (
        "driver_mutations.csv",
        "driver_mutations.tsv",
        "mutations.csv",
        "mutations.tsv",
    )
    for candidate in candidates:
        candidate_path = input_dir / candidate
        if candidate_path.exists():
            return _read_table(candidate_path)
    return None


def prepare_driver_mutations(mutation_df, tumour_id, chr_table):
    """Return driver mutations annotated with absolute coordinates for the selected tumour."""
    if mutation_df is None or mutation_df.empty:
        return None

    driver_df = mutation_df.copy()
    if "tumour_id" in driver_df.columns:
        driver_df = driver_df[driver_df["tumour_id"] == tumour_id]
        if driver_df.empty:
            return None

    if "clone" not in driver_df.columns:
        warnings.warn(
            "Driver mutation table is missing 'clone' column; skipping mutation overlay."
        )
        return None

    if "abs_position" not in driver_df.columns:
        chr_column = next(
            (c for c in ("chr", "chromosome", "chrom") if c in driver_df.columns), None
        )
        pos_column = next(
            (c for c in ("position", "pos", "start", "bp") if c in driver_df.columns),
            None,
        )
        if chr_column is None or pos_column is None:
            warnings.warn(
                "Driver mutation table requires chromosome and position columns to compute absolute coordinates; skipping mutation overlay."
            )
            return None
        chr_lookup = chr_table.set_index("chr")["shift"]
        driver_df["chr"] = (
            driver_df[chr_column]
            .astype(str)
            .apply(_normalize_chr_label)
        )
        driver_df["abs_position"] = driver_df["chr"].map(chr_lookup) + driver_df[
            pos_column
        ].astype(float)
        driver_df = driver_df.dropna(subset=["abs_position"])
        if driver_df.empty:
            return None

    return driver_df.reset_index(drop=True)


def find_alpaca_output_file(output_dir):
    matches = sorted(Path(output_dir).glob("ALPACA_output*.csv"))
    if not matches:
        raise FileNotFoundError(
            f"Could not find ALPACA output CSV in {output_dir}. Expected a file named 'ALPACA_output*.csv'."
        )
    return matches[0]


def plot_heatmap_with_tree(
    tree,
    alpaca_output,
    cp_table,
    chr_table,
    driver_mutations=None,
    allele="A",
    max_cpn_cap=8,
    heatmap_palette=_DEFAULT_HEATMAP_PALETTE,
):
    tumour_id = alpaca_output.tumour_id.iloc[0]
    mrca = tree[0][0]

    # drop diploid:
    alpaca_output = alpaca_output[alpaca_output.clone != "diploid"]
    # add coords:
    alpaca_output["chr"] = (
        alpaca_output.segment.str.split("_", expand=True)[0].apply(_normalize_chr_label)
    )
    alpaca_output["Start"] = alpaca_output.segment.str.split("_", expand=True)[
        1
    ].astype(int)
    alpaca_output["End"] = alpaca_output.segment.str.split("_", expand=True)[2].astype(
        int
    )
    # modify segment positions to absolute:
    alpaca_output = alpaca_output.merge(chr_table, on="chr", how="left")
    alpaca_output["abs_start"] = alpaca_output["Start"] + alpaca_output["shift"]
    alpaca_output["abs_end"] = alpaca_output["End"] + alpaca_output["shift"]
    alpaca_output = alpaca_output.sort_values(["abs_start"], ascending=False)
    alpaca_output.loc[alpaca_output["pred_CN_A"] > max_cpn_cap, "pred_CN_A"] = (
        max_cpn_cap
    )
    alpaca_output.loc[alpaca_output["pred_CN_B"] > max_cpn_cap, "pred_CN_B"] = (
        max_cpn_cap
    )
    max_cp_state = int(
        math.ceil(
            max(alpaca_output["pred_CN_A"].max(), alpaca_output["pred_CN_B"].max())
        )
    )
    copy_number_palette = build_copy_number_palette(
        max_cp_state, palette_name=heatmap_palette
    )

    number_of_clones = len(alpaca_output.clone.unique())
    max_levels = max([len(b) for b in tree])
    tree_with_levels = [dict(zip(b, range(0, len(b)))) for b in tree]
    tree_with_levels = (
        pd.concat(
            [
                pd.DataFrame(tree_with_levels[x], index=[0]).transpose()
                for x, _ in enumerate(tree_with_levels)
            ]
        )
        .reset_index()
        .drop_duplicates()
        .rename(columns={"index": "clone", 0: "level"})
    )

    # make empty y_loc df
    clone_y_location = dict(
        zip(alpaca_output.clone.unique(), range(0, number_of_clones))
    )
    clone_y_location = (
        pd.DataFrame(clone_y_location, index=[0])
        .transpose()
        .reset_index()
        .rename(columns={"index": "clone", 0: "y_loc"})
    )
    clone_y_location["y_loc"] = 100
    # find sections:
    # section is a part of a path that requires its own horizontal space on the final graph
    ori_tree = tree.copy()
    sections = [[tree[0][0]]]
    while len(flat_list(ori_tree)) > 1:
        for i, branch in enumerate(ori_tree):
            if (branch != []) and (
                remove_duplicates_preserve_order(flat_list(sections))
                != remove_duplicates_preserve_order(flat_list(tree))
            ):
                branching_clones = list(
                    pd.Series(flat_list(ori_tree))
                    .value_counts()[pd.Series(flat_list(ori_tree)).value_counts() > 1]
                    .index
                )
                if branching_clones == []:
                    branching_clones = [mrca]
                section_start = (
                    max([branch.index(x) for x in branching_clones if x in branch]) + 1
                )
                section = branch[section_start:]
                ori_tree[ori_tree.index(branch)] = branch[:section_start]
                ori_tree = get_unique_lists(ori_tree)
                if section != []:
                    sections.append(section)
    # order sections according to proximity
    # start with on of the longest paths
    section_termini = [x[-1] for x in sections]
    initial_node = [s for s in sections if len(s) == max([len(x) for x in sections])][
        0
    ][-1]
    # simplify tree graph:
    simple_tree = [
        [clone for clone in branch if clone in section_termini] for branch in tree
    ]

    edges = get_tree_edges(simple_tree)

    G = nx.Graph()
    for edge in edges:
        G.add_edge(edge[0], edge[1], weight=0)

    nodes = [initial_node] + [x for x in list(G.nodes) if x is not initial_node]
    distance_to_nodes = dict(nx.all_pairs_shortest_path_length(G))
    processed_nodes = [[s for s in sections if s[-1] == initial_node][0][-1]]

    while len(nodes) > 0:
        node = processed_nodes[-1]
        neighbours = pd.DataFrame(distance_to_nodes[node], index=["val"]).transpose()
        neighbours.drop(inplace=True, index=processed_nodes)
        try:
            neighbours.drop(inplace=True, index=node)
        except KeyError:
            pass
        neighbours = (
            neighbours.reset_index()
            .rename(columns={"index": "clone"})
            .merge(tree_with_levels, how="left", on="clone")
        )
        neighbours["level"] = neighbours["level"].astype(int)
        if len(neighbours) > 0:
            # if there are to neighbours equaly close, choose the one which has higher level (i.e. is deeper in the tree)
            # to do that, multiply proximity by negaive level:
            closest_candidates = neighbours[neighbours.val == min(neighbours.val)]
            closest_neighbours = closest_candidates[
                closest_candidates.level == max(closest_candidates.level)
            ]
            closest_neighbour = closest_neighbours.clone.values[0]
            if closest_neighbour not in processed_nodes:
                processed_nodes.append(closest_neighbour)
        nodes = [n for n in nodes if n != node]

    sorted_sections_raw = []
    for node in processed_nodes:
        sorted_sections_raw.append([s for s in sections if s[-1] == node][0])

    # join single-clone sections below MRCA with their descendants:
    sorted_sections = []
    below_MRCA = True
    skip_element = False
    for ss in sorted_sections_raw:
        if skip_element:
            skip_element = False
            continue
        if mrca in ss:
            below_MRCA = False
        if below_MRCA:
            if len(ss) == 1:
                if sorted_sections_raw.index(ss) < len(sorted_sections_raw) - 1:
                    joined = ss + sorted_sections_raw[sorted_sections_raw.index(ss) + 1]
                    sorted_sections.append(joined)
                else:
                    sorted_sections[-1] = sorted_sections[-1] + ss
                skip_element = True
            else:
                sorted_sections.append(ss)
        else:
            sorted_sections.append(ss)

    # assign y location on the plot to sorted sections:
    available_y_locs = range(0, number_of_clones)
    below_MRCA = True
    for section in sorted_sections:
        if mrca in section:
            below_MRCA = False
        locs_for_this_section = list(available_y_locs[: len(section)])
        if below_MRCA:
            locs_for_this_section = list(reversed(locs_for_this_section))
        locs_for_this_section_dict = dict(zip(section, locs_for_this_section))
        available_y_locs = available_y_locs[len(section) :]
        for n in locs_for_this_section_dict.keys():
            clone_y_location.loc[clone_y_location.clone == n, "y_loc"] = (
                locs_for_this_section_dict[n]
            )

    tree_graph_df = pd.merge(tree_with_levels, clone_y_location)

    total_plot_height = max(1000, 75 * number_of_clones)
    clone_prop_title = "Clone proportions in regions"

    s = [
        [
            {"type": "xy", "rowspan": number_of_clones},
            {"type": "xy", "rowspan": number_of_clones},
            {"type": "xy"},
        ]
    ]
    for c in range(number_of_clones - 1):
        s.append([None, None, {"type": "xy"}])

    fig = make_subplots(
        rows=number_of_clones,
        cols=3,
        column_widths=[0.1, 0.8, 0.1],
        specs=s,
        horizontal_spacing=0.02,
        vertical_spacing=0.01,
        subplot_titles=("", "", clone_prop_title),
    )

    for clone_pos in tree_graph_df.y_loc:
        hline = go.Scatter(
            showlegend=False,
            x=[-0.3, max_levels],
            y=[clone_pos + 0.5, clone_pos + 0.5],
            mode="lines",
            line=dict(color="Green", dash="dot"),
        )

        fig.append_trace(hline, row=1, col=1)

    # *** plot tree ***
    for branch in tree:
        branch_df = tree_graph_df[tree_graph_df.clone.isin(branch)]
        fig.append_trace(
            go.Scatter(
                showlegend=False,
                name="tree",
                x=branch_df["level"],
                y=branch_df["y_loc"],
                mode="lines+markers",
                marker=dict(
                    symbol="circle",
                    color="purple",
                    size=10,
                    line=dict(color="purple", width=2),
                ),
                text=branch_df["clone"],
            ),
            row=1,
            col=1,
        )
    fig.update_yaxes(
        showgrid=True,
        tickmode="array",
        tickvals=list(tree_graph_df.sort_values("y_loc").y_loc),
        ticktext=list(tree_graph_df.sort_values("y_loc").clone),
        range=[-0.5, number_of_clones - 0.5],
        showticklabels=True,
        zeroline=False,
        row=1,
        col=1,
    )
    fig.update_xaxes(showgrid=True, zeroline=True, row=1, col=1)
    fig.update_yaxes(
        showgrid=True,
        tickmode="array",
        tickvals=list(tree_graph_df.sort_values("y_loc").y_loc),
        ticktext=list(tree_graph_df.sort_values("y_loc").clone),
        range=[-0.5, number_of_clones - 0.5],
        showticklabels=True,
        zeroline=False,
        row=1,
        col=2,
    )
    fig.update_xaxes(showgrid=True, zeroline=True, row=1, col=2)

    # *** plot clones ***
    df = alpaca_output
    y_limit = df["pred_CN_A"].max()
    shapes = []
    chr_len = chr_table.copy()
    fig = plot_heat_map(
        alpaca_output.copy(),
        allele,
        fig,
        tree_graph_df,
        copy_number_palette,
        chr_table,
        driver_mutations,
    )
    for clone in df.clone.unique():
        i = int(tree_graph_df[tree_graph_df.clone == clone]["y_loc"].iloc[0])
        i = list(reversed(range(number_of_clones)))[i]

        fig.update_xaxes(showgrid=False, row=i + 1, col=2)
        # add proportions in regions:
        clone_cp = cp_table.loc[[clone]]

        showscale = i == len(df.clone.unique()) - 1
        clone_proportion_heatmap = go.Heatmap(
            z=clone_cp.values,
            x=clone_cp.columns,
            y=clone_cp.index,
            text=np.round(clone_cp.values, 2),
            texttemplate="%{text}",
            textfont={"size": 12},
            colorscale="Blues",
            showscale=False,
            colorbar=dict(
                tickfont=dict(size=12),
                orientation="h",
                x=0.89,
                y=-0.1,
                len=0.1,
                thickness=20,
            ),
            hoverinfo="z",
            zauto=False,
            zmin=0,
            zmax=1,
        )
        fig.add_trace(clone_proportion_heatmap, row=i + 1, col=3)

        # cleanup axes:

        fig.update_yaxes(showticklabels=False, row=i + 1, col=3)
        fig.update_yaxes(showticklabels=False, row=i + 1, col=2)
        sample_names = clone_cp.columns
        if i != len(df.clone.unique()) - 1:
            fig.update_xaxes(showticklabels=False, row=i + 1, col=1)

            fig.update_xaxes(showticklabels=False, row=i + 1, col=3)
        else:
            sample_names = clone_cp.columns
            # if sample names are in the long format, with tumour_id, split them:
            if tumour_id in sample_names[0]:
                sample_names = [x.split(f"{tumour_id}_")[1] for x in sample_names]
            fig.update_xaxes(
                tickmode="array",
                ticktext=sample_names,
                showticklabels=True,
                row=i + 1,
                col=3,
            )

        if i == 0:
            fig.update_xaxes(
                tickmode="array",
                tickvals=chr_len["cumsum"] - (chr_len["len"] / 2),
                ticktext=chr_len["chr"].astype(str).str.replace("chr", "", regex=False),
                showticklabels=True,
                row=i + 1,
                col=2,
            )
    # subtitle font size:
    fig.update_annotations(font_size=12)
    fig.update_layout(
        title=f"{tumour_id}<br>Allele: {allele}",
        plot_bgcolor="rgba(255,255,255,0)",
        autosize=False,
        width=1600,
        height=total_plot_height,
        legend_tracegroupgap=10,
        legend=dict(orientation="h", yanchor="top", y=1.4, xanchor="left", x=0.2),
    )

    # build legend:

    max_palette_state = max(copy_number_palette.keys())
    legend_items = [
        (
            str(state),
            _format_rgb(
                copy_number_palette.get(state, copy_number_palette[max_palette_state])
            ),
        )
        for state in range(0, max_palette_state + 1)
    ]
    for label, colour in legend_items:
        fig.add_trace(
            go.Scatter(
                legendgroup="copy-number",
                x=[None],
                y=[None],
                mode="markers",
                name=label,
                marker=dict(color=colour, size=10, line=dict(color="black", width=1)),
                showlegend=True,
            ),
            row=1,
            col=1,
        )

    return fig


def plot_cpn_per_clone(
    tree, alpaca_output, cp_table, chr_table, driver_mutations=None, max_cpn_cap=8
):
    tumour_id = alpaca_output.tumour_id.iloc[0]
    mrca = tree[0][0]

    # drop diploid:
    alpaca_output = alpaca_output[alpaca_output.clone != "diploid"]
    # add coords:
    alpaca_output["chr"] = (
        alpaca_output.segment.str.split("_", expand=True)[0].apply(_normalize_chr_label)
    )
    alpaca_output["Start"] = alpaca_output.segment.str.split("_", expand=True)[
        1
    ].astype(int)
    alpaca_output["End"] = alpaca_output.segment.str.split("_", expand=True)[2].astype(
        int
    )
    # modify segment positions to absolute:
    alpaca_output = alpaca_output.merge(chr_table, on="chr", how="left")
    alpaca_output["abs_start"] = alpaca_output["Start"] + alpaca_output["shift"]
    alpaca_output["abs_end"] = alpaca_output["End"] + alpaca_output["shift"]
    alpaca_output = alpaca_output.sort_values(["abs_start"], ascending=False)
    alpaca_output.loc[alpaca_output["pred_CN_A"] > max_cpn_cap, "pred_CN_A"] = (
        max_cpn_cap
    )
    alpaca_output.loc[alpaca_output["pred_CN_B"] > max_cpn_cap, "pred_CN_B"] = (
        max_cpn_cap
    )

    number_of_clones = len(alpaca_output.clone.unique())
    max_levels = max([len(b) for b in tree])
    tree_with_levels = [dict(zip(b, range(0, len(b)))) for b in tree]
    tree_with_levels = (
        pd.concat(
            [
                pd.DataFrame(tree_with_levels[x], index=[0]).transpose()
                for x, _ in enumerate(tree_with_levels)
            ]
        )
        .reset_index()
        .drop_duplicates()
        .rename(columns={"index": "clone", 0: "level"})
    )

    # make empty y_loc df
    clone_y_location = dict(
        zip(alpaca_output.clone.unique(), range(0, number_of_clones))
    )
    clone_y_location = (
        pd.DataFrame(clone_y_location, index=[0])
        .transpose()
        .reset_index()
        .rename(columns={"index": "clone", 0: "y_loc"})
    )
    clone_y_location["y_loc"] = 100
    # find sections:
    # section is a part of a path that requires its own horizontal space on the final graph
    ori_tree = tree.copy()
    sections = [[tree[0][0]]]
    while len(flat_list(ori_tree)) > 1:
        for i, branch in enumerate(ori_tree):
            if (branch != []) and (set(flat_list(sections)) != set(flat_list(tree))):
                branching_clones = list(
                    pd.Series(flat_list(ori_tree))
                    .value_counts()[pd.Series(flat_list(ori_tree)).value_counts() > 1]
                    .index
                )
                if branching_clones == []:
                    branching_clones = [mrca]
                section_start = (
                    max([branch.index(x) for x in branching_clones if x in branch]) + 1
                )
                section = branch[section_start:]
                ori_tree[ori_tree.index(branch)] = branch[:section_start]
                ori_tree = get_unique_lists(ori_tree)
                if section != []:
                    sections.append(section)
    # order sections according to proximity
    # start with on of the longest paths
    section_termini = [x[-1] for x in sections]
    initial_node = [s for s in sections if len(s) == max([len(x) for x in sections])][
        0
    ][-1]
    # simplify tree graph:
    simple_tree = [
        [clone for clone in branch if clone in section_termini] for branch in tree
    ]

    edges = get_tree_edges(simple_tree)

    G = nx.Graph()
    for edge in edges:
        G.add_edge(edge[0], edge[1], weight=0)

    nodes = [initial_node] + [x for x in list(G.nodes) if x is not initial_node]
    distance_to_nodes = dict(nx.all_pairs_shortest_path_length(G))
    processed_nodes = [[s for s in sections if s[-1] == initial_node][0][-1]]

    while len(nodes) > 0:
        node = processed_nodes[-1]
        neighbours = pd.DataFrame(distance_to_nodes[node], index=["val"]).transpose()
        neighbours.drop(inplace=True, index=processed_nodes)
        try:
            neighbours.drop(inplace=True, index=node)
        except KeyError:
            pass
        neighbours = (
            neighbours.reset_index()
            .rename(columns={"index": "clone"})
            .merge(tree_with_levels, how="left", on="clone")
        )
        neighbours["level"] = neighbours["level"].astype(int)
        if len(neighbours) > 0:
            # if there are to neighbours equaly close, choose the one which has higher level (i.e. is deeper in the tree)
            # to do that, multiply proximity by negaive level:
            closest_candidates = neighbours[neighbours.val == min(neighbours.val)]
            closest_neighbours = closest_candidates[
                closest_candidates.level == max(closest_candidates.level)
            ]
            closest_neighbour = closest_neighbours.clone.values[0]
            if closest_neighbour not in processed_nodes:
                processed_nodes.append(closest_neighbour)
        nodes = [n for n in nodes if n != node]

    sorted_sections_raw = []
    for node in processed_nodes:
        sorted_sections_raw.append([s for s in sections if s[-1] == node][0])

    # join single-clone sections below MRCA with their descendants:
    sorted_sections = []
    below_MRCA = True
    skip_element = False
    for ss in sorted_sections_raw:
        if skip_element:
            skip_element = False
            continue
        if mrca in ss:
            below_MRCA = False
        if below_MRCA:
            if len(ss) == 1:
                joined = ss + sorted_sections_raw[sorted_sections_raw.index(ss) + 1]
                sorted_sections.append(joined)
                skip_element = True
            else:
                sorted_sections.append(ss)
        else:
            sorted_sections.append(ss)

    # assign y location on the plot to sorted sections:
    available_y_locs = range(0, number_of_clones)
    below_MRCA = True
    for section in sorted_sections:
        if mrca in section:
            below_MRCA = False
        locs_for_this_section = list(available_y_locs[: len(section)])
        if below_MRCA:
            locs_for_this_section = list(reversed(locs_for_this_section))
        locs_for_this_section_dict = dict(zip(section, locs_for_this_section))
        available_y_locs = available_y_locs[len(section) :]
        for n in locs_for_this_section_dict.keys():
            clone_y_location.loc[clone_y_location.clone == n, "y_loc"] = (
                locs_for_this_section_dict[n]
            )

    tree_graph_df = pd.merge(tree_with_levels, clone_y_location)

    clone_prop_title = "Clone proportions in regions"

    s = [[{"type": "xy", "rowspan": number_of_clones}, {"type": "xy"}, {"type": "xy"}]]
    for c in range(number_of_clones - 1):
        s.append([None, {"type": "xy"}, {"type": "xy"}])

    fig = make_subplots(
        rows=number_of_clones,
        cols=3,
        column_widths=[0.1, 0.8, 0.1],
        specs=s,
        horizontal_spacing=0.02,
        vertical_spacing=0.01,
        subplot_titles=("", "", clone_prop_title),
    )

    for clone_pos in tree_graph_df.y_loc:
        hline = go.Scatter(
            showlegend=False,
            x=[-0.3, max_levels],
            y=[clone_pos + 0.5, clone_pos + 0.5],
            mode="lines",
            line=dict(color="Green", dash="dot"),
        )

        fig.append_trace(hline, row=1, col=1)

    # *** plot tree ***
    for branch in tree:
        branch_df = tree_graph_df[tree_graph_df.clone.isin(branch)]
        fig.append_trace(
            go.Scatter(
                showlegend=False,
                name="tree",
                x=branch_df["level"],
                y=branch_df["y_loc"],
                mode="lines+markers",
                marker=dict(
                    symbol="circle",
                    color="purple",
                    size=10,
                    line=dict(color="purple", width=2),
                ),
                text=branch_df["clone"],
            ),
            row=1,
            col=1,
        )
    fig.update_yaxes(
        showgrid=False,
        tickmode="array",
        tickvals=list(tree_graph_df.sort_values("y_loc").y_loc),
        ticktext=list(tree_graph_df.sort_values("y_loc").clone),
        range=[-0.3, number_of_clones - 0.7],
        showticklabels=True,
        zeroline=False,
        row=1,
        col=1,
    )
    fig.update_xaxes(showgrid=False, zeroline=False, row=1, col=1)
    total_plot_height = max(1000, 75 * number_of_clones)
    # *** plot clones ***
    df = alpaca_output
    shapes = []
    chr_len = chr_table.copy()
    for clone in df.clone.unique():
        i = int(tree_graph_df[tree_graph_df.clone == clone]["y_loc"].iloc[0])
        i = list(reversed(range(number_of_clones)))[i]
        y_limit = df["pred_CN_A"].max()
        clone_df = df[df.clone == clone]
        clone_df = clone_df[
            ["abs_start", "abs_end", "pred_CN_A", "pred_CN_B"]
        ].drop_duplicates()
        clone_df = clone_df.sort_values("abs_start", ascending=True)
        clone_df["space"] = None
        x = flat_list(
            [
                [v[1]["abs_start"], v[1]["abs_end"], v[1]["space"]]
                for v in clone_df.iterrows()
            ]
        )
        ya = flat_list(
            [
                [v[1]["pred_CN_A"], v[1]["pred_CN_A"], v[1]["space"]]
                for v in clone_df.iterrows()
            ]
        )
        yb = flat_list(
            [
                [v[1]["pred_CN_B"], v[1]["pred_CN_B"], v[1]["space"]]
                for v in clone_df.iterrows()
            ]
        )
        clone_df_sameAB = clone_df[clone_df["pred_CN_A"] == clone_df["pred_CN_B"]]
        if clone_df_sameAB.empty:
            yab = []
        else:
            yab = flat_list(
                [
                    [v[1]["pred_CN_B"], v[1]["pred_CN_B"], v[1]["space"]]
                    for v in clone_df_sameAB.iterrows()
                ]
            )
        clone_mutations = (
            driver_mutations[driver_mutations.clone == clone]
            if driver_mutations is not None
            else pd.DataFrame([])
        )

        fig.add_trace(
            go.Scatter(
                showlegend=False,
                x=x,
                y=ya,
                mode="lines",
                line=dict(color="rgb(255, 164, 0)", width=1),
            ),
            row=i + 1,
            col=2,
        )

        fig.add_trace(
            go.Scatter(
                showlegend=False,
                x=x,
                y=yb,
                mode="lines",
                line=dict(color="rgb(0, 128, 128)", width=1),
            ),
            row=i + 1,
            col=2,
        )

        fig.add_trace(
            go.Scatter(
                showlegend=False,
                x=x,
                y=yab,
                mode="lines",
                line=dict(color="rgb(187, 206, 186)", width=1),
            ),
            row=i + 1,
            col=2,
        )

        # add loss/gain areas:
        ancestor_clone = find_parent(clone, tree)
        if ancestor_clone == "diploid":
            ancestor_df = clone_df.copy()
            ancestor_df.loc[:, "pred_CN_A"] = 1
            ancestor_df.loc[:, "pred_CN_B"] = 1
        else:
            ancestor_df = df[df.clone == ancestor_clone]
            ancestor_df = ancestor_df[
                ["abs_start", "abs_end", "pred_CN_A", "pred_CN_B"]
            ].drop_duplicates()
            ancestor_df = ancestor_df.sort_values("abs_start", ascending=True)
            ancestor_df["space"] = None
        gain_loss_df = pd.merge(
            clone_df, ancestor_df, on=["abs_start", "abs_end", "space"]
        )
        gain_loss_df["gain"] = (
            gain_loss_df["pred_CN_A_x"] > gain_loss_df["pred_CN_A_y"]
        ) | (gain_loss_df["pred_CN_B_x"] > gain_loss_df["pred_CN_B_y"])
        gain_loss_df["loss"] = (
            gain_loss_df["pred_CN_A_x"] < gain_loss_df["pred_CN_A_y"]
        ) | (gain_loss_df["pred_CN_B_x"] < gain_loss_df["pred_CN_B_y"])

        gain_and_loss_df = gain_loss_df[
            (gain_loss_df["loss"] == True) & (gain_loss_df["gain"] == True)
        ]
        gain_df = gain_loss_df[
            (gain_loss_df["loss"] == False) & (gain_loss_df["gain"] == True)
        ]
        loss_df = gain_loss_df[
            (gain_loss_df["loss"] == True) & (gain_loss_df["gain"] == False)
        ]
        D = {
            "gain_and_loss_df": gain_and_loss_df,
            "gain_df": gain_df,
            "loss_df": loss_df,
        }

        xxyy = {}

        for d in D.items():
            if len(d[1] > 0):
                xxyy[f"x_{d[0]}"] = flat_list(
                    [
                        [
                            v[1]["abs_start"],
                            v[1]["abs_start"],
                            v[1]["abs_end"],
                            v[1]["abs_end"],
                            v[1]["abs_start"],
                            v[1]["space"],
                        ]
                        for v in d[1].iterrows()
                    ]
                )
                xxyy[f"y_{d[0]}"] = flat_list(
                    [
                        [0, y_limit, y_limit, 0, 0, v[1]["space"]]
                        for v in d[1].iterrows()
                    ]
                )

        try:
            fig.add_trace(
                go.Scatter(
                    showlegend=False,
                    x=xxyy["x_gain_df"],
                    y=xxyy["y_gain_df"],
                    mode="lines",
                    fill="toself",
                    line=dict(color="rgb(255, 230, 229,0.5)", width=1),
                ),
                row=i + 1,
                col=2,
            )
        except KeyError:
            pass
        try:
            fig.add_trace(
                go.Scatter(
                    showlegend=False,
                    x=xxyy["x_loss_df"],
                    y=xxyy["y_loss_df"],
                    mode="lines",
                    fill="toself",
                    line=dict(color="rgb(230, 229, 255,0.5)", width=1),
                ),
                row=i + 1,
                col=2,
            )
        except KeyError:
            pass
        try:
            fig.add_trace(
                go.Scatter(
                    showlegend=False,
                    x=xxyy["x_gain_and_loss_df"],
                    y=xxyy["y_gain_and_loss_df"],
                    mode="lines",
                    fill="toself",
                    line=dict(color="rgb(230, 207, 232,0.5)", width=1),
                ),
                row=i + 1,
                col=2,
            )
        except KeyError:
            pass
        if not clone_mutations.empty:
            fig.add_trace(
                go.Scatter(
                    x=clone_mutations.abs_position,
                    y=[0] * len(clone_mutations),
                    hovertext=clone_mutations.get("gene"),
                    mode="markers",
                    showlegend=False,
                    marker=dict(color="teal"),
                ),
                row=i + 1,
                col=2,
            )
        # add static annotations:
        """
        for index, row in clone_mutations.iterrows():
            fig.add_annotation(
                x=row['abs_position'],
                y=[0]*len(clone_mutations),
                text=row['gene'],
                font=dict(size=20),
                showarrow=False,
                yshift=0,
                xref=f"x{i+1}",
                yref=f"y2"
            )
         """
        fig.update_xaxes(showgrid=False, row=i + 1, col=2)

        # add proportions in regions:
        clone_cp = cp_table.loc[[clone]]

        showscale = i == len(df.clone.unique()) - 1
        clone_proportion_heatmap = go.Heatmap(
            z=clone_cp.values,
            x=clone_cp.columns,
            y=clone_cp.index,
            text=np.round(clone_cp.values, 2),
            texttemplate="%{text}",
            textfont={"size": 12},
            colorscale="Blues",
            showscale=False,
            colorbar=dict(
                tickfont=dict(size=12),
                orientation="h",
                x=0.89,
                y=-0.1,
                len=0.1,
                thickness=20,
            ),
            hoverinfo="z",
            zauto=False,
            zmin=0,
            zmax=1,
        )
        fig.add_trace(clone_proportion_heatmap, row=i + 1, col=3)

        # cleanup axes:
        for chromosome_line in chr_len["cumsum"]:
            fig.add_trace(
                go.Scatter(
                    x=[chromosome_line, chromosome_line],
                    y=[0, y_limit],
                    mode="lines",
                    line=dict(color="black", width=1, dash="dot"),
                    showlegend=False,
                ),
                row=i + 1,
                col=2,
            )
        fig.update_yaxes(showticklabels=False, row=i + 1, col=3)
        # fig.update_yaxes(showticklabels=False, row=i + 1, col=2)
        sample_names = clone_cp.columns
        if i != len(df.clone.unique()) - 1:
            fig.update_xaxes(showticklabels=False, row=i + 1, col=1)
            fig.update_xaxes(showticklabels=False, row=i + 1, col=2)
            fig.update_xaxes(showticklabels=False, row=i + 1, col=3)
        else:
            fig.update_xaxes(
                tickmode="array",
                tickvals=chr_len["cumsum"] - (chr_len["len"] / 2),
                ticktext=chr_len["chr"].astype(str).str.replace("chr", "", regex=False),
                showticklabels=True,
                row=i + 1,
                col=2,
            )
            sample_names = clone_cp.columns
            # if sample names are in the long format, with tumour_id, split them:
            if tumour_id in sample_names[0]:
                sample_names = [x.split(f"{tumour_id}_")[1] for x in sample_names]
            fig.update_xaxes(
                tickmode="array",
                ticktext=sample_names,
                showticklabels=True,
                row=i + 1,
                col=3,
            )

    # subtitle font size:
    fig.update_annotations(font_size=12)
    fig.update_layout(
        title=f"{tumour_id}",
        plot_bgcolor="rgb(255,255,255)",
        autosize=False,
        width=1600,
        height=total_plot_height,
        legend_tracegroupgap=10,
        legend=dict(orientation="h", yanchor="top", y=1.4, xanchor="left", x=0.2),
    )

    # build legend:
    # gain loss areas:
    fig.add_trace(
        go.Scatter(
            legendgroup="1",
            x=[None],
            y=[None],
            mode="markers",
            name="gain relative to parent",
            marker=dict(size=10, color="rgb(255, 230, 229)", symbol="square"),
        ),
        row=1,
        col=1,
    )
    fig.add_trace(
        go.Scatter(
            legendgroup="1",
            x=[None],
            y=[None],
            mode="markers",
            name="loss relative to parent",
            marker=dict(size=10, color="rgb(230, 229, 255)", symbol="square"),
        ),
        row=1,
        col=1,
    )
    fig.add_trace(
        go.Scatter(
            legendgroup="1",
            x=[None],
            y=[None],
            mode="markers",
            name="gain and loss relative to parent",
            marker=dict(size=10, color="rgb(230, 207, 232)", symbol="square"),
        ),
        row=1,
        col=1,
    )
    # cpn lines:
    fig.add_trace(
        go.Scatter(
            legendgroup="2",
            x=[None],
            y=[None],
            mode="markers+lines",
            name="A allele",
            marker=dict(size=10, color="rgb(255, 164, 0)", symbol="line-ew"),
        ),
        row=1,
        col=1,
    )
    fig.add_trace(
        go.Scatter(
            legendgroup="2",
            x=[None],
            y=[None],
            mode="markers+lines",
            name="B allele",
            marker=dict(size=10, color="rgb(0, 128, 128)", symbol="line-ew"),
        ),
        row=1,
        col=1,
    )
    fig.add_trace(
        go.Scatter(
            legendgroup="2",
            x=[None],
            y=[None],
            mode="markers+lines",
            name="A and B allele",
            marker=dict(size=10, color="rgb(187, 206, 186)", symbol="line-ew"),
        ),
        row=1,
        col=1,
    )

    return fig


def plot_sample_level_copy_numbers(
    sample_table,
    chr_table,
    allele="A",
    max_cpn_cap=8,
):
    """Plot sample-level copy numbers across the genome for a single allele."""

    if sample_table is None or sample_table.empty:
        raise ValueError("Sample-level copy number table is empty.")

    allele = (allele or "").upper()
    if allele not in {"A", "B"}:
        raise ValueError("Allele must be 'A' or 'B'.")

    cpn_col = f"cpn{allele}"
    required_columns = {"segment", "sample", cpn_col}
    missing = required_columns - set(sample_table.columns)
    if missing:
        raise ValueError(
            "Sample-level table is missing required columns: "
            + ", ".join(sorted(missing))
        )

    df = sample_table.copy()
    segment_parts = df["segment"].astype(str).str.split("_", expand=True)
    if segment_parts.shape[1] < 3:
        raise ValueError(
            "Segment column must be formatted as <chr>_<start>_<end>."
        )

    df["chr"] = segment_parts[0].apply(_normalize_chr_label)
    df["Start"] = pd.to_numeric(segment_parts[1], errors="coerce")
    df["End"] = pd.to_numeric(segment_parts[2], errors="coerce")
    df = df.dropna(subset=["Start", "End", "chr", cpn_col])
    if df.empty:
        raise ValueError("No valid segments found after parsing the sample table.")

    df["Start"] = df["Start"].astype(int)
    df["End"] = df["End"].astype(int)
    df[cpn_col] = pd.to_numeric(df[cpn_col], errors="coerce")
    df = df.dropna(subset=[cpn_col])
    if df.empty:
        raise ValueError("No valid copy-number values found in the sample table.")

    for col in ("cpnA", "cpnB"):
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce").round(2)

    df = df.merge(chr_table, on="chr", how="left")
    df = df.dropna(subset=["shift"])
    if df.empty:
        raise ValueError("Chromosome table did not match any segments in input.")

    df["abs_start"] = df["Start"] + df["shift"]
    df["abs_end"] = df["End"] + df["shift"]
    df = df.sort_values(["abs_start", "abs_end"])  # genome order

    if max_cpn_cap is not None:
        df.loc[df[cpn_col] > max_cpn_cap, cpn_col] = max_cpn_cap

    sample_order = remove_duplicates_preserve_order(df["sample"].astype(str).tolist())
    if not sample_order:
        raise ValueError("No samples found in the sample table.")

    tumour_id = df["tumour_id"].iloc[0] if "tumour_id" in df.columns else None
    sample_labels = sample_order
    if tumour_id:
        stripped = []
        for sample in sample_labels:
            prefix = f"{tumour_id}_"
            stripped.append(sample[len(prefix) :] if sample.startswith(prefix) else sample)
        sample_labels = stripped

    number_of_samples = len(sample_order)
    total_plot_height = max(600, 120 * number_of_samples)
    line_color = "rgb(255, 164, 0)" if allele == "A" else "rgb(0, 128, 128)"

    fig = make_subplots(
        rows=number_of_samples,
        cols=1,
        shared_xaxes=True,
        horizontal_spacing=0.01,
        vertical_spacing=0.08,
        subplot_titles=sample_labels,
    )

    genome_start = float(chr_table["shift"].min())
    genome_end = float(chr_table["cumsum"].max())

    for idx, sample in enumerate(sample_order, start=1):
        sample_df = df[df["sample"].astype(str) == sample]
        sample_df = sample_df[["abs_start", "abs_end", cpn_col]].drop_duplicates()
        sample_df = sample_df.sort_values("abs_start", ascending=True)
        sample_df["space"] = None
        epsilon = 0.03
        sample_df["plot_cpn"] = sample_df[cpn_col].where(
            sample_df[cpn_col] != 0, epsilon
        )
        sample_y_limit = float(sample_df[cpn_col].max())
        if math.isnan(sample_y_limit):
            sample_y_limit = 0.0
        sample_y_limit *= 1.1
        sample_y_limit = max(sample_y_limit, 1.0)

        x_vals = flat_list(
            [
                [row["abs_start"], row["abs_end"], row["space"]]
                for _, row in sample_df.iterrows()
            ]
        )
        y_vals = flat_list(
            [
                [row["plot_cpn"], row["plot_cpn"], row["space"]]
                for _, row in sample_df.iterrows()
            ]
        )

        fig.add_trace(
            go.Scatter(
                showlegend=False,
                x=x_vals,
                y=y_vals,
                mode="lines",
                line=dict(color=line_color, width=1),
            ),
            row=idx,
            col=1,
        )

        for chromosome_line in chr_table["cumsum"]:
            fig.add_trace(
                go.Scatter(
                    x=[chromosome_line, chromosome_line],
                    y=[0, sample_y_limit],
                    mode="lines",
                    line=dict(color="black", width=1, dash="dot"),
                    showlegend=False,
                ),
                row=idx,
                col=1,
            )

        fig.add_shape(
            type="rect",
            x0=genome_start,
            x1=genome_end,
            y0=0,
            y1=sample_y_limit,
            line=dict(color="black", width=1),
            fillcolor="rgba(0,0,0,0)",
            layer="below",
            row=idx,
            col=1,
        )

        fig.update_yaxes(range=[0, sample_y_limit], zeroline=False, row=idx, col=1)
        fig.update_xaxes(zeroline=False, row=idx, col=1)

        if idx != number_of_samples:
            fig.update_xaxes(showticklabels=False, row=idx, col=1)
        else:
            fig.update_xaxes(
                tickmode="array",
                tickvals=chr_table["cumsum"] - (chr_table["len"] / 2),
                ticktext=chr_table["chr"].str.replace("chr", "", regex=False),
                title_text="Genome",
                showticklabels=True,
                row=idx,
                col=1,
            )

    title_prefix = f"{tumour_id} " if tumour_id else ""
    fig.update_layout(
        title=f"{title_prefix}Sample-level copy numbers (Allele {allele})",
        plot_bgcolor="rgba(255,255,255,0)",
        autosize=False,
        width=1600,
        height=total_plot_height,
        showlegend=False,
    )
    fig.update_yaxes(title_text="Copy numbers", row=1, col=1)

    return fig


def plot_heat_map(
    patient_output,
    allele,
    fig,
    tree_graph_df,
    cpn_palette,
    chr_table,
    driver_mutations=None,
):
    clones = tree_graph_df.sort_values("y_loc", ascending=True).clone
    clone_number = len(clones)
    chromosome_table = chr_table.copy()
    tumour_id = patient_output.tumour_id.iloc[0]
    patient_output = clean_output(patient_output)
    patient_output["predicted_cpn"] = patient_output[f"pred_CN_{allele}"]
    # patient_output['fractional_cpn'] = patient_output[f'mphase{allele}_allele']
    patient_output = patient_output.merge(tree_graph_df)
    patient_output = patient_output.sort_values("y_loc")
    max_palette_state = max(cpn_palette.keys())

    def getColour(cp_state):
        try:
            cp_key = int(round(cp_state))
        except (TypeError, ValueError):
            cp_key = 0
        cp_key = max(0, min(cp_key, max_palette_state))
        colour = cpn_palette.get(cp_key, cpn_palette[max_palette_state])
        return _format_rgb(colour)

    # fig = go.Figure(layout_xaxis_range=[0, 3.2 * 10 ** 9], layout_yaxis_range=[0, len(clones) + 1])
    for clone_index, clone_name in enumerate(clones):
        clone_df = patient_output[patient_output["clone"] == clone_name]
        for cp_state in clone_df.predicted_cpn.unique():
            clone_df_cp_state = clone_df[clone_df.predicted_cpn == cp_state]
            segments_predicted = [
                [
                    tuple([row[1]["abs_start"], clone_index - 0.5]),
                    tuple([row[1]["abs_end"], clone_index + 0.5]),
                ]
                for row in clone_df_cp_state.iterrows()
            ]
            segments_predicted_unique = []
            for x in segments_predicted:
                if x not in segments_predicted_unique:
                    segments_predicted_unique.append(x)
            cpn_color = getColour(cp_state)
            for rectangle in segments_predicted_unique:
                segment = clone_df_cp_state[
                    clone_df_cp_state["abs_start"] == rectangle[0][0]
                ].segment.unique()[0]
                fig.add_trace(
                    go.Scatter(
                        showlegend=False,
                        x=[
                            rectangle[0][0],
                            rectangle[0][0],
                            rectangle[1][0],
                            rectangle[1][0],
                        ],
                        y=[
                            rectangle[0][1],
                            rectangle[1][1],
                            rectangle[1][1],
                            rectangle[0][1],
                        ],
                        # y=[rectangle[0][1], rectangle[0][1] + 1, rectangle[0][1] + 1, rectangle[0][1]],
                        fill="toself",
                        mode="lines",
                        fillcolor=cpn_color,
                        line_color=cpn_color,
                        name=f"clone: {rectangle[0][1]}, seg: {segment}",
                    ),
                    row=1,
                    col=2,
                )
        if driver_mutations is not None:
            clone_mutations = driver_mutations[driver_mutations.clone == clone_name]
            if not clone_mutations.empty:
                fig.add_trace(
                    go.Scatter(
                        x=clone_mutations.abs_position,
                        y=[clone_index] * len(clone_mutations),
                        mode="markers",
                        showlegend=False,
                        marker=dict(color="teal"),
                    ),
                    row=1,
                    col=2,
                )

    chromosomes = [
        [tuple([row[1]["cumsum"], 0]), tuple([row[1]["cumsum"], clone_number])]
        for row in chromosome_table.iterrows()
    ]
    for chromosome_line in chromosomes:
        fig.add_trace(
            go.Scatter(
                x=[chromosome_line[0][0], chromosome_line[0][0]],
                y=[-0.5, clone_number + 0.5],
                mode="lines",
                line=dict(color="black", width=1, dash="dash"),
                showlegend=False,
            ),
            row=1,
            col=2,
        )

    return fig


def _load_plot_inputs(
    tumour_input_dir,
    tumour_output_dir,
    alpaca_output_path=None,
    chr_table_override=None,
    mutation_table_override=None,
    genome_build=_DEFAULT_GENOME_BUILD,
):
    input_dir = Path(tumour_input_dir).expanduser().resolve()
    output_dir = Path(tumour_output_dir).expanduser().resolve()
    if alpaca_output_path:
        alpaca_path = Path(alpaca_output_path).expanduser().resolve()
    else:
        alpaca_path = find_alpaca_output_file(output_dir)

    if not alpaca_path.exists():
        raise FileNotFoundError(f"ALPACA output not found: {alpaca_path}")

    cp_table_path = input_dir / "cp_table.csv"
    if not cp_table_path.exists():
        raise FileNotFoundError(
            f"Missing clone proportion table at {cp_table_path}. Expected cp_table.csv in tumour input directory."
        )

    tree_json_path = input_dir / "tree_paths.json"
    tree_display_path = (
        tree_json_path
        if tree_json_path.exists()
        else tree_json_path.with_suffix(".nwk")
    )
    if not tree_display_path.exists():
        raise FileNotFoundError(
            f"Missing tumour tree. Expected {tree_json_path} or {tree_json_path.with_suffix('.nwk')}"
        )

    chr_table = load_chr_table(chr_table_override, genome_build=genome_build)
    tree = read_tree_json(str(tree_json_path))
    alpaca_output = pd.read_csv(alpaca_path)
    if alpaca_output.empty:
        raise ValueError("ALPACA output is empty; cannot produce plots.")
    cp_table = pd.read_csv(cp_table_path).set_index("clone")

    input_table_path = input_dir / "ALPACA_input_table.csv"
    sample_table = None
    if input_table_path.exists():
        sample_table = pd.read_csv(input_table_path)
    else:
        warnings.warn(
            f"Sample-level input table not found at {input_table_path}; skipping sample-level plots."
        )

    mutation_table = load_mutation_table(input_dir, mutation_table_override)
    tumour_id = alpaca_output.tumour_id.iloc[0]
    driver_mutations = prepare_driver_mutations(mutation_table, tumour_id, chr_table)

    return {
        "input_dir": input_dir,
        "output_dir": output_dir,
        "alpaca_output_path": alpaca_path,
        "cp_table_path": cp_table_path,
        "input_table_path": input_table_path,
        "tree_path_for_config": tree_display_path,
        "tree_json_path": tree_json_path,
        "tree": tree,
        "alpaca_output": alpaca_output,
        "cp_table": cp_table,
        "sample_table": sample_table,
        "chr_table": chr_table,
        "driver_mutations": driver_mutations,
        "tumour_id": tumour_id,
        "genome_build": genome_build,
    }


def _write_pdf_plots(plot_inputs, output_dir, heatmap_palette):
    tumour_id = plot_inputs["tumour_id"]
    output_dir = Path(output_dir).expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    generated = []
    for allele in ("A", "B"):
        fig = plot_heatmap_with_tree(
            tree=plot_inputs["tree"],
            alpaca_output=plot_inputs["alpaca_output"].copy(),
            cp_table=plot_inputs["cp_table"],
            chr_table=plot_inputs["chr_table"],
            driver_mutations=plot_inputs["driver_mutations"],
            allele=allele,
            heatmap_palette=heatmap_palette,
        )
        pdf_path = output_dir / f"{tumour_id}_{allele}_heatmap.pdf"
        fig.write_image(str(pdf_path))
        generated.append(pdf_path)

    cn_changes = plot_cpn_per_clone(
        tree=plot_inputs["tree"],
        alpaca_output=plot_inputs["alpaca_output"].copy(),
        cp_table=plot_inputs["cp_table"],
        chr_table=plot_inputs["chr_table"],
        driver_mutations=plot_inputs["driver_mutations"],
        max_cpn_cap=8,
    )
    cn_path = output_dir / f"{tumour_id}_cn_changes_per_clone.pdf"
    cn_changes.write_image(str(cn_path))
    generated.append(cn_path)

    sample_table = plot_inputs.get("sample_table")
    if sample_table is not None and not sample_table.empty:
        for allele in ("A", "B"):
            sample_fig = plot_sample_level_copy_numbers(
                sample_table=sample_table.copy(),
                chr_table=plot_inputs["chr_table"],
                allele=allele,
                max_cpn_cap=8,
            )
            sample_path = (
                output_dir / f"{tumour_id}_sample_copy_numbers_{allele}.pdf"
            )
            sample_fig.write_image(str(sample_path))
            generated.append(sample_path)
    else:
        warnings.warn(
            "Sample-level input table missing or empty; skipping sample-level plots."
        )

    return generated


def _format_source_block(text: str) -> list[str]:
    stripped = text.strip("\n")
    if not stripped:
        return ["\n"]
    lines = stripped.split("\n")
    return [f"{line}\n" for line in lines]


def _make_code_cell(source_text: str) -> dict:
    return {
        "cell_type": "code",
        "execution_count": None,
        "metadata": {"language": "python"},
        "outputs": [],
        "source": _format_source_block(source_text),
    }


def _build_plotting_notebook(plot_inputs, heatmap_palette):
    input_dir_literal = plot_inputs["input_dir"].as_posix()
    output_dir_literal = plot_inputs["output_dir"].as_posix()
    tree_path_literal = plot_inputs["tree_path_for_config"].as_posix()
    cp_table_literal = plot_inputs["cp_table_path"].as_posix()
    alpaca_output_literal = plot_inputs["alpaca_output_path"].as_posix()
    heatmap_palette_literal = json.dumps(
        heatmap_palette or _DEFAULT_HEATMAP_PALETTE
    )
    genome_build_literal = json.dumps(
        plot_inputs.get("genome_build", _DEFAULT_GENOME_BUILD)
    )

    imports_code = """from pathlib import Path

import pandas as pd
import plotly.io as pio

from alpaca.plotting import (
    load_chr_table,
    load_mutation_table,
    prepare_driver_mutations,
    plot_cpn_per_clone,
    plot_heatmap_with_tree,
    plot_sample_level_copy_numbers,
)
from alpaca.utils import read_tree_json

pio.renderers.default = "notebook_connected"
"""

    config_code = f"""INPUT_DIR = Path(r"{input_dir_literal}")
OUTPUT_DIR = Path(r"{output_dir_literal}")
TREE_PATH = Path(r"{tree_path_literal}")
CP_TABLE_PATH = Path(r"{cp_table_literal}")
ALPACA_OUTPUT_PATH = Path(r"{alpaca_output_literal}")
SAMPLE_TABLE_PATH = INPUT_DIR / "ALPACA_input_table.csv"

HEATMAP_PALETTE = {heatmap_palette_literal}
GENOME_BUILD = {genome_build_literal}

chr_table = load_chr_table(genome_build=GENOME_BUILD)
tree = read_tree_json(str(TREE_PATH))
alpaca_output = pd.read_csv(ALPACA_OUTPUT_PATH)
cp_table = pd.read_csv(CP_TABLE_PATH).set_index("clone")

sample_table = None
if SAMPLE_TABLE_PATH.exists():
    sample_table = pd.read_csv(SAMPLE_TABLE_PATH)
else:
    print("Sample-level input table not found; sample-level plots will be skipped.")

mutation_table = load_mutation_table(INPUT_DIR)
tumour_id = alpaca_output.tumour_id.iloc[0]
driver_mutations = prepare_driver_mutations(mutation_table, tumour_id, chr_table)


if driver_mutations is None:
    print("Driver mutation annotations unavailable; plots will omit mutation markers.")
else:
    print(f"Driver mutation annotations loaded")
"""

    heatmap_a_code = """heatmap_A = plot_heatmap_with_tree(
    tree=tree,
    alpaca_output=alpaca_output.copy(),
    cp_table=cp_table,
    chr_table=chr_table,
    driver_mutations=driver_mutations,
    allele="A",
    heatmap_palette=HEATMAP_PALETTE,
)
heatmap_A.show()
"""

    heatmap_b_code = """heatmap_B = plot_heatmap_with_tree(
    tree=tree,
    alpaca_output=alpaca_output.copy(),
    cp_table=cp_table,
    chr_table=chr_table,
    driver_mutations=driver_mutations,
    allele="B",
    heatmap_palette=HEATMAP_PALETTE,
)
heatmap_B.show()
"""

    cn_changes_code = """cn_changes = plot_cpn_per_clone(
    tree=tree,
    alpaca_output=alpaca_output.copy(),
    cp_table=cp_table,
    chr_table=chr_table,
    driver_mutations=driver_mutations,
    max_cpn_cap=8,
)
cn_changes.show()
"""

    sample_cpn_code = """if sample_table is None or sample_table.empty:
    print("Sample-level input table missing; skipping sample-level copy-number plots.")
else:
    sample_cpn_A = plot_sample_level_copy_numbers(
        sample_table=sample_table,
        chr_table=chr_table,
        allele="A",
        max_cpn_cap=8,
    )
    sample_cpn_A.show()

    sample_cpn_B = plot_sample_level_copy_numbers(
        sample_table=sample_table,
        chr_table=chr_table,
        allele="B",
        max_cpn_cap=8,
    )
    sample_cpn_B.show()
"""

    cells = [
        _make_code_cell(imports_code),
        _make_code_cell(config_code),
        _make_code_cell(heatmap_a_code),
        _make_code_cell(heatmap_b_code),
        _make_code_cell(cn_changes_code),
        _make_code_cell(sample_cpn_code),
    ]

    return {
        "cells": cells,
        "metadata": {
            "kernelspec": {
                "display_name": "Python 3",
                "language": "python",
                "name": "python3",
            },
            "language_info": {
                "name": "python",
                "pygments_lexer": "ipython3",
            },
        },
        "nbformat": 4,
        "nbformat_minor": 5,
    }


def _write_plot_notebook(plot_inputs, output_dir, notebook_name, heatmap_palette):
    notebook_path = Path(output_dir).expanduser().resolve() / notebook_name
    if notebook_path.suffix.lower() != ".ipynb":
        notebook_path = notebook_path.with_suffix(".ipynb")
    notebook_path.parent.mkdir(parents=True, exist_ok=True)
    notebook = _build_plotting_notebook(plot_inputs, heatmap_palette)
    with notebook_path.open("w", encoding="utf-8") as handle:
        json.dump(notebook, handle, indent=2)
    return notebook_path


def export_plot_outputs(
    mode,
    tumour_input_dir,
    tumour_output_dir,
    alpaca_output_path=None,
    chr_table_override=None,
    mutation_table_override=None,
    notebook_name=None,
    heatmap_palette=_DEFAULT_HEATMAP_PALETTE,
    genome_build=_DEFAULT_GENOME_BUILD,
):
    """Generate ALPACA visualisations as PDFs, notebooks, or skip entirely.

    The optional heatmap_palette argument allows callers to switch the
    copy-number heatmap colouring without re-implementing plotting logic. The
    genome_build parameter controls which UCSC chromosome lengths table to use
    when a custom --chr-table path is not provided.
    """

    normalized_mode = (mode or "notebook").lower()
    if normalized_mode not in {"pdf", "notebook", "none"}:
        raise ValueError(
            f"Unsupported plot output mode '{mode}'. Expected one of pdf, notebook, none."
        )
    if normalized_mode == "none":
        return []

    palette_choice = heatmap_palette or _DEFAULT_HEATMAP_PALETTE

    plot_inputs = _load_plot_inputs(
        tumour_input_dir=tumour_input_dir,
        tumour_output_dir=tumour_output_dir,
        alpaca_output_path=alpaca_output_path,
        chr_table_override=chr_table_override,
        mutation_table_override=mutation_table_override,
        genome_build=genome_build,
    )

    notebook_target = notebook_name or f"{plot_inputs['tumour_id']}_plots.ipynb"

    if normalized_mode == "notebook":
        notebook_path = _write_plot_notebook(
            plot_inputs=plot_inputs,
            output_dir=tumour_output_dir,
            notebook_name=notebook_target,
            heatmap_palette=palette_choice,
        )
        return [notebook_path]

    pdf_paths = _write_pdf_plots(plot_inputs, tumour_output_dir, palette_choice)
    return pdf_paths


def main():
    parser = argparse.ArgumentParser(description="Plot ALPACA heatmap with tree")
    parser.add_argument(
        "--output_directory", type=str, required=True, help="Output directory"
    )
    parser.add_argument(
        "--input_directory", type=str, required=True, help="Input directory"
    )
    parser.add_argument(
        "--chr-table",
        dest="chr_table",
        type=str,
        default=None,
        help="Optional path to a chromosome length table (defaults to the packaged hg19 table).",
    )
    parser.add_argument(
        "--mutation-table",
        dest="mutation_table",
        type=str,
        default=None,
        help="Optional path to a driver mutation table (CSV/TSV).",
    )
    parser.add_argument(
        "--plot_output_mode",
        dest="plot_output_mode",
        type=str,
        choices=["pdf", "notebook", "none"],
        default="notebook",
        help="Choose whether to save PDF figures, write an interactive notebook (default), or skip plotting.",
    )
    parser.add_argument(
        "--heatmap_palette",
        dest="heatmap_palette",
        type=str,
        default=_DEFAULT_HEATMAP_PALETTE,
        choices=_SUPPORTED_HEATMAP_CHOICES,
        help=(
            "Colour palette for copy-number gains (>=2). Use 'classic' for the legacy look or choose one of {options}.".format(
                options=", ".join(_SUPPORTED_HEATMAP_CHOICES)
            )
        ),
    )
    parser.add_argument(
        "--genome_build",
        dest="genome_build",
        type=str,
        default=_DEFAULT_GENOME_BUILD,
        choices=SUPPORTED_GENOME_BUILDS,
        help="Genome reference build used for chromosome lengths when plotting (default: hg19).",
    )
    args = parser.parse_args()
    output_directory = Path(args.output_directory).expanduser().resolve()
    input_directory = Path(args.input_directory).expanduser().resolve()

    alpaca_output_file = find_alpaca_output_file(output_directory)

    export_plot_outputs(
        mode=args.plot_output_mode,
        tumour_input_dir=input_directory,
        tumour_output_dir=output_directory,
        alpaca_output_path=alpaca_output_file,
        chr_table_override=args.chr_table,
        mutation_table_override=args.mutation_table,
        heatmap_palette=args.heatmap_palette,
        genome_build=args.genome_build,
    )


if __name__ == "__main__":
    main()

import pandas as pd

from pathlib import Path

from alpaca.plotting import _build_plotting_notebook, _normalize_chr_label, plot_segment_fit
from alpaca.plotting_helpers import get_chr_table


def test_normalize_chr_label_maps_sex_chromosomes():
    assert _normalize_chr_label("23") == "chrX"
    assert _normalize_chr_label("24") == "chrY"
    assert _normalize_chr_label("chr23") == "chrX"
    assert _normalize_chr_label("chr24") == "chrY"
    assert _normalize_chr_label("X") == "chrX"
    assert _normalize_chr_label("chrY") == "chrY"


def test_normalized_chr_labels_join_chr_table_shifts():
    chr_table = pd.DataFrame(
        {
            "chr": ["chr1", "chrX"],
            "shift": [0, 249250621],
        }
    )
    seg_tokens = pd.Series(["1", "23"]).apply(_normalize_chr_label)
    merged = pd.DataFrame({"chr": seg_tokens}).merge(chr_table, on="chr", how="left")
    assert merged["shift"].notna().all()


def test_get_chr_table_keeps_sex_chromosomes(tmp_path):
    table_path = tmp_path / "chr_lengths.csv"
    pd.DataFrame(
        {
            "chr": ["chr1", "chr2", "chrX", "chrY"],
            "len": [10, 20, 30, 40],
        }
    ).to_csv(table_path, index=False)

    chr_table = get_chr_table(table_path)
    assert chr_table["chr"].tolist() == ["chr1", "chr2", "chrX", "chrY"]
    assert chr_table.set_index("chr").loc["chrX", "shift"] == 30
    assert chr_table.set_index("chr").loc["chrY", "shift"] == 60


def test_plot_segment_fit_returns_expected_scores_and_traces():
    sample_table = pd.DataFrame(
        {
            "tumour_id": ["T1", "T1"],
            "sample": ["T1_S1", "T1_S2"],
            "segment": ["1_10_20", "1_10_20"],
            "cpnA": [1.2, 1.8],
            "cpnB": [0.8, 1.3],
        }
    )
    ci_table = pd.DataFrame(
        {
            "tumour_id": ["T1", "T1"],
            "sample": ["T1_S1", "T1_S2"],
            "segment": ["1_10_20", "1_10_20"],
            "lower_CI_A": [1.0, 1.7],
            "upper_CI_A": [1.4, 1.9],
            "lower_CI_B": [0.7, 1.0],
            "upper_CI_B": [0.9, 1.2],
        }
    )
    alpaca_output = pd.DataFrame(
        {
            "clone": ["clone1", "clone2"],
            "pred_CN_A": [1, 2],
            "pred_CN_B": [1, 1],
            "segment": ["1_10_20", "1_10_20"],
            "D_score": [0.6, 0.6],
            "CI_score": [1, 1],
        }
    )
    cp_table = pd.DataFrame(
        {
            "clone": ["clone1", "clone2"],
            "T1_S1": [0.8, 0.2],
            "T1_S2": [0.2, 0.8],
        }
    ).set_index("clone")

    fig = plot_segment_fit(
        sample_table=sample_table,
        ci_table=ci_table,
        alpaca_output=alpaca_output,
        cp_table=cp_table,
        segment="1_10_20",
    )

    assert len(fig.data) == 4
    assert fig.layout.meta["segment"] == "1_10_20"
    assert fig.layout.meta["D_score"] == 0.6
    assert fig.layout.meta["CI_score"] == 1
    assert list(fig.data[0].x) == ["S1", "S2"]
    assert list(fig.data[1].y) == [1.2, 1.8]
    assert list(fig.data[3].y) == [1.0, 1.0]


def test_plot_segment_fit_overlays_enlarged_ci_from_report():
    sample_table = pd.DataFrame(
        {
            "tumour_id": ["T1", "T1"],
            "sample": ["T1_S1", "T1_S2"],
            "segment": ["1_10_20", "1_10_20"],
            "cpnA": [1.2, 1.8],
            "cpnB": [0.8, 1.3],
        }
    )
    ci_table = pd.DataFrame(
        {
            "tumour_id": ["T1", "T1"],
            "sample": ["T1_S1", "T1_S2"],
            "segment": ["1_10_20", "1_10_20"],
            "lower_CI_A": [1.1, 1.7],
            "upper_CI_A": [1.3, 1.9],
            "lower_CI_B": [0.7, 1.0],
            "upper_CI_B": [0.9, 1.2],
        }
    )
    alpaca_output = pd.DataFrame(
        {
            "clone": ["clone1", "clone2"],
            "pred_CN_A": [1, 2],
            "pred_CN_B": [1, 1],
            "segment": ["1_10_20", "1_10_20"],
            "D_score": [0.6, 0.6],
            "CI_score": [1, 1],
        }
    )
    cp_table = pd.DataFrame(
        {
            "clone": ["clone1", "clone2"],
            "T1_S1": [0.8, 0.2],
            "T1_S2": [0.2, 0.8],
        }
    ).set_index("clone")
    ci_modified_report = pd.DataFrame(
        {
            "tumour_id": ["T1"],
            "segment": ["1_10_20"],
            "affected_sample": ["T1_S1"],
            "affected_allele": ["A"],
            "min_ci": [1.0],
        }
    )

    fig = plot_segment_fit(
        sample_table=sample_table,
        ci_table=ci_table,
        alpaca_output=alpaca_output,
        cp_table=cp_table,
        segment="1_10_20",
        ci_modified_report=ci_modified_report,
    )

    assert len(fig.data) == 5
    assert fig.data[1].name == "Enlarged CI A"
    assert fig.data[0].error_y.arrayminus[0] == 0.1
    assert fig.data[1].error_y.arrayminus[0] == 0.5


def test_build_plotting_notebook_appends_segment_fit_cell():
    notebook = _build_plotting_notebook(
        {
            "input_dir": Path("/tmp/input"),
            "output_dir": Path("/tmp/output"),
            "tree_path_for_config": Path("/tmp/input/tree_paths.json"),
            "cp_table_path": Path("/tmp/input/cp_table.csv"),
            "alpaca_output_path": Path("/tmp/output/ALPACA_output_T1.csv"),
            "tumour_id": "T1",
            "genome_build": "hg19",
        },
        "classic",
    )

    assert len(notebook["cells"]) == 7
    final_cell_source = "".join(notebook["cells"][-1]["source"])
    assert "TARGET_SEGMENT = None" in final_cell_source
    assert "plot_segment_fit(" in final_cell_source
    assert "ci_modified_report=ci_modified_report" in final_cell_source

    config_cell_source = "".join(notebook["cells"][1]["source"])
    assert "CI_MODIFIED_REPORT_PATH" in config_cell_source

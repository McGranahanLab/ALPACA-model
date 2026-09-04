import pandas as pd

from alpaca.plotting import _normalize_chr_label
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

"""Per-segment report tables are merged whether they were written as csv or parquet."""

import pandas as pd

from alpaca.utils import process_monoclonal_reports, process_run_summary_reports


def test_monoclonal_reports_merge_csv_and_parquet_and_delete_all(tmp_path):
    rows = {"tumour_id": ["T1"], "sample": ["S1"], "cpnA": [2.1], "cpnB": [1.0],
            "distance_to_integer_A": [0.1], "distance_to_integer_B": [0.0]}
    pd.DataFrame({**rows, "segment": ["1_1_10"]}).to_parquet(tmp_path / "T1_1_1_10_monoclonal_samples_report.parquet")
    pd.DataFrame({**rows, "segment": ["2_1_10"]}).to_parquet(tmp_path / "T1_2_1_10_monoclonal_samples_report.parquet")
    pd.DataFrame({**rows, "segment": ["3_1_10"]}).to_csv(tmp_path / "T1_3_1_10_monoclonal_samples_report.csv", index=False)
    out = tmp_path / "combined.csv"

    process_monoclonal_reports(str(tmp_path), delete=True, outpath=str(out))

    assert sorted(pd.read_csv(out)["segment"]) == ["1_1_10", "2_1_10", "3_1_10"]
    assert not list(tmp_path.glob("*_monoclonal_samples_report.*"))


def test_run_gap_summaries_merge_parquet(tmp_path):
    pd.DataFrame({"tumour_id": ["T1"], "segment": ["1_1_10"], "max_gap": [0.2], "gap_reason": ["time_limit"],
                  "runtime": [60], "optimal_complexity": [3], "strict_gap_enabled": [True]}).to_parquet(
        tmp_path / "T1_1_1_10_run_gap_summary.parquet"
    )
    out = tmp_path / "run_gap_summary.csv"
    process_run_summary_reports(str(tmp_path), outpath=str(out))
    assert pd.read_csv(out)["segment"].tolist() == ["1_1_10"]

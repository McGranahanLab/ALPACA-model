import os
import pandas as pd
import importlib

from alpaca.ALPACA_segment_solution_class import SegmentSolution
from alpaca.utils import read_dataframe


def make_minimal_instance(tmp_path):
    inst = object.__new__(SegmentSolution)
    inst.input_file_name = "ALPACA_input_table_T1_seg_name_x.csv"
    inst.tumour_id = "T1"
    inst.segment = "seg_name"
    inst.config = {"preprocessing_config": {"output_directory": str(tmp_path)}, "model_config": {}}
    inst.optimal_solution_index = 2
    inst.s_type = "s_strictly_decreasing"
    inst.solutions_combined = pd.DataFrame(
        {
            "complexity": [0, 1, 2, 3],
            "D_score": [5.0, 3.2, 1.5, 1.4],
            "CI_score": [0.1, 0.2, 0.15, 0.12],
            "allowed_complexity": [0, 1, 2, 3],
            "clone": ["diploid", "A", "B", "C"],
            "pred_CN_A": [1.0, 2.0, 3.0, 4.0],
            "pred_CN_B": [1.0, 2.0, 3.0, 4.0],
        }
    )
    inst.optimal_solution = inst.solutions_combined.query("allowed_complexity == 2").copy()
    return inst


def test_save_all_solutions_writes_csv(tmp_path):
    inst = make_minimal_instance(tmp_path)
    all_solutions = inst.solutions_combined[["clone", "pred_CN_A", "pred_CN_B", "complexity"]].copy()
    all_solutions["tumour_id"] = inst.tumour_id
    all_solutions["segment"] = inst.segment

    all_dir = tmp_path / "all_solutions"
    path = inst._save_all_solutions(str(all_dir), all_solutions)
    assert os.path.exists(path)
    df = read_dataframe(path)
    assert set(["clone", "pred_CN_A", "pred_CN_B", "complexity", "tumour_id", "segment"]).issubset(df.columns)
    assert len(df) == len(all_solutions)


def test_save_elbow_table_and_plot(tmp_path, monkeypatch):
    inst = make_minimal_instance(tmp_path)
    inst.elbow = {"s_min": 2, "s_raw": 2, "s_strictly_decreasing": 2, "raw_code": "default", "dec_code": "default"}
    inst.optimal_solution_index = 2

    all_dir = str(tmp_path / "all_solutions")
    os.makedirs(all_dir, exist_ok=True)

    elbow_path, elbow_df = inst._save_elbow_table(all_dir)
    assert elbow_path is not None
    assert os.path.exists(elbow_path)
    assert "knee_s_min" in elbow_df.columns
    assert elbow_df["knee_s_min"].iloc[0] == inst.elbow["s_min"]

    # fake plotting backend
    class FakeFig:
        def savefig(self, path, bbox_inches=None):
            open(path, "wb").close()

    class FakeAx:
        def plot(self, *args, **kwargs):
            pass

        def axvline(self, *args, **kwargs):
            pass

        def set_xlabel(self, *a, **k):
            pass

        def set_ylabel(self, *a, **k):
            pass

        def set_title(self, *a, **k):
            pass

        def legend(self):
            pass

    class FakePlt:
        def subplots(self):
            return FakeFig(), FakeAx()

    mod = importlib.import_module("alpaca.ALPACA_segment_solution_class")
    monkeypatch.setattr(mod, "plt", FakePlt())

    plot_path = inst._plot_elbow(all_dir, elbow_df)
    # If plotting returns a path it should exist; if None, plotting was skipped
    if plot_path:
        assert os.path.exists(plot_path)

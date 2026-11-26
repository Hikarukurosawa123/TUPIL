"""Rerun the original Final_100 models with 5 folds and 10 repeats."""

import argparse
import runpy
import sys
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
import pandas as pd
from scipy.stats import norm
from sklearn.metrics import roc_auc_score

REPO = Path(__file__).resolve().parents[2]
RESULTS = REPO / "Results_Rough" / "Final_100"
OUTPUT = RESULTS / "auc_5fold_10repeats"
SOURCE = (
    REPO.parent / "POCUS_Kidney_Analysis" / "Analysis_scripts" / "Redo_analysis_1106_2026"
    / "ci_biopsy_spatial_analysis_cv_nested" / "redo_analysis_clean_1122_2026"
    / "kp_cases" / "run_spatial_94_simple.py"
)
FOLDS, REPEATS = 5, 10
CV_LABEL = f"{FOLDS} folds × {REPEATS} repeats"


def export_results(workbook):
    with pd.ExcelFile(workbook) as source:
        auc = pd.read_excel(source, "AUC")
        scores = pd.read_excel(source, "Patient scores")
        paired = pd.read_excel(source, "ROI comparisons DeLong")
        timing = pd.read_excel(source, "Early vs late DeLong")

    registered = pd.read_excel(RESULTS / "spatial_Final_100.xlsx", sheet_name="Patient scores")
    pd.testing.assert_frame_equal(scores.iloc[:, :6], registered.iloc[:, :6])
    assert len(scores) == 100 and scores.scan_id.is_unique
    for row in auc.itertuples(index=False):
        cases = scores if row.half == "whole" else scores[scores.half.eq(row.half)]
        prefix = "frozen" if row.model == "frozen" else f"{row.model}_withinhalf"
        truth = cases.ci_grade.gt(0)
        assert (len(cases), truth.sum()) == (row.n, row.n_ci12)
        np.testing.assert_allclose(roc_auc_score(truth, cases[f"{prefix}__{row.spatial}"]), row.AUC)

    OUTPUT.mkdir(parents=True, exist_ok=True)
    for name, frame in {
        "AUC": auc,
        "ROI_comparisons": paired,
        "early_vs_late": timing,
    }.items():
        if "z" in frame:
            np.testing.assert_allclose(frame.p_DeLong_1sided, norm.sf(frame.z))
            frame = frame.drop(columns="p_DeLong_2sided").assign(multiple_comparison_correction="none")
        frame.assign(cv_folds=FOLDS, cv_repeats=REPEATS).to_csv(OUTPUT / f"{name}.csv", index=False)

    plot = runpy.run_path(str(REPO / "Scripts" / "final_100_auc_by_roi" / "plot_auc.py"))
    plot["main"](workbook=workbook, output_dir=OUTPUT, cv_label=CV_LABEL, include_binary=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fit-source", type=Path, default=SOURCE)
    source = parser.parse_args().fit_source.resolve()
    with TemporaryDirectory(prefix="tupil_5x10_") as work_dir:
        work = Path(work_dir)
        sys.path.insert(0, str(source.parent))
        args = sys.argv
        try:
            sys.argv = [str(source), "--cohort101", "--exclude", "291", "--plain", "--out", str(work)]
            analysis = runpy.run_path(str(source))
            analysis["S"].FOLDS = FOLDS
            analysis["S"].REPEATS = REPEATS
            analysis["main"]()
        finally:
            sys.argv = args
            sys.path.pop(0)
        export_results(work / f"spatial_{work.name}.xlsx")


if __name__ == "__main__":
    main()

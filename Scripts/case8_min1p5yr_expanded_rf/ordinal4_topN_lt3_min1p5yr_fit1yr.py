"""Ordinal 4-category top-N radiomics models (LOOCV only), eGFR slope < -3, >= 1.5 yr expanded-RF cohort, kidney-ROI
patients (n = 130), with the slope taken from the LMM fitted on the >= 1 yr cohort (lmm_fit1yr_restrict1p5yr.py ->
case8_min1p5yr_expanded_rf_slope_groups_fit1yr.xlsx).

Pipeline / features / scoring identical to ordinal4_topN_lt3_min1p5yr.py (cv_scores, boot_ci, c_index imported):
  Sets   SS, ESD, H-scan alone; ESD + SS, SS + H-scan, ESD + SS + H-scan with top-N chosen overall (pooled)
  Top N  1, 2, 3
  CV     LOOCV only; AUC for slope < -3 with 95 % CI (2000 stratified bootstrap resamples), 4-cat C-index.
Original results (slope fitted on >= 1.5 yr) are read from ordinal4_topN_lt3_min1p5yr.xlsx for comparison.
Run: C:/Users/tupil/venvs/radiomics_env/python.exe ordinal4_topN_lt3_min1p5yr_fit1yr.py
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import LeaveOneOut

sys.path.insert(0, str(Path(__file__).resolve().parent))
import ordinal4_topN_lt3_min1p5yr as O  # noqa: E402

RES = O.RES
NEW = RES / "case8_min1p5yr_expanded_rf_slope_groups_fit1yr.xlsx"
SETS = [(["SS"], "single"), (["ESD"], "single"), (["HScan"], "single"),
        (["ESD", "SS"], "overall"), (["SS", "HScan"], "overall"), (["ESD", "SS", "HScan"], "overall")]
TOPN = (1, 2, 3)


def run(X, slope, mods, mode, topn):
    names = [c for c in X.columns if c.split("__", 1)[0] in mods]
    Xm = X[names].to_numpy(float)
    y4 = sum((slope < c).astype(int) for c in O.CUTS4)
    y3 = (slope < -3).astype(int)
    s, chosen = O.cv_scores(Xm, y4, names, topn, True, LeaveOneOut().split(Xm))
    lo, hi = O.boot_ci(y3, s)
    freq = pd.Series(chosen).value_counts() / len(y4)
    return dict(feature_set=O.set_name(mods, mode), top_n=topn, n_candidates=len(names),
                LOOCV_AUC_lt3=roc_auc_score(y3, s), CI_low=lo, CI_high=hi, LOOCV_Cindex_4cat=O.c_index(y4, s),
                selected=", ".join(f"{f} ({v:.2f})" for f, v in freq.head(max(3, topn + 1)).items())), s


if __name__ == "__main__":
    X, slope_old = O.load()
    slope = pd.read_excel(NEW, "Cohort").set_index("Patient Case ID")[O.S].reindex(X.index).to_numpy(float)
    assert not np.isnan(slope).any()
    y3 = (slope < -3).astype(int)
    print(f"n = {len(y3)}; < -3: {y3.sum()} (was {(slope_old < -3).sum()})")
    jobs = [(m, md, n) for m, md in SETS for n in TOPN]
    out = Parallel(n_jobs=len(jobs))(delayed(run)(X, slope, m, md, n) for m, md, n in jobs)
    res = pd.DataFrame([o[0] for o in out])
    old = pd.read_excel(RES / "ordinal4_topN_lt3_min1p5yr.xlsx", "Performance")
    old = old[["feature_set", "top_n", "LOOCV_AUC_lt3", "CI_low", "CI_high"]].rename(
        columns={"LOOCV_AUC_lt3": "old_LOOCV_AUC", "CI_low": "old_CI_low", "CI_high": "old_CI_high"})
    both = res.merge(old, on=["feature_set", "top_n"], how="left")
    both["delta (new - old)"] = both.LOOCV_AUC_lt3 - both.old_LOOCV_AUC
    scores = pd.DataFrame({f"{O.set_name(m, md)} top{n}": o[1] for (m, md, n), o in zip(jobs, out)}, index=X.index)

    with pd.ExcelWriter(RES / "ordinal4_topN_lt3_min1p5yr_fit1yr.xlsx") as w:
        pd.DataFrame({"Item": ["Method"], "Value": [__doc__]}).to_excel(w, sheet_name="Read_Me", index=False)
        both.to_excel(w, sheet_name="Performance", index=False)
        scores.assign(slope_fit1yr=slope, slope_fit1p5yr=slope_old, lt_minus3=y3).reset_index() \
            .to_excel(w, sheet_name="LOOCV scores", index=False)
    pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 160)
    print(both.drop(columns="selected").round(3).to_string(index=False))
    print(both[["feature_set", "top_n", "selected"]].to_string(index=False))

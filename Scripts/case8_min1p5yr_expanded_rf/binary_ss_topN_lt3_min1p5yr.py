"""SS top-N (1, 2, 3) radiomics -> BINARY logistic regression for eGFR slope < -3, >= 1.5 yr expanded-RF cohort,
kidney-ROI patients (n = 130). Same features / CV as ordinal4_topN_lt3_min1p5yr.py; model as the 01 'Binary < -3' task
(run_model_grid.py).

Inside each training fold: median imputation -> Spearman-priority |Pearson r| >= 0.9 removal -> top-N SS features by
|Spearman rho| -> standardise -> LogisticRegression(penalty=None). Held-out score = x.beta without the intercept.
Ranking target (both the correlation-removal priority and the top-N ranking):
  'binary'  Spearman with the binary outcome (slope < -3)        <- the 01 binary task
  '4-cat'   Spearman with the 4-category slope group (as the ordinal models), model still trained on < -3
Reported: LOOCV AUC (95 % CI, 2000 stratified bootstrap resamples) and stratified 5-fold x 10 mean AUC.
Run: C:/Users/tupil/venvs/radiomics_env/python.exe binary_ss_topN_lt3_min1p5yr.py
"""
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import LeaveOneOut, StratifiedKFold

sys.path.insert(0, str(Path(__file__).resolve().parent))
import ordinal4_topN_lt3_min1p5yr as O  # noqa: E402


def scores(X, y3, yrank, names, topn, splits):
    sys.path.insert(0, str(O.MODEL_CODE))
    warnings.filterwarnings("ignore")
    from nested_loocv_radiomics import Preprocessor, select_top_spearman_per_modality
    s = np.full(len(y3), np.nan); chosen = []
    for tr, te in splits:
        prep, xtr = Preprocessor.fit(X[tr], names, 0.90, yrank[tr], "spearman")
        xte = prep.transform(X[te])
        idx, _ = select_top_spearman_per_modality(xtr, yrank[tr], prep.selected_features, topn, True)
        m = LogisticRegression(penalty=None, max_iter=10000).fit(xtr[:, idx], y3[tr])
        s[te] = xte[:, idx] @ m.coef_[0]
        chosen += [prep.selected_features[i] for i in idx]
    return s, chosen


def run(X, y3, y4, names, rank, topn):
    yrank = y3 if rank == "binary" else y4
    s, chosen = scores(X, y3, yrank, names, topn, LeaveOneOut().split(X))
    auc = roc_auc_score(y3, s)
    rng = np.random.default_rng(O.SEED)
    i1, i0 = np.where(y3 == 1)[0], np.where(y3 == 0)[0]
    bs = [roc_auc_score(y3[b], s[b]) for b in (np.r_[rng.choice(i1, len(i1)), rng.choice(i0, len(i0))] for _ in range(O.NBOOT))]
    kf = []
    for r in range(O.REPS):
        for tr, te in StratifiedKFold(5, shuffle=True, random_state=O.SEED + r).split(X, y3):
            sk, _ = scores(X, y3, yrank, names, topn, [(tr, te)])
            kf.append(roc_auc_score(y3[te], sk[te]))
    freq = pd.Series(chosen).value_counts() / len(y3)
    return dict(ranking=rank, top_n=topn, LOOCV_AUC_lt3=auc, CI_low=np.percentile(bs, 2.5), CI_high=np.percentile(bs, 97.5),
                kfold5x10_mean_AUC_lt3=np.mean(kf), kfold5x10_AUC_IQR=f"{np.percentile(kf, 25):.2f}–{np.percentile(kf, 75):.2f}",
                selected=", ".join(f"{f.split('__', 1)[1]} ({v:.2f})" for f, v in freq.head(topn + 2).items())), s


if __name__ == "__main__":
    X, slope = O.load()
    y3 = (slope < -3).astype(int)
    y4 = sum((slope < c).astype(int) for c in O.CUTS4)
    names = [c for c in X.columns if c.startswith("SS__")]
    Xs = X[names].to_numpy(float)
    jobs = [(r, n) for r in ("binary", "4-cat") for n in (1, 2, 3)]
    out = Parallel(n_jobs=6)(delayed(run)(Xs, y3, y4, names, r, n) for r, n in jobs)
    res = pd.DataFrame([o[0] for o in out])
    with pd.ExcelWriter(O.RES / "binary_ss_topN_lt3_min1p5yr.xlsx") as w:
        pd.DataFrame({"Item": ["Method"], "Value": [__doc__]}).to_excel(w, sheet_name="Read_Me", index=False)
        res.to_excel(w, sheet_name="Performance", index=False)
        pd.DataFrame({f"rank {r} top{n}": o[1] for (r, n), o in zip(jobs, out)}, index=X.index) \
            .assign(slope=slope, lt_minus3=y3).to_excel(w, sheet_name="LOOCV scores")
    pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 200)
    print(f"n = {len(y3)}, < -3: {y3.sum()}")
    print(res.round(3).to_string(index=False))

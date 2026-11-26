"""SS top-1 models for eGFR slope < -3, >= 1.5 yr expanded-RF cohort, kidney-ROI patients (n = 130), with the slope
taken from the LMM fitted on the >= 1 yr cohort (lmm_fit1yr_restrict1p5yr.py ->
case8_min1p5yr_expanded_rf_slope_groups_fit1yr.xlsx) instead of the >= 1.5 yr cohort's own fit.

Everything else is unchanged (same patients, features, CV, seeds, scoring code - functions imported from the original
scripts):
  Binary      binary_ss_topN_lt3_min1p5yr.run, top N = 1, ranking by the binary target and by the 4-cat group
  Ordinal 4   ordinal4_ss1_acr_lt3_min1p5yr.scores, cuts -5/-3/-1: ACR alone, SS top 1, SS top 1 + ACR
  Ordinal 3   same, cuts -3/-1
AUC for slope < -3 (LOOCV, 95 % CI from 2000 stratified bootstrap resamples) and stratified 5-fold x 10 mean AUC.
The original result workbooks are read (not written) for the side-by-side comparison.
Run: C:/Users/tupil/venvs/radiomics_env/python.exe ss_top1_lt3_min1p5yr_fit1yr.py
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from joblib import Parallel, delayed
from sklearn.metrics import roc_auc_score, roc_curve
from sklearn.model_selection import LeaveOneOut, StratifiedKFold

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import ordinal4_topN_lt3_min1p5yr as O  # noqa: E402
import binary_ss_topN_lt3_min1p5yr as B  # noqa: E402
import ordinal4_ss1_acr_lt3_min1p5yr as A  # noqa: E402

RES = O.RES
NEW = RES / "case8_min1p5yr_expanded_rf_slope_groups_fit1yr.xlsx"
S = O.S


def ordinal(Xss, names, acr, slope, scheme):
    cuts = O.SCHEMES[scheme][0]
    y4 = sum((slope < c).astype(int) for c in cuts)
    y3 = (slope < -3).astype(int)
    loo, chosen = A.scores(Xss, names, acr, y4, LeaveOneOut().split(Xss))
    kf = {m: [] for m in A.MODELS}
    for r in range(O.REPS):
        for tr, te in StratifiedKFold(5, shuffle=True, random_state=O.SEED + r).split(Xss, y4):
            s, _ = A.scores(Xss, names, acr, y4, [(tr, te)])
            for m in A.MODELS:
                kf[m].append(roc_auc_score(y3[te], s[m][te]))
    rng = np.random.default_rng(O.SEED)
    i1, i0 = np.where(y3 == 1)[0], np.where(y3 == 0)[0]
    boots = {m: [] for m in A.MODELS}
    for _ in range(O.NBOOT):
        b = np.r_[rng.choice(i1, len(i1)), rng.choice(i0, len(i0))]
        for m in A.MODELS:
            boots[m].append(roc_auc_score(y3[b], loo[m][b]))
    boots = {m: np.array(v) for m, v in boots.items()}
    rows = []
    for m in A.MODELS:
        lo, hi = np.percentile(boots[m], [2.5, 97.5])
        rows.append(dict(model=f"Ordinal {scheme}: {m}", LOOCV_AUC_lt3=roc_auc_score(y3, loo[m]), CI_low=lo, CI_high=hi,
                         LOOCV_Cindex=O.c_index(y4, loo[m]), kfold5x10_mean_AUC_lt3=np.mean(kf[m]),
                         kfold5x10_AUC_IQR=f"{np.percentile(kf[m], 25):.2f}–{np.percentile(kf[m], 75):.2f}",
                         selected=", ".join(f"{f.split('__', 1)[1]} ({v / len(y3):.2f})"
                                            for f, v in pd.Series(chosen).value_counts().head(3).items())
                         if m != "ACR alone" else ""))
    comp = []
    for a, b in [("SS top 1 + ACR", "ACR alone"), ("SS top 1 + ACR", "SS top 1"), ("SS top 1", "ACR alone")]:
        d = boots[a] - boots[b]
        comp.append(dict(comparison=f"Ordinal {scheme}: {a} − {b}",
                         delta_AUC=roc_auc_score(y3, loo[a]) - roc_auc_score(y3, loo[b]),
                         CI_low=np.percentile(d, 2.5), CI_high=np.percentile(d, 97.5),
                         p_two_sided=min(1.0, 2 * min((d <= 0).mean(), (d >= 0).mean()))))
    return rows, comp, {f"Ordinal {scheme}: {m}": loo[m] for m in A.MODELS}


def binary(Xss, names, slope, rank):
    y3 = (slope < -3).astype(int)
    y4 = sum((slope < c).astype(int) for c in O.CUTS4)
    r, s = B.run(Xss, y3, y4, names, rank, 1)
    name = f"Binary: SS top 1 (rank {rank})"
    return [dict(model=name, LOOCV_AUC_lt3=r["LOOCV_AUC_lt3"], CI_low=r["CI_low"], CI_high=r["CI_high"],
                 LOOCV_Cindex=np.nan, kfold5x10_mean_AUC_lt3=r["kfold5x10_mean_AUC_lt3"],
                 kfold5x10_AUC_IQR=r["kfold5x10_AUC_IQR"], selected=r["selected"])], [], {name: s}


def old_results():
    """Same models from the original (slope fitted on >= 1.5 yr) workbooks."""
    rows = []
    b = pd.read_excel(RES / "binary_ss_topN_lt3_min1p5yr.xlsx", "Performance")
    for _, r in b[b.top_n == 1].iterrows():
        rows.append(dict(model=f"Binary: SS top 1 (rank {r.ranking})", old_LOOCV_AUC=r.LOOCV_AUC_lt3,
                         old_CI=f"{r.CI_low:.2f}–{r.CI_high:.2f}", old_kfold5x10=r.kfold5x10_mean_AUC_lt3))
    for scheme, pre in (("4cat", "ordinal4"), ("3cat", "ordinal3")):
        o = pd.read_excel(RES / f"{pre}_ss1_acr_lt3_min1p5yr.xlsx", "Performance")
        for _, r in o.iterrows():
            rows.append(dict(model=f"Ordinal {scheme}: {r.model}", old_LOOCV_AUC=r.LOOCV_AUC_lt3,
                             old_CI=f"{r.CI_low:.2f}–{r.CI_high:.2f}", old_kfold5x10=r.kfold5x10_mean_AUC_lt3))
    return pd.DataFrame(rows)


if __name__ == "__main__":
    X, slope_old = O.load()
    new = pd.read_excel(NEW, "Cohort").set_index("Patient Case ID")
    assert set(X.index) <= set(new.index)
    slope = new[S].reindex(X.index).to_numpy(float)
    assert not np.isnan(slope).any()
    y3 = (slope < -3).astype(int)
    y3_old = (slope_old < -3).astype(int)
    print(f"n = {len(y3)}; < -3: {y3.sum()} (was {y3_old.sum()} with the >= 1.5 yr fit); "
          f"changed: {int((y3 != y3_old).sum())}")
    for scheme in ("4cat", "3cat"):
        cuts = O.SCHEMES[scheme][0]
        print(f"  {scheme} counts new {np.bincount(sum((slope < c).astype(int) for c in cuts)).tolist()}, "
              f"old {np.bincount(sum((slope_old < c).astype(int) for c in cuts)).tolist()}")

    pp = pd.read_excel(A.BASE, "Per patient")
    pp = pp[pp.cohort == "≥ 1.5 yr"].set_index("Patient Case ID")
    acr = np.log(pp.acr.reindex(X.index).to_numpy(float))
    names = [c for c in X.columns if c.startswith("SS__")]
    Xss = X[names].to_numpy(float)

    tasks = [delayed(binary)(Xss, names, slope, r) for r in ("binary", "4-cat")] + \
            [delayed(ordinal)(Xss, names, acr, slope, s) for s in ("4cat", "3cat")]
    out = Parallel(n_jobs=4)(tasks)
    res = pd.DataFrame([r for o in out for r in o[0]])
    comp = pd.DataFrame([c for o in out for c in o[1]])
    loo = {k: v for o in out for k, v in o[2].items()}
    both = res.merge(old_results(), on="model", how="left")
    both["delta_LOOCV_AUC (new - old)"] = both.LOOCV_AUC_lt3 - both.old_LOOCV_AUC

    # ROC: SS top 1 (and ordinal + ACR) with the new slope
    show = ["Binary: SS top 1 (rank binary)", "Ordinal 4cat: SS top 1", "Ordinal 3cat: SS top 1",
            "Ordinal 4cat: SS top 1 + ACR", "Ordinal 3cat: SS top 1 + ACR", "Ordinal 4cat: ACR alone"]
    cols = ["#2a78d6", "#1baf7a", "#6250d6", "#0d5c3f", "#3b2c8f", "#eda100"]
    fig, ax = plt.subplots(figsize=(6.4, 7.6))
    for m, c in zip(show, cols):
        fpr, tpr, _ = roc_curve(y3, loo[m])
        r = res[res.model == m].iloc[0]
        ax.step(fpr, tpr, where="post", color=c, lw=1.8,
                label=f"{m}: {r.LOOCV_AUC_lt3:.2f} ({r.CI_low:.2f}–{r.CI_high:.2f}); 5-fold {r.kfold5x10_mean_AUC_lt3:.2f}")
    ax.plot([0, 1], [0, 1], color="#bdbcb6", lw=1, ls="--")
    ax.set_xlabel("1 − specificity", color=O.INK); ax.set_ylabel("Sensitivity", color=O.INK)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1.01); ax.set_aspect("equal")
    ax.grid(color="#e6e5e1", lw=0.6); ax.set_axisbelow(True)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.tick_params(colors=O.INK2)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.12), fontsize=8.5, frameon=False,
              title="LOOCV AUC (95% CI); 5-fold×10 mean", title_fontsize=8.5)
    ax.set_title(f"SS top 1, slope < −3 — LMM fitted on ≥ 1 yr, reported on ≥ 1.5 yr\n"
                 f"kidney-ROI patients (n = {len(y3)}; < −3: {y3.sum()})", fontsize=11, color=O.INK)
    fig.savefig(RES / "ss_top1_lt3_min1p5yr_fit1yr_roc.png", dpi=300, bbox_inches="tight")
    fig.savefig(RES / "ss_top1_lt3_min1p5yr_fit1yr_roc.pdf", bbox_inches="tight")

    with pd.ExcelWriter(RES / "ss_top1_lt3_min1p5yr_fit1yr.xlsx") as w:
        pd.DataFrame({"Item": ["Method"], "Value": [__doc__]}).to_excel(w, sheet_name="Read_Me", index=False)
        both.to_excel(w, sheet_name="Performance", index=False)
        comp.to_excel(w, sheet_name="Paired differences", index=False)
        pd.DataFrame({**loo, "slope_fit1yr": slope, "slope_fit1p5yr": slope_old, "lt_minus3_fit1yr": y3,
                      "lt_minus3_fit1p5yr": y3_old}, index=X.index).to_excel(w, sheet_name="LOOCV scores")
    pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 120)
    print(both[["model", "LOOCV_AUC_lt3", "CI_low", "CI_high", "kfold5x10_mean_AUC_lt3", "old_LOOCV_AUC", "old_CI",
                "old_kfold5x10", "delta_LOOCV_AUC (new - old)"]].round(3).to_string(index=False))
    print(comp.round(3).to_string(index=False))
    print(res[["model", "selected"]].to_string(index=False))

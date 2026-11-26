"""SS top-1 radiomics + baseline ACR (ordinal 4-category, proportional odds), scored for eGFR slope < -3,
>= 1.5 yr expanded-RF cohort, kidney-ROI patients (n = 130).

Models (same CV / pipeline / score as ordinal4_topN_lt3_min1p5yr.py):
  ACR alone        log(baseline ACR)
  SS top 1         top-1 SS feature by |Spearman| with the 4-cat group, chosen inside each training fold
  SS top 1 + ACR   both; ACR is a forced covariate (not subject to selection)
ACR = baseline urine ACR (mg/mmol) from baseline_egfr_acr_kfre_by_slope.xlsx (sheet 'Per patient', >= 1.5 yr cohort;
      closest to imaging within +-90 d, transplant 3-month rule), log-transformed; training-fold median imputation for
      missing values and training-fold standardisation (nested_loocv_radiomics.forced_block).
Score = held-out linear predictor x.beta without thresholds. AUC for slope < -3 with 95 % CI (2000 stratified
bootstrap resamples) and paired bootstrap differences between models (same resamples). LOOCV + 5-fold x 10 repeats.
Run: C:/Users/tupil/venvs/radiomics_env/python.exe ordinal4_ss1_acr_lt3_min1p5yr.py
"""
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.metrics import roc_auc_score, roc_curve
from sklearn.model_selection import LeaveOneOut, StratifiedKFold

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import ordinal4_topN_lt3_min1p5yr as O  # noqa: E402

sys.path.insert(0, str(O.MODEL_CODE))
from nested_loocv_radiomics import Preprocessor, select_top_spearman_per_modality, forced_block  # noqa: E402

RES = O.RES
BASE = RES / "baseline_egfr_acr_kfre_by_slope.xlsx"
MODELS = ["ACR alone", "SS top 1", "SS top 1 + ACR"]
COL = {"ACR alone": "#eda100", "SS top 1": "#1baf7a", "SS top 1 + ACR": "#2a78d6"}  # slots 4, 3, 1


def scores(Xss, names, acr, y4, splits):
    warnings.filterwarnings("ignore")
    out = {m: np.full(len(y4), np.nan) for m in MODELS}
    chosen = []
    for tr, te in splits:
        prep, xtr = Preprocessor.fit(Xss[tr], names, 0.90, y4[tr], "spearman")
        xte = prep.transform(Xss[te])
        idx, _ = select_top_spearman_per_modality(xtr, y4[tr], prep.selected_features, 1, True)
        chosen.append(prep.selected_features[idx[0]])
        a_tr, a_te = forced_block(acr[tr, None], acr[te, None])
        out["ACR alone"][te] = O.fit_score(a_tr, y4[tr], a_te)
        out["SS top 1"][te] = O.fit_score(xtr[:, idx], y4[tr], xte[:, idx])
        out["SS top 1 + ACR"][te] = O.fit_score(np.hstack([xtr[:, idx], a_tr]), y4[tr], np.hstack([xte[:, idx], a_te]))
    return out, chosen


if __name__ == "__main__":
    scheme = sys.argv[1] if len(sys.argv) > 1 else "4cat"   # 4cat (default) | 3cat
    cuts, prefix, slabel = O.SCHEMES[scheme]
    X, slope = O.load()
    y4 = sum((slope < c).astype(int) for c in cuts)
    y3 = (slope < -3).astype(int)
    pp = pd.read_excel(BASE, "Per patient")
    pp = pp[pp.cohort == "≥ 1.5 yr"].set_index("Patient Case ID")
    acr_raw = pp.acr.reindex(X.index)
    acr = np.log(acr_raw.to_numpy(float))
    print(f"n = {len(y3)}, < -3: {y3.sum()}, ACR missing: {int(np.isnan(acr).sum())}")
    names = [c for c in X.columns if c.startswith("SS__")]
    Xss = X[names].to_numpy(float)

    loo, chosen = scores(Xss, names, acr, y4, LeaveOneOut().split(Xss))
    kf = {m: [] for m in MODELS}
    for r in range(O.REPS):
        for tr, te in StratifiedKFold(5, shuffle=True, random_state=O.SEED + r).split(Xss, y4):
            s, _ = scores(Xss, names, acr, y4, [(tr, te)])
            for m in MODELS:
                kf[m].append(roc_auc_score(y3[te], s[m][te]))

    rng = np.random.default_rng(O.SEED)
    i1, i0 = np.where(y3 == 1)[0], np.where(y3 == 0)[0]
    boots = {m: [] for m in MODELS}
    for _ in range(O.NBOOT):
        b = np.r_[rng.choice(i1, len(i1)), rng.choice(i0, len(i0))]
        for m in MODELS:
            boots[m].append(roc_auc_score(y3[b], loo[m][b]))
    boots = {m: np.array(v) for m, v in boots.items()}
    rows = []
    for m in MODELS:
        lo, hi = np.percentile(boots[m], [2.5, 97.5])
        rows.append(dict(model=m, LOOCV_AUC_lt3=roc_auc_score(y3, loo[m]), CI_low=lo, CI_high=hi,
                         **{f"LOOCV_Cindex_{scheme}": O.c_index(y4, loo[m])}, kfold5x10_mean_AUC_lt3=np.mean(kf[m]),
                         kfold5x10_AUC_IQR=f"{np.percentile(kf[m], 25):.2f}–{np.percentile(kf[m], 75):.2f}"))
    res = pd.DataFrame(rows)
    comp = []
    for a, b in [("SS top 1 + ACR", "ACR alone"), ("SS top 1 + ACR", "SS top 1"), ("SS top 1", "ACR alone")]:
        d = boots[a] - boots[b]
        comp.append(dict(comparison=f"{a} − {b}", delta_AUC=roc_auc_score(y3, loo[a]) - roc_auc_score(y3, loo[b]),
                         CI_low=np.percentile(d, 2.5), CI_high=np.percentile(d, 97.5),
                         p_two_sided=min(1.0, 2 * min((d <= 0).mean(), (d >= 0).mean()))))
    comp = pd.DataFrame(comp)
    feat = pd.Series(chosen).value_counts()

    fig, ax = plt.subplots(figsize=(6.2, 6.8))
    for m in MODELS:
        fpr, tpr, _ = roc_curve(y3, loo[m])
        r = res[res.model == m].iloc[0]
        ax.step(fpr, tpr, where="post", color=COL[m], lw=2,
                label=f"{m}: {r.LOOCV_AUC_lt3:.2f} ({r.CI_low:.2f}–{r.CI_high:.2f}); 5-fold {r.kfold5x10_mean_AUC_lt3:.2f}")
    ax.plot([0, 1], [0, 1], color="#bdbcb6", lw=1, ls="--")
    ax.set_xlabel("1 − specificity", color=O.INK); ax.set_ylabel("Sensitivity", color=O.INK)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1.01); ax.set_aspect("equal")
    ax.grid(color="#e6e5e1", lw=0.6); ax.set_axisbelow(True)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.tick_params(colors=O.INK2)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.12), fontsize=9, frameon=False,
              title="LOOCV AUC (95% CI); 5-fold×10 mean", title_fontsize=9)
    ax.set_title(f"Ordinal {slabel.split(' (')[0]} model scored for slope < −3\n≥ 1.5 yr, kidney-ROI patients "
                 f"(n = {len(y3)}; < −3: {y3.sum()})", fontsize=11, color=O.INK)
    fig.savefig(RES / f"{prefix}_ss1_acr_lt3_min1p5yr_roc.png", dpi=300, bbox_inches="tight")
    fig.savefig(RES / f"{prefix}_ss1_acr_lt3_min1p5yr_roc.pdf", bbox_inches="tight")

    with pd.ExcelWriter(RES / f"{prefix}_ss1_acr_lt3_min1p5yr.xlsx") as w:
        pd.DataFrame({"Item": ["Method"], "Value": [__doc__]}).to_excel(w, sheet_name="Read_Me", index=False)
        res.to_excel(w, sheet_name="Performance", index=False)
        comp.to_excel(w, sheet_name="Paired differences", index=False)
        feat.rename("LOOCV folds").to_frame().to_excel(w, sheet_name="SS feature chosen")
        pd.DataFrame({**{m: loo[m] for m in MODELS}, "acr_mg_mmol": acr_raw.to_numpy(), "slope": slope,
                      "group": y4, "lt_minus3": y3}, index=X.index).to_excel(w, sheet_name="LOOCV scores")
    pd.set_option("display.width", 200)
    print(res.round(3).to_string(index=False)); print(comp.round(3).to_string(index=False)); print(feat.to_string())

#!/usr/bin/env python3
"""Clinical add-on (eGFR, ACR, KFRE: all 7 non-empty combinations) for the best all-combinations radiomics models,
>= 1.5 yr expanded-RF cohort (n = 130), LOOCV, focus slope < -3.

Radiomics models: the 10 best unpenalised top-N proportional-odds runs of results_<scheme>/summary_all_runs.csv by
  held-out AUC for slope < -3 (runs with identical held-out scores, e.g. the same single feature reached through
  different sets, are kept once).
Per LOOCV training fold: radiomics features exactly as the grid (reference run_all_combinations.prepare: median
  imputation -> |r| >= 0.9 removal by |Spearman| -> standardise -> top-N by |Spearman rho| with the ordinal group);
  clinical covariates forced in (nested_loocv_radiomics.forced_block: training-fold median imputation + standardisation).
Clinical (as 01 clinical_addon/run_clinical_addon.py): ACR = ln(baseline ACR mg/mmol, floored 0.1);
  eGFR = baseline CKD-EPI 2021; KFRE = Tangri 4-variable linear predictor ln(ln(1 - risk_2yr) / ln 0.975).
Model: unpenalised proportional-odds ordinal LR (statsmodels OrderedModel) on the scheme's groups.
Score: held-out x.beta (no intercepts). AUC for slope < -3, stratified bootstrap 95 % CI (2,000), paired bootstrap
  differences vs radiomics only and vs the same clinical set alone (same resamples).
Usage: C:/Users/tupil/venvs/radiomics_env/python.exe clinical_addon.py 4cat|3cat
"""
import itertools
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.model_selection import LeaveOneOut

ROOT = Path(__file__).resolve().parents[2]
REF = ROOT / "Results_Rough" / "egfr_slope_all_combinations"
sys.path.insert(0, str(REF))
import run_all_combinations as R  # noqa: E402
from nested_loocv_radiomics import forced_block  # noqa: E402

OUTROOT = ROOT / "Results_Rough" / "case8_min1p5yr_all_combinations"
SCHEMES = {"4cat": (-5, -3, -1), "3cat": (-3, -1)}
METHODS = {"top1_LR": (1, False), "top2_LR": (2, False), "top3_LR": (3, False), "top5_LR": (5, False),
           "top1overall_LR": (1, True), "top2overall_LR": (2, True), "top3overall_LR": (3, True), "top5overall_LR": (5, True)}
CLIN_VARS = ["eGFR", "ACR", "KFRE"]
CLIN_SETS = [list(c) for k in (1, 2, 3) for c in itertools.combinations(CLIN_VARS, k)]
TOP_K, N_BOOT, SEED = 10, 2000, 2026


def clin_name(s):
    return " + ".join(s)


def main():
    warnings.simplefilter("ignore"); np.seterr(all="ignore")
    scheme = sys.argv[1]
    cuts = SCHEMES[scheme]
    out = OUTROOT / f"results_{scheme}"
    F = pd.read_csv(OUTROOT / "inputs" / "patient_features_6mod.csv")
    names_all = [c for c in F.columns if "__" in c]
    X = F[names_all].to_numpy(float)
    slope = F.slope.to_numpy(float)
    y = sum((slope < c).astype(int) for c in cuts)
    lab = slope < -3
    risk = (F.kfre2_pct / 100).clip(1e-9, 1 - 1e-9)
    C = pd.DataFrame({"eGFR": F.egfr, "ACR": np.log(F.acr_mg_mmol.clip(lower=0.1)),
                      "KFRE": np.log(np.log(1 - risk) / np.log(0.975))}).to_numpy(float)

    summ = pd.read_csv(out / "summary_all_runs.csv")
    preds = pd.read_csv(out / "heldout_predictions_all_runs.csv.gz")
    cand = summ[summ.method.isin(METHODS)].sort_values("AUC slope < -3", ascending=False)
    chosen, seen = [], []
    for _, r in cand.iterrows():
        s = preds[(preds.combination == r.combination) & (preds.method == r.method)].heldout_score.to_numpy()
        if any(np.allclose(s, t, equal_nan=True) for t in seen):
            continue
        seen.append(s); chosen.append((r.combination, r.method))
        if len(chosen) == TOP_K:
            break
    print(f"[{scheme}] radiomics models: {chosen}")

    scores = {}
    loo = list(LeaveOneOut().split(X))
    for cs in CLIN_SETS:  # clinical only
        idx = [CLIN_VARS.index(v) for v in cs]
        s = np.full(len(y), np.nan)
        for tr, te in loo:
            a, b = forced_block(C[tr][:, idx], C[te][:, idx])
            s[te] = R.fit_po(a, y[tr], b)[0]
        scores[("Clinical only", clin_name(cs))] = s
    for combo, method in chosen:
        top_n, overall = METHODS[method]
        cols = [i for i, n in enumerate(names_all) if n.split("__", 1)[0] in combo.split("+")]
        x, names = X[:, cols], [names_all[i] for i in cols]
        model = f"{combo.replace('HScan', 'H-scan').replace('+', ' + ')} | {method}"
        S = {v: np.full(len(y), np.nan) for v in ["Radiomics only"] + [clin_name(c) for c in CLIN_SETS]}
        for tr, te in loo:
            xt, xs, _, _, _ = R.prepare(x[tr], y[tr], x[te], names, top_n, overall)
            S["Radiomics only"][te] = R.fit_po(xt, y[tr], xs)[0]
            for cs in CLIN_SETS:
                idx = [CLIN_VARS.index(v) for v in cs]
                a, b = forced_block(C[tr][:, idx], C[te][:, idx])
                S[clin_name(cs)][te] = R.fit_po(np.hstack([xt, a]), y[tr], np.hstack([xs, b]))[0]
        ref = preds[(preds.combination == combo) & (preds.method == method)].heldout_score.to_numpy()
        assert np.allclose(S["Radiomics only"], ref, equal_nan=True), f"radiomics-only differs from grid for {model}"
        for v, s in S.items():
            scores[(model, v)] = s

    rng = np.random.default_rng(SEED)
    i1, i0 = np.where(lab)[0], np.where(~lab)[0]
    keys = list(scores)
    B = np.empty((N_BOOT, len(keys)))
    for b in range(N_BOOT):
        i = np.r_[rng.choice(i1, i1.size), rng.choice(i0, i0.size)]
        B[b] = [R.auc(lab[i], scores[k][i]) for k in keys]
    rows = []
    for j, (model, var) in enumerate(keys):
        s = scores[(model, var)]
        row = dict(model=model, variant=var, AUC_lt3=R.auc(lab, s), ci95_low=np.percentile(B[:, j], 2.5),
                   ci95_high=np.percentile(B[:, j], 97.5), c_index=R.c_index(y, s))
        if model != "Clinical only" and var != "Radiomics only":
            for refk, tag in (((model, "Radiomics only"), "vs radiomics only"), (("Clinical only", var), "vs clinical only")):
                d = B[:, j] - B[:, keys.index(refk)]
                row[f"delta {tag}"] = row["AUC_lt3"] - R.auc(lab, scores[refk])
                row[f"delta {tag} ci95"] = f"{np.percentile(d, 2.5):.3f} to {np.percentile(d, 97.5):.3f}"
                row[f"p {tag}"] = min(1.0, 2 * min((d <= 0).mean(), (d >= 0).mean()))
        rows.append(row)
    res = pd.DataFrame(rows)
    grid = res.pivot(index="model", columns="variant", values="AUC_lt3")
    order_rows = ["Clinical only"] + [m for m in dict.fromkeys(res.model) if m != "Clinical only"]
    order_cols = ["Radiomics only"] + [clin_name(c) for c in CLIN_SETS]
    grid = grid.reindex(index=order_rows, columns=order_cols)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    M = grid.to_numpy(float)
    fig, ax = plt.subplots(figsize=(13, 0.5 * len(order_rows) + 2.5))
    im = ax.imshow(M, cmap="RdBu_r", vmin=0.3, vmax=0.85, aspect="auto")
    for i in range(M.shape[0]):
        for jj in range(M.shape[1]):
            if np.isfinite(M[i, jj]):
                ax.text(jj, i, f"{M[i, jj]:.2f}", ha="center", va="center", fontsize=8,
                        color="white" if abs(M[i, jj] - 0.575) > 0.2 else "black", fontweight="bold" if M[i, jj] >= 0.75 else "normal")
    ax.set_xticks(range(len(order_cols)), [c.replace(" + ", "\n+ ") for c in order_cols], fontsize=8); ax.xaxis.tick_top()
    ax.set_yticks(range(len(order_rows)), order_rows, fontsize=8)
    ax.axhline(0.5, color="black", lw=1.5)
    fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02).set_label("held-out AUC for slope < −3", fontsize=8)
    ax.set_title(f"Clinical add-on (forced in), ordinal {scheme} model, LOOCV AUC for eGFR slope < −3 "
                 f"(n = {len(y)}; < −3: {lab.sum()}); bold ≥ 0.75", fontsize=10, pad=30)
    fig.savefig(out / f"clinical_addon_auc_lt3_{scheme}.png", dpi=200, bbox_inches="tight")
    fig.savefig(out / f"clinical_addon_auc_lt3_{scheme}.pdf", bbox_inches="tight")
    with pd.ExcelWriter(out / f"clinical_addon_{scheme}.xlsx") as w:
        pd.DataFrame({"Item": ["Description"], "Value": [__doc__]}).to_excel(w, sheet_name="Read_Me", index=False)
        grid.to_excel(w, sheet_name="Grid_AUC_lt-3")
        res.to_excel(w, sheet_name="All_with_CI", index=False)
        pd.DataFrame({f"{m} | {v}": s for (m, v), s in scores.items()}, index=F.patient_id) \
            .assign(slope=slope, group=y).to_excel(w, sheet_name="Heldout_scores")
    pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 60)
    print(grid.round(3).to_string())


if __name__ == "__main__":
    main()

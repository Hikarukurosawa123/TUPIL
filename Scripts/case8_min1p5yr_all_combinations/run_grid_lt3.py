#!/usr/bin/env python3
"""All-combinations grid (6 QUS modalities, 63 combinations), >= 1.5 yr expanded-RF Case 8 cohort (n = 130), LOOCV.
Focus: eGFR slope < -3.

The model code is the reference Results_Rough/egfr_slope_all_combinations/run_all_combinations.py, imported and called
unchanged (run_one: in-fold median imputation -> |r| >= 0.9 removal keeping the higher training |Spearman rho| ->
standardisation -> top-N by |Spearman rho| with the ordinal group -> model; held-out score = x.beta without intercepts).
Differences from the reference run:
  * cohort / features: inputs/patient_features_6mod.csv (prepare_inputs.py, n = 130, 4,463 features)
  * target scheme (ranking target for Spearman AND ordinal model):
      4cat: >= -1 | -3 to -1 | -5 to -3 | < -5      3cat: >= -1 | -3 to -1 | < -3
    (N_CLASSES is set inside every worker for the penalised cumulative models and the bootstrap)
  * methods: the reference 8 + top2_LR (per modality) and top2overall_LR
  * main metric: AUC of the held-out score for slope < -3 with a stratified bootstrap 95 % CI (2,000 resamples);
    ordinal C-index (reference metric) and AUC < -1 / < -5 also reported.
Usage: C:/Users/tupil/venvs/radiomics_env/python.exe run_grid_lt3.py 4cat|3cat   (env N_JOBS, default 20;
       env ONLY="SS:top1_LR" runs a single combination:method as a check and prints its summary)
"""
import itertools
import json
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
REF = ROOT / "Results_Rough" / "egfr_slope_all_combinations"
HERE = Path(__file__).resolve().parent
OUTROOT = ROOT / "Results_Rough" / "case8_min1p5yr_all_combinations"
FEATURES = OUTROOT / "inputs" / "patient_features_6mod.csv"
for p in (str(HERE), str(REF)):   # REF ends up first: the reference module wins any name lookup
    if p not in sys.path:
        sys.path.insert(0, p)
os.environ["PYTHONPATH"] = os.pathsep.join([str(REF), str(HERE), os.environ.get("PYTHONPATH", "")])
import run_all_combinations as R  # noqa: E402  (reference module)

SCHEMES = {"4cat": (-5, -3, -1), "3cat": (-3, -1)}
EXTRA_METHODS = {"top2_LR": dict(top_n=2, penalty=None, overall=False),
                 "top2overall_LR": dict(top_n=2, penalty=None, overall=True)}
METHOD_ORDER = ["top1_LR", "top2_LR", "top3_LR", "top5_LR", "top1overall_LR", "top2overall_LR", "top3overall_LR",
                "top5overall_LR", "top20_LASSO", "top20_RIDGE"]
N_BOOT, SEED = 2000, 2026


def auc_ci(label, score, seed):
    label = np.asarray(label, bool)
    rng = np.random.default_rng(seed)
    i1, i0 = np.where(label)[0], np.where(~label)[0]
    v = [R.auc(label[i], score[i]) for i in (np.r_[rng.choice(i1, i1.size), rng.choice(i0, i0.size)] for _ in range(N_BOOT))]
    return tuple(np.percentile(v, [2.5, 97.5]))


def task(combo, method, x, names, y, slope, ids, n_classes):
    import run_all_combinations as RR  # re-imported in the worker
    RR.N_CLASSES = n_classes
    RR.METHODS.update(EXTRA_METHODS)
    summary, pred, sel, tune = RR.run_one(combo, method, x, names, y, slope, ids)
    s = pred.heldout_score.to_numpy(float)
    ok = np.isfinite(s)
    lo, hi = auc_ci(slope[ok] < -3, s[ok], SEED + sum(map(ord, "+".join(combo) + method)))
    summary["AUC < -3 ci95_low"], summary["AUC < -3 ci95_high"] = lo, hi
    return summary, pred, sel, tune


def heatmap(summary, combos, value, path, title, vmin=0.3, vmax=0.8):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    rows = ["+".join(c) for c in combos]
    M = summary.pivot(index="combination", columns="method", values=value).loc[rows, METHOD_ORDER].to_numpy()
    fig, ax = plt.subplots(figsize=(12, 0.26 * len(rows) + 2.5))
    im = ax.imshow(M, cmap="RdBu_r", vmin=vmin, vmax=vmax, aspect="auto")
    best = np.nanmax(M)
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            v = M[i, j]
            ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=6.5,
                    color="white" if abs(v - (vmin + vmax) / 2) > 0.3 * (vmax - vmin) else "black",
                    fontweight="bold" if v >= 0.70 else "normal")
            if v == best:
                ax.add_patch(plt.Rectangle((j - .5, i - .5), 1, 1, fill=False, ec="k", lw=1.5))
    ax.set_xticks(range(len(METHOD_ORDER)))
    ax.set_xticklabels([m.replace("_", "\n") for m in METHOD_ORDER], fontsize=8)
    ax.xaxis.tick_top()
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([r.replace("HScan", "H-scan").replace("+", " + ") for r in rows], fontsize=7)
    for b in np.where(np.diff([len(c) for c in combos]))[0]:
        ax.axhline(b + 0.5, color="black", lw=1)
    ax.axvline(3.5, color="black", lw=1); ax.axvline(7.5, color="black", lw=1)
    cb = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    cb.set_label(value, fontsize=8)
    ax.set_title(title, fontsize=9, pad=28)
    fig.savefig(path.with_suffix(".png"), dpi=200, bbox_inches="tight")
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def main():
    from joblib import Parallel, delayed
    scheme = sys.argv[1]
    cuts = SCHEMES[scheme]
    n_classes = len(cuts) + 1
    out = OUTROOT / f"results_{scheme}"
    out.mkdir(parents=True, exist_ok=True)
    F = pd.read_csv(FEATURES)
    names = [c for c in F.columns if "__" in c]
    x = F[names].to_numpy(float)
    slope = F["slope"].to_numpy(float)
    y = sum((slope < c).astype(int) for c in cuts)
    ids = F["patient_id"].astype(str).to_numpy()
    combos = [c for k in range(1, len(R.MODALITIES) + 1) for c in itertools.combinations(R.MODALITIES, k)]
    tasks = [(c, m) for c in combos for m in METHOD_ORDER]
    tasks.sort(key=lambda t: ("LASSO" not in t[1] and "RIDGE" not in t[1], -len(t[0])))  # slow nested runs first
    only = os.environ.get("ONLY")
    if only:
        cn, mn = only.split(":")
        t1 = time.time()
        r = task(tuple(cn.split("+")), mn, x, names, y, slope, ids, n_classes)
        print(pd.Series(r[0]).to_string()); print(f"{time.time() - t1:.1f} s"); return
    print(f"[{scheme}] {len(F)} patients, {len(names)} features, groups {np.bincount(y).tolist()}; "
          f"{len(combos)} combinations x {len(METHOD_ORDER)} methods -> {len(tasks)} runs", flush=True)
    t0 = time.time()
    res = Parallel(n_jobs=int(os.environ.get("N_JOBS", 20)), verbose=5)(
        delayed(task)(c, m, x, names, y, slope, ids, n_classes) for c, m in tasks)
    summary = pd.DataFrame([r[0] for r in res])
    summary["_c"] = summary["combination"].map({"+".join(c): i for i, c in enumerate(combos)})
    summary["_m"] = summary["method"].map({m: i for i, m in enumerate(METHOD_ORDER)})
    summary = summary.sort_values(["_c", "_m"]).drop(columns=["_c", "_m"]).reset_index(drop=True)
    summary.insert(0, "scheme", scheme)
    preds = pd.concat([r[1] for r in res], ignore_index=True)
    sel = pd.concat([r[2] for r in res], ignore_index=True)
    tune = pd.concat([r[3] for r in res if len(r[3])], ignore_index=True)
    summary.to_csv(out / "summary_all_runs.csv", index=False)
    preds.to_csv(out / "heldout_predictions_all_runs.csv.gz", index=False)
    sel.to_csv(out / "selected_features_coefficients_all_runs.csv.gz", index=False)
    tune.to_csv(out / "inner_cv_tuning_nested_runs.csv.gz", index=False)
    rows = ["+".join(c) for c in combos]
    with pd.ExcelWriter(out / f"all_combinations_{scheme}_results.xlsx") as w:
        pd.DataFrame({"Item": ["Description"], "Value": [__doc__]}).to_excel(w, sheet_name="Read_Me", index=False)
        summary.sort_values("AUC slope < -3", ascending=False).to_excel(w, sheet_name="Ranked_by_AUC_lt-3", index=False)
        for v, sh in (("AUC slope < -3", "AUC_lt-3"), ("c_index", "C_index"), ("AUC slope < -1", "AUC_lt-1"),
                      ("AUC slope < -5", "AUC_lt-5"), ("c_index_rscore", "C_index_Rscore")):
            summary.pivot(index="combination", columns="method", values=v).loc[rows, METHOD_ORDER] \
                .dropna(axis=1, how="all").to_excel(w, sheet_name=f"Grid_{sh}")
        summary.to_excel(w, sheet_name="All_runs", index=False)
    n, lab = len(y), {"4cat": "4 groups (≥−1 | −3 to −1 | −5 to −3 | <−5)", "3cat": "3 groups (≥−1 | −3 to −1 | <−3)"}[scheme]
    heatmap(summary, combos, "AUC slope < -3", out / f"heatmap_auc_lt3_{scheme}",
            f"Held-out AUC for eGFR slope < −3 (LOOCV, n = {n}; < −3: {(slope < -3).sum()}); ordinal model on {lab}\n"
            "In-fold: |r| ≥ 0.9 removal (higher |Spearman ρ| kept) → top-N by |ρ| with the ordinal group; "
            "PO ordinal LR (top-N) or nested 5-fold LASSO/ridge (top 20 per modality); bold ≥ 0.70; box = best")
    heatmap(summary, combos, "c_index", out / f"heatmap_c_index_{scheme}",
            f"Held-out ordinal C-index, {lab} (LOOCV, n = {n}); bold ≥ 0.70; box = best")
    settings = {"scheme": scheme, "cuts": cuts, "n_classes": n_classes, "modalities": R.MODALITIES,
                "methods": {**R.METHODS, **EXTRA_METHODS}, "outer_cv": f"LeaveOneOut ({n} folds)",
                "inner_cv_nested": f"StratifiedKFold({R.INNER_FOLDS}, shuffle=True, random_state={R.SEED}+outer_fold)",
                "c_grid": R.C_GRID, "corr_removal": f"|Pearson r| >= {R.CORR_THRESHOLD}, greedy by training |Spearman rho|",
                "n_boot": N_BOOT, "total_minutes": round((time.time() - t0) / 60, 1)}
    (out / "settings.json").write_text(json.dumps(settings, indent=2, default=str), encoding="utf-8")
    print(summary.sort_values("AUC slope < -3", ascending=False)[
        ["combination", "method", "AUC slope < -3", "AUC < -3 ci95_low", "AUC < -3 ci95_high", "c_index"]].head(25).to_string(index=False))
    print(f"done in {settings['total_minutes']} min", flush=True)


if __name__ == "__main__":
    main()

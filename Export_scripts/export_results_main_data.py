"""Data + statistics for the Results_Main figures (drawn by make_results_main_figures.m).

    python export_results_main_data.py

Output: Results_Main/<section>/*.csv  (figures are written next to them by the MATLAB script)

01_Rscore                  held-out R-score by ci (box, Spearman rho)
                           src: Results_Rough/02_ci_Rscore_radiomics_top3overall_SI_HScan
                           quadratic weighted kappa (confusion): SI + H-scan top 3 proportional-odds model, CV test
                           src: POCUS_Kidney_Analysis/.../Rscore_ordinal_SI3_HScan3_results.xlsx (kidney_pr_all), as in
                                v10_figures/rscore_banff/weighted_kappa
02_Rscore_all_combinations held-out Spearman rho / QWK / AUC ci0 vs ci1-2 / AUC ci0 vs ci1-2+HD, 63 combinations x 8
                           selection strategies (stratified 5-fold, n = 174)
                           src: Results_Rough/ci_R_score_all_combinations/results/pivot_*.csv
03_eGFR_slope              LOOCV held-out C-index and AUC (< -1, < -3, < -5), 63 combinations x 8 strategies;
                           best model (H-scan + ESD + SS, top 1 per modality) held-out score, radiomics only and
                           radiomics + ln ACR; baseline ACR / eGFR / KFRE by slope group (4-group, binary < -1/-3/-5)
                           src: Results_Rough/egfr_slope_all_combinations (+ addon_acr),
                                Results_Rough/01_eGFR_slope_model3_min1yr/.../baseline_clinical_model3_91.csv
04_ci_spatial_median_gap   AUC (ci 0 vs 1-2) by ROI x early / late half (median gap 44 d), frozen + ordinal logistic;
                           cross-ROI AUC (whole-cohort CV) and cross-time transfer (early <-> late) per half
                           src: Results_Rough/Final_100/spatial_Final_100.xlsx,
                                Results_Rough/Final_100/cross_analysis_median_gap/spatial_cross_Final_100_median.xlsx
05_ci_spatial_time_windows AUC by ROI x interval / cumulative window; cross-ROI AUC (whole-cohort CV) for the
                           intervals and the cumulative windows
                           src: Results_Rough/Final_100/ci_spatial_0_1_3_5_months (+ cross_analysis),
                                Results_Rough/Final_100/cross_analysis_median_gap

Box plots: Welch one-way ANOVA + pairwise comparisons (>= 3 groups) or Welch t-test (2 groups).
Pairwise: Tukey-Kramer for 01 (R-score); unadjusted Welch t-tests for 03 (eGFR slope).
ACR and KFRE are tested on log10 (ACR floored at 0.1 mg/mmol, as in the ln-ACR model term).
Spatial AUC tests are taken as reported in the source workbooks (DeLong, one-sided H1 first > second; ROI
comparisons paired, between-window comparisons unpaired; uncorrected).
"""
from itertools import combinations
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy import stats

REPO = Path(__file__).resolve().parents[1]
RR = REPO / "Results_Rough"
OUT = REPO / "Results_Main"
KAPPA_MODEL = (REPO.parents[0] / "POCUS_Kidney_Analysis" / "Analysis_scripts" / "Redo_analysis_1106_2026" /
               "ci_spearman_radiomics_final" / "cohort174_firstVisit_2p5MHz_20261113" /
               "Rscore_ordinal_SI_HScan_top3overall" / "Rscore_ordinal_SI3_HScan3_results.xlsx")
O1, O2, O3, O4, O5 = (OUT / d for d in ("01_Rscore", "02_Rscore_all_combinations", "03_eGFR_slope",
                                       "04_ci_spatial_median_gap", "05_ci_spatial_time_windows"))

STRATS = [("top1", "Per modality", "1"), ("top3", "Per modality", "3"), ("top5", "Per modality", "5"),
          ("top1overall", "Overall", "1"), ("top3overall", "Overall", "3"), ("top5overall", "Overall", "5"),
          ("top20_LASSO", "Top 20", "LASSO"), ("top20_RIDGE", "Top 20", "Ridge")]
MOD_NAME = {"HScan": "H-scan", "H-scan": "H-scan"}


# ---------------------------------------------------------------- statistics
def welch_anova(groups):
    g = [np.asarray(x, float) for x in groups]
    k = len(g); n = np.array([len(x) for x in g]); m = np.array([x.mean() for x in g])
    w = n / np.array([x.var(ddof=1) for x in g]); W = w.sum(); mw = (w * m).sum() / W
    lam = ((1 - w / W) ** 2 / (n - 1)).sum()
    F = ((w * (m - mw) ** 2).sum() / (k - 1)) / (1 + 2 * (k - 2) / (k ** 2 - 1) * lam)
    df2 = (k ** 2 - 1) / (3 * lam)
    return F, k - 1, df2, stats.f.sf(F, k - 1, df2)


def group_tests(groups, labels, adjust=True):
    """>= 3 groups: Welch ANOVA + Tukey-Kramer (adjust=True) or pairwise Welch t-tests without
    multiple-comparison adjustment (adjust=False); 2 groups: Welch t-test. Column 'p' drives the brackets."""
    groups = [np.asarray(x, float) for x in groups]
    pairs = list(combinations(range(len(groups)), 2))
    mw = [stats.mannwhitneyu(groups[a], groups[b]).pvalue for a, b in pairs]
    if len(groups) >= 3:
        F, d1, d2, p_all = welch_anova(groups)
        if adjust:
            tk = stats.tukey_hsd(*groups).pvalue
            pp = [tk[a, b] for a, b in pairs]; test = "Tukey-Kramer"
        else:
            pp = [stats.ttest_ind(groups[a], groups[b], equal_var=False).pvalue for a, b in pairs]
            test = "Welch t-test (unadjusted)"
        overall = {"Welch ANOVA F": F, "df1": d1, "df2": d2, "Welch ANOVA p": p_all,
                   "Kruskal-Wallis p": stats.kruskal(*groups).pvalue}
    else:
        t = stats.ttest_ind(groups[0], groups[1], equal_var=False)
        pp = [t.pvalue]; test = "Welch t-test"
        overall = {"Welch t": t.statistic, "Welch t-test p": t.pvalue, "Mann-Whitney p": mw[0]}
    rows = [{"group A": labels[a], "group B": labels[b], "pos A": a + 1, "pos B": b + 1,
             "n A": len(groups[a]), "n B": len(groups[b]),
             "median A": np.median(groups[a]), "median B": np.median(groups[b]),
             "test": test, "p": pp[i], "Mann-Whitney p": mw[i]} for i, (a, b) in enumerate(pairs)]
    return overall, pd.DataFrame(rows)


def c_index(y, s):
    num = den = 0.0
    for lo in np.unique(y):
        for hi in np.unique(y):
            if hi > lo:
                d = s[y == hi][:, None] - s[y == lo][None, :]
                num += (d > 0).sum() + 0.5 * (d == 0).sum(); den += d.size
    return num / den


def nice_combo(c):
    return " + ".join(MOD_NAME.get(m.strip(), m.strip()) for m in str(c).replace(" + ", "+").split("+"))


def long_grid(pivots, metric_labels):
    """pivots: {metric key: DataFrame(combination x strategy columns)} -> tidy rows in fixed order."""
    rows = []
    for key, P in pivots.items():
        for r, c in enumerate(P.combination):
            for j, (col, grp, lab) in enumerate(STRATS):
                rows.append({"metric": key, "metric label": metric_labels[key], "row": r + 1, "col": j + 1,
                             "combination": nice_combo(c), "n modalities": len(nice_combo(c).split(" + ")),
                             "strategy": col, "strategy group": grp, "strategy label": lab,
                             "value": float(P.loc[P.combination == c, col].iloc[0])})
    return pd.DataFrame(rows)


# ---------------------------------------------------------------- 01 R-score
K = 4


def confusion(y, x):
    O = np.zeros((K, K)); np.add.at(O, (y, x), 1); return O


def kappa(y, x, weights="quadratic"):
    O = confusion(y, x); E = np.outer(O.sum(1), O.sum(0)) / O.sum()
    i, j = np.indices((K, K))
    W = {"quadratic": (i - j) ** 2, "linear": np.abs(i - j), None: (i != j).astype(float)}[weights]
    return 1 - (W * O).sum() / (W * E).sum()


def kappa_null_test(y, x):
    O = confusion(y, x); n = O.sum(); p = O / n; pr, pc = p.sum(1), p.sum(0)
    i, j = np.indices((K, K)); w = 1 - (i - j) ** 2 / (K - 1) ** 2
    po, pe = (w * p).sum(), (w * np.outer(pr, pc)).sum(); k = (po - pe) / (1 - pe)
    wr, wc = w @ pc, pr @ w
    var0 = ((np.outer(pr, pc) * (w - (wr[:, None] + wc[None, :])) ** 2).sum() - pe ** 2) / (n * (1 - pe) ** 2)
    z = k / np.sqrt(var0)
    return z, 2 * stats.norm.sf(abs(z))


def section01():
    O1.mkdir(parents=True, exist_ok=True)
    # held-out R-score of the SI + H-scan top 3 overall proportional-odds model (kidney_pr_all, 10 x 5-fold CV);
    # same model and condition as the weighted-kappa confusion matrix below.
    pm = pd.read_excel(KAPPA_MODEL, sheet_name="R-score per visit")
    pm = pm[pm["condition"] == "kidney_pr_all"].reset_index(drop=True)
    pm = pm.rename(columns={"ci": "true_index", "R-score CV test": "r_score"})
    y = pm.true_index.to_numpy(int); s = pm.r_score.to_numpy(float)
    assert len(pm) == 174 and set(np.unique(y)) == {0, 1, 2, 3} and np.isfinite(s).all()
    pm[["patient", "scan_id", "true_index", "r_score", "predicted ci (CV)"]].to_csv(
        O1 / "rscore_heldout_per_patient.csv", index=False)

    labels = ["ci 0", "ci 1", "ci 2", "HD"]
    overall, pairs = group_tests([s[y == k] for k in range(4)], labels)
    rho, prho = stats.spearmanr(s, y)
    pairs.to_csv(O1 / "rscore_heldout_by_ci_pairwise.csv", index=False)
    pd.DataFrame([{"score": "held-out R-score (proportional-odds model, SI + H-scan top 3 overall, kidney_pr_all, "
                            "10 x 5-fold CV)", "n": len(y),
                   **{f"n {l}": int((y == k).sum()) for k, l in enumerate(labels)}, "C-index": c_index(y, s),
                   "Spearman rho": rho, "Spearman p": prho, **overall}]).to_csv(O1 / "rscore_heldout_by_ci_stats.csv",
                                                                                  index=False)

    # weighted kappa: SI + H-scan top 3 (overall) proportional-odds model, cohort174, kidney_pr_all, 10 x 5-fold CV;
    # predicted class = argmax of the class probabilities averaged over the 10 held-out predictions ('CV test').
    # Same source and rule as v10_figures/rscore_banff/weighted_kappa/export_weighted_kappa.py.
    pv = pd.read_excel(KAPPA_MODEL, sheet_name="R-score per visit")
    pv = pv[pv["condition"] == "kidney_pr_all"].reset_index(drop=True)
    summ = pd.read_excel(KAPPA_MODEL, sheet_name="Summary per condition").set_index("condition").loc["kidney_pr_all"]
    yk = pv["ci"].to_numpy(int); xk = pv["predicted ci (CV)"].to_numpy(int)
    qwk = kappa(yk, xk)
    assert abs(qwk - summ["QW kappa [CV test]"]) < 1e-9 and abs(qwk - 0.7780286811929968) < 1e-9
    O = confusion(yk, xk).astype(int)
    pd.DataFrame(O, index=["true " + l for l in labels], columns=["predicted " + l for l in labels]).to_csv(
        O1 / "rscore_weighted_kappa_confusion.csv")
    pv[["patient", "scan_id", "ci", "predicted ci (CV)"]].to_csv(O1 / "rscore_weighted_kappa_per_patient.csv",
                                                                 index=False)
    rng = np.random.default_rng(42)
    bs = [kappa(yk[i], xk[i]) for i in (rng.integers(0, len(yk), len(yk)) for _ in range(2000))]
    z, p = kappa_null_test(yk, xk)
    pd.DataFrame([{"model": "SI + H-scan top 3 overall, proportional odds, kidney_pr_all, 10 x 5-fold CV",
                   "predictions": "CV test: argmax of the class probabilities averaged over the 10 held-out predictions",
                   "n": len(yk), "quadratic weighted kappa": qwk, "95% CI low": np.percentile(bs, 2.5),
                   "95% CI high": np.percentile(bs, 97.5), "z (H0 kappa = 0)": z, "p (H0 kappa = 0)": p,
                   "linear weighted kappa": kappa(yk, xk, "linear"), "unweighted kappa": kappa(yk, xk, None),
                   "exact accuracy": np.mean(yk == xk)}]).to_csv(O1 / "rscore_weighted_kappa_stats.csv", index=False)
    print(f"01 R-score n={len(y)}  rho={rho:.3f}  Welch p={overall['Welch ANOVA p']:.2e}  QWK={qwk:.3f} p={p:.2e}")


# ---------------------------------------------------------------- 02 R-score all combinations
def section02():
    O2.mkdir(parents=True, exist_ok=True)
    d = RR / "ci_R_score_all_combinations" / "results"
    files = {"spearman": "pivot_spearman_rho_0to3.csv", "qwk": "pivot_qwk_quadratic.csv",
             "auc_ci0_vs_ci12": "pivot_auc_ci0_vs_ci1-2.csv", "auc_ci0_vs_ci123": "pivot_auc_ci0_vs_ci1-3.csv"}
    labels = {"spearman": "Spearman ρ", "qwk": "Quadratic weighted κ",
              "auc_ci0_vs_ci12": "AUC ci 0 vs 1–2", "auc_ci0_vs_ci123": "AUC ci 0 vs 1–2 + HD"}
    piv = {k: pd.read_csv(d / f) for k, f in files.items()}
    ref = piv["spearman"].combination.tolist()
    assert all(p.combination.tolist() == ref for p in piv.values()) and len(ref) == 63
    g = long_grid(piv, labels)
    g.to_csv(O2 / "rscore_all_combinations_grid.csv", index=False, encoding="utf-8-sig")
    print("02 grid rows", len(g))


# ---------------------------------------------------------------- 03 eGFR slope
SLOPE4 = ["≥ −1", "−3 to −1", "−5 to −3", "< −5"]
SLOPE3 = ["≥ −1", "−3 to −1", "< −3"]


def section03():
    O3.mkdir(parents=True, exist_ok=True)
    d = RR / "egfr_slope_all_combinations"
    xl = d / "results" / "egfr_slope_all_combinations_results.xlsx"
    sheets = {"c_index": "Grid_C_index", "auc_lt1": "Grid_AUC_lt-1", "auc_lt3": "Grid_AUC_lt-3",
              "auc_lt5": "Grid_AUC_lt-5"}
    labels = {"c_index": "C-index", "auc_lt1": "AUC slope < −1", "auc_lt3": "AUC slope < −3",
              "auc_lt5": "AUC slope < −5"}
    piv = {}
    for k, s in sheets.items():
        P = pd.read_excel(xl, sheet_name=s)
        P.columns = [c.replace("_LR", "") for c in P.columns]
        piv[k] = P
    g = long_grid(piv, labels)
    g.to_csv(O3 / "egfr_all_combinations_grid.csv", index=False, encoding="utf-8-sig")
    best = g[g.metric == "c_index"].sort_values("value", ascending=False).iloc[0]
    print("03 best C-index:", best.combination, best.strategy, round(best.value, 3))

    # best model, held-out scores with / without ACR
    H = pd.read_excel(d / "addon_acr" / "hscan_esd_ss_top1_acr_results.xlsx", sheet_name="Heldout_scores")
    assert len(H) == 91
    H.to_csv(O3 / "egfr_best_model_heldout_scores.csv", index=False)
    y = H.true_class.to_numpy(int)
    st_rows, pr = [], []
    for v in ["Radiomics only", "Radiomics + ACR"]:
        x = H[f"score_{v}"].to_numpy(float)
        overall, pairs = group_tests([x[y == k] for k in range(4)], SLOPE4, adjust=False)
        pairs.insert(0, "variant", v); pr.append(pairs)
        rho, prho = stats.spearmanr(x, H.slope)
        st_rows.append({"variant": v, "n": len(x), "C-index": c_index(y, x), "Spearman rho vs slope": rho,
                        "Spearman p": prho, **overall})
    pd.concat(pr).to_csv(O3 / "egfr_best_model_pairwise.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(st_rows).to_csv(O3 / "egfr_best_model_stats.csv", index=False)
    print(pd.DataFrame(st_rows)[["variant", "C-index", "Welch ANOVA p"]].round(4).to_string(index=False))

    # same 4-group model and held-out scores, regrouped into 3 groups (-5 to -3 and < -5 combined); no refit
    H3 = H.copy()
    H3["true_class"] = np.minimum(H3.true_class.to_numpy(int), 2)
    H3.to_csv(O3 / "egfr_best_model_ordinal3_heldout_scores.csv", index=False, encoding="utf-8-sig")
    y3 = H3.true_class.to_numpy(int)
    st_rows, pr = [], []
    for v in ["Radiomics only", "Radiomics + ACR"]:
        x = H3[f"score_{v}"].to_numpy(float)
        overall, pairs = group_tests([x[y3 == k] for k in range(3)], SLOPE3, adjust=False)
        pairs.insert(0, "variant", v); pr.append(pairs)
        rho, prho = stats.spearmanr(x, H3.slope)
        st_rows.append({"variant": v, "n": len(x), "C-index (3 groups)": c_index(y3, x),
                        **{f"AUC < {c}": c_index((H3.slope < c).astype(int).to_numpy(), x) for c in (-1, -3, -5)},
                        "Spearman rho vs slope": rho, "Spearman p": prho, **overall})
    pd.concat(pr).to_csv(O3 / "egfr_best_model_ordinal3_pairwise.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(st_rows).to_csv(O3 / "egfr_best_model_ordinal3_stats.csv", index=False, encoding="utf-8-sig")
    print(pd.DataFrame(st_rows)[["variant", "C-index (3 groups)", "AUC < -1", "AUC < -3", "AUC < -5",
                                 "Welch ANOVA p"]].round(4).to_string(index=False))

    # baseline clinical markers by slope group
    C = pd.read_csv(RR / "01_eGFR_slope_model3_min1yr" / "Heat_map" / "model_cv" / "clinical_addon" /
                    "baseline_clinical_model3_91.csv")
    M = H[["patient_id", "slope", "true_class"]].merge(C, on="patient_id", how="left", validate="1:1")
    assert M.baseline_ACR_mg_per_mmol.notna().all() and len(M) == 91
    M["ACR (mg/mmol, floored 0.1)"] = M.baseline_ACR_mg_per_mmol.clip(lower=0.1)
    M["eGFR"] = M.baseline_eGFR_ckdepi2021
    M["KFRE 2-yr risk (%)"] = 100 * M.KFRE_2yr_risk
    for c in (1, 3, 5):
        M[f"binary_lt{c}"] = (M.slope < -c).astype(int)
    M[["patient_id", "slope", "true_class", "binary_lt1", "binary_lt3", "binary_lt5", "baseline_ACR_mg_per_mmol",
       "ACR (mg/mmol, floored 0.1)", "eGFR", "KFRE 2-yr risk (%)"]].to_csv(
        O3 / "egfr_clinical_markers_per_patient.csv", index=False)
    markers = [("ACR", "ACR (mg/mmol, floored 0.1)", np.log10), ("eGFR", "eGFR", lambda v: v),
               ("KFRE", "KFRE 2-yr risk (%)", np.log10)]
    groupings = {"ordinal4": ("true_class", SLOPE4)}
    for c in (1, 3, 5):
        groupings[f"binary_lt{c}"] = (f"binary_lt{c}", [f"≥ −{c}", f"< −{c}"])
    pr, st_rows = [], []
    for gk, (col, labs) in groupings.items():
        gy = M[col].to_numpy(int)
        for mk, mc, tf in markers:
            v = tf(M[mc].to_numpy(float))
            overall, pairs = group_tests([v[gy == k] for k in range(len(labs))], labs, adjust=False)
            pairs.insert(0, "marker", mk); pairs.insert(0, "grouping", gk); pr.append(pairs)
            st_rows.append({"grouping": gk, "marker": mk, "scale tested": "log10" if mk != "eGFR" else "linear",
                            **{f"n {l}": int((gy == k).sum()) for k, l in enumerate(labs)}, **overall})
    pd.concat(pr).to_csv(O3 / "egfr_clinical_markers_pairwise.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(st_rows).to_csv(O3 / "egfr_clinical_markers_stats.csv", index=False, encoding="utf-8-sig")
    print(pd.concat(pr).query("p < 0.05")[["grouping", "marker", "group A", "group B", "p"]].to_string(index=False))


# ---------------------------------------------------------------- 04 / 05 spatial
ROIS = ["all_images", "matched", "opposite"]
ROI_NAME = {"All images": "all_images", "Matched third": "matched", "Opposite third": "opposite"}


def xpos(roi, j, W, width=0.8):
    bw = width / W
    return ROIS.index(roi) + 1 + (j - (W - 1) / 2) * bw


def bar_tables(auc, windows, model_col, win_col, roi_tests, win_tests, out_prefix, folder):
    """auc: rows model / window / spatial / AUC / CI. Writes <prefix>_auc.csv and <prefix>_tests.csv (with x)."""
    W = len(windows)
    A = auc[auc[model_col].isin(["frozen", "polr"]) & auc[win_col].isin(windows)].copy()
    A["window index"] = A[win_col].map({w: i + 1 for i, w in enumerate(windows)})
    A["x"] = [xpos(r, j - 1, W) for r, j in zip(A.spatial, A["window index"])]
    A["CI_high"] = A.CI_high.clip(upper=1)
    A.to_csv(folder / f"{out_prefix}_auc.csv", index=False)
    T = []
    for r in roi_tests.itertuples():
        a, b = r.comparison.split(" vs ")
        j = windows.index(r.window)
        T.append({"model": r.model, "type": "ROI (paired DeLong)", "group A": f"{a}, {r.window}",
                  "group B": f"{b}, {r.window}", "x A": xpos(ROI_NAME[a], j, W), "x B": xpos(ROI_NAME[b], j, W),
                  "AUC A": r.AUC_1, "AUC B": r.AUC_2, "p 1-sided": r.p_1, "p 2-sided": r.p_2})
    for r in (win_tests.itertuples() if win_tests is not None else []):
        ja, jb = windows.index(r.w1), windows.index(r.w2)
        T.append({"model": r.model, "type": getattr(r, "test_type", "window (unpaired DeLong)"),
                  "group A": f"{r.spatial}, {r.w1}",
                  "group B": f"{r.spatial}, {r.w2}", "x A": xpos(r.spatial, ja, W), "x B": xpos(r.spatial, jb, W),
                  "AUC A": r.AUC_1, "AUC B": r.AUC_2, "p 1-sided": r.p_1, "p 2-sided": r.p_2})
    T = pd.DataFrame(T)
    T = T[T.model.isin(["frozen", "polr"])]
    T.to_csv(folder / f"{out_prefix}_tests.csv", index=False)
    print(out_prefix, "significant (1-sided):")
    print(T[T["p 1-sided"] < 0.05][["model", "group A", "group B", "p 1-sided"]].to_string(index=False))


def subroi_pairs(P, windows, prefix, folder):
    """Every pair among the matched / opposite bars (all-images bars excluded). P: model, w1, roi1, w2, roi2, test,
    p_1 (one-sided, H1 first > second, as reported), p_2, status. Untestable pairs (nested windows) keep p = NaN."""
    W = len(windows)
    P = P[P.model.isin(["frozen", "polr"])].copy()
    P["x A"] = [xpos(r, windows.index(w), W) for r, w in zip(P.roi1, P.w1)]
    P["x B"] = [xpos(r, windows.index(w), W) for r, w in zip(P.roi2, P.w2)]
    P["group A"] = P.roi1 + ", " + P.w1
    P["group B"] = P.roi2 + ", " + P.w2
    P = P.rename(columns={"p_1": "p 1-sided", "p_2": "p 2-sided"})
    P[["model", "group A", "group B", "test", "x A", "x B", "p 1-sided", "p 2-sided", "status"]].to_csv(
        folder / f"{prefix}_subroi_pairs.csv", index=False)
    print(prefix, "sub-ROI pairs:", len(P), "| testable:", int(P["p 1-sided"].notna().sum()),
          "| p < 0.05:", int((P["p 1-sided"] < 0.05).sum()))


def section04():
    O4.mkdir(parents=True, exist_ok=True)
    x = RR / "Final_100" / "spatial_Final_100.xlsx"
    A = pd.read_excel(x, sheet_name="AUC").rename(columns={"half": "window"})
    R = pd.read_excel(x, sheet_name="ROI comparisons DeLong").rename(
        columns={"half": "window", "p_DeLong_1sided": "p_1", "p_DeLong_2sided": "p_2"})
    E = pd.read_excel(x, sheet_name="Early vs late DeLong").rename(
        columns={"AUC_early": "AUC_1", "AUC_late": "AUC_2", "p_DeLong_1sided": "p_1", "p_DeLong_2sided": "p_2"})
    E["w1"], E["w2"] = "early", "late"
    R = R[R.window.isin(["early", "late"])]
    bar_tables(A, ["early", "late"], "model", "window", R, E, "auc_by_roi_early_late", O4)
    # all 6 pairs among matched / opposite x early / late (four_group_onesided.xlsx; codes ME, OE, ML, OL)
    G = pd.read_excel(RR / "Final_100" / "four_group_onesided.xlsx", sheet_name="Pairwise one-sided")
    dec = lambda c: ({"M": "matched", "O": "opposite"}[c[0]], {"E": "early", "L": "late"}[c[1]])
    rows = []
    for r in G.itertuples():
        a, b = r.comparison.split(" vs ")
        (ra, wa), (rb, wb) = dec(a), dec(b)
        rows.append({"model": r.model, "w1": wa, "roi1": ra, "w2": wb, "roi2": rb, "test": r.test,
                     "p_1": r.p_1sided, "p_2": r.p_2sided, "status": r.H1})
    subroi_pairs(pd.DataFrame(rows), ["early", "late"], "auc_by_roi_early_late", O4)
    w = A[A.window != "whole"].drop_duplicates("window")[["window", "half_label", "n", "n_ci12"]]
    w.to_csv(O4 / "windows.csv", index=False)
    mg = RR / "Final_100" / "cross_analysis_median_gap" / "spatial_cross_Final_100_median.xlsx"
    cross_tables(mg, "Median gap", "cross_roi_median_gap", O4)
    cross_tables(mg, "Median gap", "cross_time_median_gap", O4, sheet="Cross time AUC", fit="Disjoint full fit")
    cross_median_gap_tests(mg, O4)


def cross_tables(xlsx, window_set, prefix, folder, sheet="Cross ROI AUC", fit="Whole-cohort CV"):
    """Train ROI x test ROI AUC per test window. Cross ROI: models fit by whole-cohort CV, scored within each
    window. Cross time: fit on one disjoint window, tested on the other (fit_design 'Disjoint full fit')."""
    c = pd.read_excel(xlsx, sheet_name=sheet)
    c = c[(c.window_set == window_set) & (c.fit_design == fit)]
    if sheet == "Cross ROI AUC":
        assert set(c.train_window) == {"whole"}
    c.to_csv(folder / f"{prefix}.csv", index=False)
    print(prefix, len(c), "rows; test windows", list(c.test_window.unique()))


def cross_median_gap_tests(xlsx, folder):
    """Test one-factor AUC contrasts among the ordinal logistic bars using saved scores."""
    sys.path.insert(0, str(REPO / "Shared Scripts" / "pocus_kidney_main_figures"))
    try:
        from export_main_figure_data import delong_auc_var, paired_delong
    finally:
        sys.path.pop(0)

    bounds = pd.read_excel(xlsx, sheet_name="Windows").set_index("group")
    window_order = ["early", "late"]
    factors = ["test_window", "train_ROI", "test_ROI"]
    roi_order = {roi: i for i, roi in enumerate(ROIS)}
    designs = [
        ("cross_roi_median_gap", "Whole cohort scores", "Test window"),
        ("cross_time_median_gap", "Median gap transfer scores", "Transfer direction"),
    ]
    results = []
    for prefix, score_sheet, time_factor in designs:
        scores = pd.read_excel(xlsx, sheet_name=score_sheet).set_index("scan_id", verify_integrity=True)
        cohorts = {
            window: scores[scores.gap.between(bounds.loc[window, "gap_min_days"],
                                              bounds.loc[window, "gap_max_days"])]
            for window in window_order
        }
        assert cohorts["early"].index.intersection(cohorts["late"].index).empty
        auc_table = pd.read_csv(folder / f"{prefix}.csv")
        auc_table = auc_table[auc_table.model.eq("polr")].sort_values(
            factors, key=lambda col: col.map(roi_order) if col.name != "test_window" else col)
        assert len(auc_table) == 18 and not auc_table.duplicated(factors).any()

        conditions = []
        for row in auc_table.itertuples(index=False):
            cohort = cohorts[row.test_window]
            fit = "whole" if row.train_window == "whole" else f"{row.train_window} to {row.test_window}"
            predictions = cohort[f"{fit}/polr/{row.train_ROI}/{row.test_ROI}"].to_numpy(float)
            outcomes = cohort.ci_grade.isin([1, 2]).to_numpy()
            assert np.isfinite(predictions).all() and cohort.ci_grade.isin([0, 1, 2]).all()
            assert len(cohort) == row.n_test and outcomes.sum() == row.n_test_ci12
            auc, variance = delong_auc_var(outcomes, predictions)
            assert np.isclose(auc, row.AUC, rtol=0, atol=1e-12), f"AUC mismatch for {row}"
            conditions.append((row, outcomes, predictions, auc, variance))

        for (a, y_a, score_a, auc_a, var_a), (b, _, score_b, auc_b, var_b) in combinations(conditions, 2):
            changes = [factor for factor in factors if getattr(a, factor) != getattr(b, factor)]
            if len(changes) != 1:
                continue
            factor = changes[0]
            if factor == "test_window":
                z = (auc_a - auc_b) / np.sqrt(var_a + var_b)
                test = "Unpaired DeLong"
            else:
                _, _, z = paired_delong(y_a, score_a, score_b)
                test = "Paired DeLong"
            assert np.isfinite(z), f"Undefined test statistic for {a} and {b}"
            results.append({
                "fit_design": a.fit_design, "model": a.model,
                "factor": {"test_window": time_factor, "train_ROI": "Training ROI", "test_ROI": "Testing ROI"}[factor],
                "test": test,
                **{f"{name}_{side}": getattr(row, name) for side, row in [("A", a), ("B", b)]
                   for name in ["train_window", "test_window", "train_ROI", "test_ROI", "n_test"]},
                "AUC_A": auc_a, "AUC_B": auc_b, "AUC_difference_A_minus_B": auc_a - auc_b,
                "z": z, "p_two_sided": 2 * stats.norm.sf(abs(z)),
            })

    tests = pd.DataFrame(results)
    counts = tests.groupby(["fit_design", "factor"]).size()
    assert sorted(counts) == [9, 9, 18, 18, 18, 18]
    order = np.argsort(tests.p_two_sided.to_numpy())
    scaled = tests.p_two_sided.to_numpy()[order] * np.arange(len(tests), 0, -1)
    tests.loc[order, "p_Holm"] = np.minimum(1, np.maximum.accumulate(scaled))
    tests["significant_Holm_0_05"] = tests.p_Holm.lt(0.05)
    tests["significance"] = np.select(
        [tests.p_Holm.lt(0.001), tests.p_Holm.lt(0.01), tests.significant_Holm_0_05],
        ["***", "**", "*"], default="ns",
    )
    tests.to_csv(folder / "cross_roi_time_median_gap_pairwise_tests.csv", index=False)
    print(f"Cross-ROI and cross-time pairwise tests: {tests.significant_Holm_0_05.sum()} / {len(tests)} significant after Holm adjustment")


def all_images_cumulative_tests(scores, auc, windows):
    """Bootstrap fixed-score AUC differences, sharing each patient's draw across nested windows."""
    masks = np.stack([scores[f"in_{window}"].to_numpy(bool) for window in windows])
    assert scores.scan_id.is_unique and scores.ci_grade.isin([0, 1, 2]).all()
    assert np.all(~masks[:-1] | masks[1:])
    outcome = scores.ci_grade.isin([1, 2]).to_numpy()
    membership = masks.sum(axis=0)
    rng = np.random.default_rng(42)
    n_boot = 20_000
    weights = np.zeros((n_boot, len(scores)))
    for count in range(1, len(windows) + 1):
        for positive in (False, True):
            idx = np.flatnonzero((membership == count) & (outcome == positive))
            if idx.size:
                weights[:, idx] = rng.multinomial(idx.size, np.full(idx.size, 1 / idx.size), size=n_boot)

    rows = []
    for model in ("frozen", "polr"):
        estimates = []
        replicates = []
        for window, mask in zip(windows, masks):
            column = "frozen__all_images" if model == "frozen" else f"polr_within_{window}__all_images"
            values = scores[column].to_numpy(float)
            pos, neg = mask & outcome, mask & ~outcome
            diff = values[pos, None] - values[None, neg]
            assert np.isfinite(diff).all()
            wins = (diff > 0) + 0.5 * (diff == 0)
            estimate = wins.mean()
            reference = auc.loc[auc.model.eq(model) & auc.window.eq(window) & auc.spatial.eq("all_images"), "AUC"]
            assert len(reference) == 1 and np.isclose(estimate, reference.iloc[0], rtol=0, atol=1e-12)
            estimates.append(estimate)
            # Stratification fixes both denominators and each patient's multiplicity is reused across windows.
            replicates.append(((weights[:, pos] @ wins) * weights[:, neg]).sum(axis=1) / wins.size)

        for i, j in combinations(range(len(windows)), 2):
            se = np.std(replicates[i] - replicates[j], ddof=1)
            z = (estimates[i] - estimates[j]) / se
            rows.append({"model": model, "spatial": "all_images", "w1": windows[i], "w2": windows[j],
                         "AUC_1": estimates[i], "AUC_2": estimates[j], "p_1": stats.norm.sf(z),
                         "p_2": 2 * stats.norm.sf(abs(z)), "test_type": "window (overlap-aware stratified bootstrap)"})
    return pd.DataFrame(rows)


def section05():
    O5.mkdir(parents=True, exist_ok=True)
    d = RR / "Final_100" / "ci_spatial_0_1_3_5_months"
    x = d / "spatial_time_windows_Final_100.xlsx"
    A = pd.read_excel(x, sheet_name="AUC")
    R = pd.read_excel(x, sheet_name="ROI comparisons DeLong").rename(
        columns={"p_DeLong_1sided": "p_1", "p_DeLong_2sided": "p_2"})
    V = pd.read_excel(x, sheet_name="Interval vs interval DeLong").rename(
        columns={"window_earlier": "w1", "window_later": "w2", "AUC_earlier": "AUC_1", "AUC_later": "AUC_2",
                 "p_DeLong_1sided": "p_1", "p_DeLong_2sided": "p_2"})
    pd.read_excel(x, sheet_name="Windows").to_csv(O5 / "windows.csv", index=False)
    iv = ["day0", ">0-1m", ">1-3m", ">3-5m"]; cu = ["0-1m", "0-3m", "0-5m"]
    bar_tables(A[A.set == "interval"], iv, "model", "window", R[R.set == "interval"], V, "auc_by_roi_interval", O5)
    cumulative = A[A.set == "cumulative"]
    cumulative_tests = all_images_cumulative_tests(pd.read_excel(x, sheet_name="Patient scores"), cumulative, cu)
    bar_tables(cumulative, cu, "model", "window", R[R.set == "cumulative"], cumulative_tests,
               "auc_by_roi_cumulative", O5)
    cw = d / "cross_analysis" / "spatial_cross_Final_100_windows.xlsx"
    # all pairs among matched / opposite bars (ROI by time DeLong: paired within window, unpaired between disjoint
    # intervals; nested cumulative windows are not testable and keep p = NaN)
    RT = pd.read_excel(cw, sheet_name="ROI by time DeLong").rename(
        columns={"window_1": "w1", "ROI_1": "roi1", "window_2": "w2", "ROI_2": "roi2", "p_1sided": "p_1",
                 "p_2sided": "p_2"})
    subroi_pairs(RT[RT.window_set == "Interval"], iv, "auc_by_roi_interval", O5)
    subroi_pairs(RT[RT.window_set == "Cumulative"], cu, "auc_by_roi_cumulative", O5)
    cross_tables(cw, "Interval", "cross_roi_interval", O5)
    cross_tables(cw, "Cumulative", "cross_roi_cumulative", O5)


if __name__ == "__main__":
    pd.set_option("display.width", 220)
    section01(); section02(); section03(); section04(); section05()

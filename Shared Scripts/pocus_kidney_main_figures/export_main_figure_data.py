"""Data + statistics for the POCUS kidney main figures (drawn by make_main_figures.m).

    python export_main_figure_data.py

Source: Results_Rough/POCUS_Kidney_Summary_2026-11-24 (01 eGFR slope, 02 ci R-score, 03 ci spatial).
Output: Results_Rough/pocus_kidney_main_figures/<01_eGFR_slope | 02_ci_Rscore | 03_ci_spatial>/*.csv

01  heat map    : held-out C-index of the ordinal (proportional-odds) models, feature set x top N, 4- and 3-group,
                  LOOCV and 5-fold (model_grid_results.xlsx)
    box plot    : LOOCV held-out linear predictor of the best C-index model (4-group, SS top 3), by slope group;
                  Welch ANOVA + Tukey-Kramer (style guide); Kruskal-Wallis / Mann-Whitney-Holm kept in the CSV
02  R-score     : held-out R-score = expected ci class, mean of 10 outer-test predictions per patient (n = 174)
    box plot    : Welch ANOVA + Tukey-Kramer; Spearman rho
    features    : selection frequency over the 50 outer folds (top 3)
    kappa       : predicted class = argmax of the 10-repeat mean class probabilities; quadratic-weighted kappa,
                  p = z-test of H0 kappa = 0 (Fleiss, Cohen & Everitt 1969), 95% CI = bootstrap (2000, seed 42)
    ROC         : 4 binary tasks of the held-out R-score, AUC with DeLong 95% CI
03  six groups  : AUC (ci 0 vs 1-2) of all images / matched / opposite third x early / late half, frozen and polr;
                  DeLong (paired within half, unpaired across halves)
    cross-cond. : AUC train ROI x test ROI, models fit on all 98 (cross_condition_auc.xlsx, design "all")
"""
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

REPO = Path(__file__).resolve().parents[2]
SRC = REPO / "Results_Rough" / "POCUS_Kidney_Summary_2026-11-24"
OUT = REPO / "Results_Rough" / "pocus_kidney_main_figures"
S01, S02, S03 = (SRC / "01_eGFR_slope_model3_min1yr", SRC / "02_ci_Rscore_radiomics_top3overall_SI_HScan",
                 SRC / "03_ci_spatial_kp_Final_98")
O01, O02, O03 = OUT / "01_eGFR_slope", OUT / "02_ci_Rscore", OUT / "03_ci_spatial"


# ---------------------------------------------------------------- statistics
def welch_anova(groups):
    g = [np.asarray(x, float) for x in groups]
    k = len(g); n = np.array([len(x) for x in g]); m = np.array([x.mean() for x in g])
    w = n / np.array([x.var(ddof=1) for x in g]); W = w.sum(); mw = (w * m).sum() / W
    lam = ((1 - w / W) ** 2 / (n - 1)).sum()
    F = ((w * (m - mw) ** 2).sum() / (k - 1)) / (1 + 2 * (k - 2) / (k ** 2 - 1) * lam)
    df2 = (k ** 2 - 1) / (3 * lam)
    return F, k - 1, df2, stats.f.sf(F, k - 1, df2)


def group_tests(groups, labels):
    F, d1, d2, p = welch_anova(groups)
    tk = stats.tukey_hsd(*groups).pvalue
    pairs = list(combinations(range(len(groups)), 2))
    mw = [stats.mannwhitneyu(groups[a], groups[b]).pvalue for a, b in pairs]
    order = np.argsort(mw); m = len(mw); holm = np.empty(m); run = 0.0
    for r, i in enumerate(order):
        run = max(run, (m - r) * mw[i]); holm[i] = min(1.0, run)
    rows = [{"group A": labels[a], "group B": labels[b], "pos A": a + 1, "pos B": b + 1,
             "n A": len(groups[a]), "n B": len(groups[b]),
             "median A": np.median(groups[a]), "median B": np.median(groups[b]),
             "Tukey-Kramer p": tk[a, b], "Mann-Whitney p": mw[i], "Mann-Whitney Holm p": holm[i]}
            for i, (a, b) in enumerate(pairs)]
    kw = stats.kruskal(*groups)
    overall = {"Welch ANOVA F": F, "df1": d1, "df2": d2, "Welch ANOVA p": p,
               "Kruskal-Wallis H": kw.statistic, "Kruskal-Wallis p": kw.pvalue}
    return overall, pd.DataFrame(rows)


def c_index(y, s):
    num = den = 0.0
    for lo in np.unique(y):
        for hi in np.unique(y):
            if hi > lo:
                d = s[y == hi][:, None] - s[y == lo][None, :]
                num += (d > 0).sum() + 0.5 * (d == 0).sum(); den += d.size
    return num / den


def _midrank(x):
    return stats.rankdata(x)


def delong_components(y, s):
    """AUC and DeLong structural components (Sun & Xu 2014)."""
    y = np.asarray(y).astype(bool); s = np.asarray(s, float)
    X, Y = s[y], s[~y]; m, n = len(X), len(Y)
    tx, ty, tz = _midrank(X), _midrank(Y), _midrank(np.r_[X, Y])
    auc = (tz[:m].sum() - m * (m + 1) / 2) / (m * n)
    v10 = (tz[:m] - tx) / n
    v01 = 1 - (tz[m:] - ty) / m
    return auc, v10, v01


def delong_auc_var(y, s):
    auc, v10, v01 = delong_components(y, s)
    return auc, np.var(v10, ddof=1) / len(v10) + np.var(v01, ddof=1) / len(v01)


def paired_delong(y, s1, s2):
    a1, x1, y1 = delong_components(y, s1)
    a2, x2, y2 = delong_components(y, s2)
    sx = np.cov(np.vstack([x1, x2])); sy = np.cov(np.vstack([y1, y2]))
    S = sx / len(x1) + sy / len(y1)
    var = S[0, 0] + S[1, 1] - 2 * S[0, 1]
    z = (a1 - a2) / np.sqrt(var)
    return a1, a2, z


# ---------------------------------------------------------------- 01 eGFR slope
def section01():
    O01.mkdir(parents=True, exist_ok=True)
    mc = S01 / "Heat_map" / "model_cv"
    summ = pd.read_excel(mc / "model_grid_results.xlsx", sheet_name="Summary")
    o = summ[summ.task.str.startswith("Ordinal")].copy()
    o["groups"] = o.task.str.extract(r"(\d)-cat").astype(int)
    o[["feature_set", "top_n", "groups", "cv", "main", "selected features (fraction of folds)"]].rename(
        columns={"main": "C-index"}).to_csv(O01 / "egfr_ordinal_cindex_grid.csv", index=False)

    best = o[(o.groups == 4) & (o.cv == "LOOCV")].sort_values("main", ascending=False).iloc[0]
    run = f"{best.feature_set}|top{int(best.top_n)}|{best.task}|LOOCV"
    pred = pd.read_csv(mc / "all_heldout_predictions.csv")
    d = pred[pred.run == run].copy()
    assert len(d) == 91
    x, g = d.heldout_score.to_numpy(float), d.true_class.to_numpy(int)
    assert abs(c_index(g, x) - best.main) < 1e-9
    labels = ["≥ −1", "−3 to −1", "−5 to −3", "< −5"]
    groups = [x[g == i] for i in range(4)]
    overall, pairs = group_tests(groups, labels)
    rho, prho = stats.spearmanr(x, g)
    d.rename(columns={"true_class": "slope group (0 = ≥-1 ... 3 = <-5)", "heldout_score": "held-out score"}).to_csv(
        O01 / "egfr_best_cindex_heldout_scores.csv", index=False)
    pairs.to_csv(O01 / "egfr_best_cindex_pairwise.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame([{"model": f"Ordinal 4-group, {best.feature_set} top {int(best.top_n)}, LOOCV", "n": len(x),
                   **{f"n {l}": len(v) for l, v in zip(labels, groups)}, "C-index": best.main,
                   "Spearman rho": rho, "Spearman p": prho, **overall,
                   "selected features (fraction of folds)": best["selected features (fraction of folds)"]}]
                 ).to_csv(O01 / "egfr_best_cindex_stats.csv", index=False, encoding="utf-8-sig")
    print("01 best:", run, round(best.main, 3), "Welch p", overall["Welch ANOVA p"])
    print(pairs[["group A", "group B", "Tukey-Kramer p", "Mann-Whitney Holm p"]].round(4).to_string(index=False))


# ---------------------------------------------------------------- 02 ci R-score
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


def section02():
    O02.mkdir(parents=True, exist_ok=True)
    pm = pd.read_csv(S02 / "outer_predictions_per_sample_mean.csv")
    assert len(pm) == 174
    y = pm.true_index.to_numpy(int); s = pm.r_score.to_numpy(float); yhat = pm.predicted_index.to_numpy(int)
    assert (pm[[f"probability_{k}" for k in range(4)]].to_numpy().argmax(1) == yhat).all()
    pm.to_csv(O02 / "rscore_heldout_per_patient.csv", index=False)

    # box plot stats
    labels = ["ci 0", "ci 1", "ci 2", "HD"]
    overall, pairs = group_tests([s[y == k] for k in range(4)], labels)
    rho, prho = stats.spearmanr(s, y)
    ci = c_index(y, s)
    pairs.to_csv(O02 / "rscore_by_ci_pairwise.csv", index=False)
    pd.DataFrame([{"score": "held-out R-score (expected ci class, mean of 10 outer-test predictions)", "n": len(y),
                   **{f"n {l}": int((y == k).sum()) for k, l in enumerate(labels)}, "C-index": ci,
                   "Spearman rho": rho, "Spearman p": prho, **overall}]).to_csv(O02 / "rscore_by_ci_stats.csv",
                                                                                  index=False)
    print("02 box: C", round(ci, 3), "rho", round(rho, 3), "Welch p", overall["Welch ANOVA p"])

    # feature selection frequency
    co = pd.read_csv(S02 / "outer_fold_coefficients.csv")
    nf = co.outer_fold.nunique()
    fr = co.groupby("feature").outer_fold.nunique().sort_values(ascending=False)
    nice = {"HScan__firstorder_90Percentile": "H-scan 90th percentile",
            "SI__firstorder_90Percentile": "SI 90th percentile", "SI__firstorder_10Percentile": "SI 10th percentile"}
    sel = pd.DataFrame({"feature": fr.index, "label": [nice.get(f, f) for f in fr.index], "folds selected": fr.values,
                        "outer folds": nf, "percent": 100 * fr.values / nf})
    sel.to_csv(O02 / "rscore_feature_selection_frequency.csv", index=False)
    print(sel.head(5).to_string(index=False))

    # weighted kappa
    O = confusion(y, yhat).astype(int)
    cm = pd.DataFrame(O, index=["true " + l for l in labels], columns=["predicted " + l for l in labels])
    cm.to_csv(O02 / "rscore_weighted_kappa_confusion.csv")
    rng = np.random.default_rng(42)
    bs = [kappa(y[i], yhat[i]) for i in (rng.integers(0, len(y), len(y)) for _ in range(2000))]
    z, p = kappa_null_test(y, yhat)
    qwk = kappa(y, yhat)
    assert abs(qwk - 0.7395269570023721) < 1e-9
    pd.DataFrame([{"predictions": "argmax of the 10-repeat mean class probabilities (held-out)", "n": len(y),
                   "quadratic weighted kappa": qwk, "95% CI low": np.percentile(bs, 2.5),
                   "95% CI high": np.percentile(bs, 97.5), "z (H0 kappa = 0)": z, "p (H0 kappa = 0)": p,
                   "linear weighted kappa": kappa(y, yhat, "linear"), "unweighted kappa": kappa(y, yhat, None),
                   "exact accuracy": np.mean(y == yhat)}]).to_csv(O02 / "rscore_weighted_kappa_stats.csv", index=False)
    print("02 kappa", round(qwk, 3), "p", p)

    # ROC tasks
    tasks = [("ci0_vs_ci12", {0}, {1, 2}, "r_score_auc_class_0_vs_1_to_2_excluding_3"),
             ("ci0_vs_ci13", {0}, {1, 2, 3}, "r_score_auc_class_0_to_0_vs_1_to_3"),
             ("ci01_vs_ci23", {0, 1}, {2, 3}, "r_score_auc_class_0_to_1_vs_2_to_3"),
             ("ci012_vs_ci3", {0, 1, 2}, {3}, "r_score_auc_class_0_to_2_vs_3_to_3")]
    import json
    ref = json.loads((S02 / "testing_metrics_per_sample_mean.json").read_text())
    rows, auc_rows = [], []
    for key, neg, pos, rk in tasks:
        keep = np.isin(y, list(neg | pos)); yy = np.isin(y[keep], list(pos)).astype(int); ss = s[keep]
        auc, var = delong_auc_var(yy, ss)
        assert abs(auc - ref[rk]) < 1e-9, key
        se = np.sqrt(var)
        auc_rows.append({"Task key": key, "n negative": int((yy == 0).sum()), "n positive": int(yy.sum()), "AUC": auc,
                         "CI low": max(0, auc - 1.96 * se), "CI high": min(1, auc + 1.96 * se)})
        rows.append(pd.DataFrame({"Task key": key, "sample_id": pm.sample_id[keep].to_numpy(), "ci (HD = 3)": y[keep],
                                  "y": yy, "R-score": ss}))
    pd.concat(rows).to_csv(O02 / "rscore_roc_scores.csv", index=False)
    pd.DataFrame(auc_rows).to_csv(O02 / "rscore_roc_auc.csv", index=False)
    print(pd.DataFrame(auc_rows).round(3).to_string(index=False))


# ---------------------------------------------------------------- 03 ci spatial
def section03():
    O03.mkdir(parents=True, exist_ok=True)
    ps = pd.read_excel(S03 / "spatial_Final_98.xlsx", sheet_name="Patient scores")
    y = (ps.ci_grade >= 1).to_numpy(int); half = ps.half.to_numpy()
    rois = ["all_images", "matched", "opposite"]
    groups = [(f"{r[0].upper()}{h[0].upper()}", r, h) for h in ("early", "late") for r in rois]
    models = {"frozen": "frozen__{}", "polr": "polr_withinhalf__{}"}
    # pre-specified one-sided directions (matched > opposite, early > late); comparisons with all images two-sided
    auc_rows, test_rows = [], []
    for m, col in models.items():
        sc = {r: ps[col.format(r)].to_numpy(float) for r in rois}
        av = {}
        for g, r, h in groups:
            k = half == h
            a, v = delong_auc_var(y[k], sc[r][k]); av[g] = (a, v)
            auc_rows.append({"model": m, "group": g, "ROI": r, "half": h, "n": int(k.sum()), "n_ci12": int(y[k].sum()),
                             "AUC": a, "CI_low": a - 1.96 * np.sqrt(v), "CI_high": min(1, a + 1.96 * np.sqrt(v))})
        for (g1, r1, h1), (g2, r2, h2) in combinations(groups, 2):
            if h1 == h2:
                k = half == h1
                a1, a2, z = paired_delong(y[k], sc[r1][k], sc[r2][k]); test = "paired DeLong"
            else:
                (a1, v1), (a2, v2) = av[g1], av[g2]; z = (a1 - a2) / np.sqrt(v1 + v2); test = "unpaired DeLong"
            prespec = "all_images" not in (r1, r2)   # matched / opposite x early / late: original hypotheses
            # all-images comparisons: only same half (vs each third) or same ROI (early vs late)
            if not prespec and h1 != h2 and r1 != r2:
                continue
            p1 = stats.norm.sf(z)                     # H1 first > second (pairs are ordered matched->opposite, early->late)
            p2 = 2 * stats.norm.sf(abs(z))
            test_rows.append({"model": m, "group 1": g1, "group 2": g2, "test": test, "AUC 1": a1, "AUC 2": a2,
                              "z": z, "p 1-sided (H1: 1 > 2)": p1, "p 2-sided": p2,
                              "p reported": p1 if prespec else p2,
                              "sidedness reported": "1-sided (pre-specified)" if prespec else "2-sided"})
    A = pd.DataFrame(auc_rows); T = pd.DataFrame(test_rows)
    # check against the original four-group results
    ref = pd.read_excel(S03 / "four_group_onesided.xlsx", sheet_name="Pairwise one-sided")
    for r in ref.itertuples():
        g1, g2 = r.comparison.split(" vs ")
        t = T[(T.model == r.model) & (T["group 1"] == g1) & (T["group 2"] == g2)].iloc[0]
        assert abs(t["p 1-sided (H1: 1 > 2)"] - r.p_1sided) < 1e-6, (r.model, r.comparison)
    A.to_csv(O03 / "spatial_six_group_auc.csv", index=False)
    T.to_csv(O03 / "spatial_six_group_tests.csv", index=False)
    print(A.round(3).to_string(index=False))
    print(T[T["p reported"] < 0.05][["model", "group 1", "group 2", "p reported", "sidedness reported"]].to_string())

    cc = pd.read_excel(S03 / "cross_condition_auc.xlsx", sheet_name="AUC with CI")
    cc = cc[cc.design == "all"]
    cc.to_csv(O03 / "spatial_cross_condition_all_auc.csv", index=False)
    print("cross-condition rows", len(cc))


if __name__ == "__main__":
    pd.set_option("display.width", 220)
    section01(); section02(); section03()

"""Ordinal (4-category) top-N radiomics models, evaluated for eGFR slope < -3 only, >= 1.5 yr expanded-RF cohort.

Mirrors Results_Rough/01_eGFR_slope_model3_min1yr/Heat_map/model_cv/run_model_grid.py ('Ordinal 4-cat' task):
  Cohort    case8_min1p5yr_expanded_rf_slope_groups.xlsx, patients with a first-visit kidney ROI (allowed RF size).
            Slope = that cohort's own LMM fit (>= 1.5 yr cohort, 139 patients).
  Target    4 groups: slope >= -1 | -3 <= s < -1 | -5 <= s < -3 | s < -5 (0..3).
  Features  same as 01 radiomics_slope_category_heatmap.load(): kidney minus pelvis, mean over the first-visit
            kidney-ROI images; SI bw 9, H-scan 3, ESD 1.9, SS 1.9; Original, LoG sigma 1/2/3, Wavelet LL/LH/HL/HH;
            first-order + GLCM/GLRLM/GLSZM/GLDM/NGTDM; features with < 90 % non-missing or zero SD dropped.
            Per-image pickles from the shared 13_multimap_filters_binwidth cache (read only) plus
            multimap_cache_extra/ (extract_multimap_missing.py, same extract_one; verified identical on 3 cached images).
  Pipeline  inside each training fold (nested_loocv_radiomics.Preprocessor / select_top_spearman_per_modality):
            median imputation -> Spearman-priority |Pearson r| >= 0.9 removal -> top-N by |Spearman rho| with the
            4-cat group -> standardise -> proportional-odds ordinal logistic (statsmodels OrderedModel, logit, BFGS).
  Sets      ESD, SI, SS, H-scan alone; ESD + SS, SS + H-scan, ESD + H-scan, ESD + SS + H-scan with top-N chosen
            per modality (N from each modality) or overall (N pooled across the set); top N = 1, 2, 3, 5.
  Score     held-out linear predictor x.beta (no thresholds), higher = steeper decline.
  Reported  AUC of the score for slope < -3 (main), 95 % CI (2000 stratified bootstrap resamples of the held-out
            scores), and the 4-cat C-index for reference.
  CV        LOOCV (main, as 01) and stratified 5-fold x 10 repeats (seeds 2026 + r; mean of per-fold AUCs) as a check.
Run: C:/Users/tupil/venvs/radiomics_env/python.exe ordinal4_topN_lt3_min1p5yr.py
"""
import pickle
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from joblib import Parallel, delayed
from sklearn.metrics import roc_auc_score, roc_curve
from sklearn.model_selection import LeaveOneOut, StratifiedKFold

ROOT = Path(__file__).resolve().parents[2]
RES = ROOT / "Results_Rough" / "case8_min1p5yr_expanded_rf"
MODEL_CODE = ROOT / "Results_Rough" / "02_ci_Rscore_radiomics_top3overall_SI_HScan"
SHARED = Path(r"\\Tupil\tupil data\Hikaru\POCUS_Kidney_Analysis\Analysis_scripts\Redo_analysis_1106_2026"
              r"\v10_figures\egfr_slope\13_multimap_filters_binwidth\cache")
EXTRA = RES / "multimap_cache_extra"
sys.path.insert(0, str(Path(__file__).resolve().parent))
from extract_multimap_missing import first_visit_images, COHORT  # noqa: E402

BIN_WIDTHS = {"SI": "9", "HScan": "3", "SS": "1.9", "ESD": "1.9"}
FILTERS = {"original", "log-sigma-1-0-mm-3D", "log-sigma-2-0-mm-3D", "log-sigma-3-0-mm-3D",
           "wavelet-LL", "wavelet-LH", "wavelet-HL", "wavelet-HH"}
CLASSES = {"firstorder", "glcm", "glrlm", "glszm", "gldm", "ngtdm"}
LBL = {"ESD": "ESD", "SI": "SI", "SS": "SS", "HScan": "H-scan"}
SETS = [(["ESD"], "single"), (["SI"], "single"), (["SS"], "single"), (["HScan"], "single")]
for combo in (["ESD", "SS"], ["SS", "HScan"], ["ESD", "HScan"], ["ESD", "SS", "HScan"]):
    SETS += [(combo, "per modality"), (combo, "overall")]
TOPN = (1, 2, 3, 5)
CUTS4 = [-5, -3, -1]
SCHEMES = {"4cat": ([-5, -3, -1], "ordinal4", "4-category (≥ −1 | −3 to −1 | −5 to −3 | < −5)"),
           "3cat": ([-3, -1], "ordinal3", "3-category (≥ −1 | −3 to −1 | < −3)")}
SEED, NBOOT, REPS = 2026, 2000, 10
COL = {"ESD": "#2a78d6", "SI": "#eb6834", "SS": "#1baf7a", "HScan": "#6250d6"}
INK, INK2 = "#0b0b0b", "#52514e"
S = "LMM eGFR Slope (mL/min/1.73m2/year)"


def set_name(mods, mode):
    base = " + ".join(LBL[m] for m in mods)
    return base if mode == "single" else f"{base} ({mode})"


def load():
    coh = pd.read_excel(COHORT, "Cohort")
    coh = coh[coh["Kidney ROI (first visit, allowed size)"]].set_index("Patient Case ID")
    imgs = first_visit_images(COHORT)
    wanted = tuple(f"{m}|{b}|" for m, b in BIN_WIDTHS.items())
    recs, missing = [], []
    for pid, stem in zip(imgs["Patient Case ID"], imgs.stem):
        path = SHARED / stem if (SHARED / stem).exists() else EXTRA / stem
        if not path.exists():
            missing.append(stem); continue
        with open(path, "rb") as fh:
            r = pickle.load(fh)
        recs.append({"pid": pid, "src": "shared" if path.parent == SHARED else "extra",
                     **{k: v for k, v in r.items() if k.startswith(wanted)}})
    if missing:
        raise FileNotFoundError(f"{len(missing)} images without features, e.g. {missing[:5]}; run extract_multimap_missing.py")
    I = pd.DataFrame(recs)
    X = I.drop(columns="src").groupby("pid").mean()
    keep = []
    for c in X.columns:
        mod, bw, name = c.split("|", 2)
        flt, cls, _ = name.split("_", 2)
        if BIN_WIDTHS.get(mod) != bw or flt not in FILTERS or cls not in CLASSES:
            continue
        v = X[c]
        if v.notna().sum() < 0.9 * len(v) or v.std() == 0:
            continue
        keep.append(c)
    X = X[keep].reindex(coh.index)
    X.columns = [f"{c.split('|')[0]}__{c.split('|', 2)[2]}" for c in X.columns]
    print(f"{len(I)} images ({(I.src == 'extra').sum()} from extra cache), {len(X)} patients, {X.shape[1]} features")
    return X, coh[S].to_numpy(float)


def c_index(y, s):
    num = den = 0.0
    for lo in np.unique(y):
        for hi in np.unique(y):
            if hi <= lo:
                continue
            d = s[y == hi][:, None] - s[y == lo][None, :]
            num += (d > 0).sum() + 0.5 * (d == 0).sum(); den += d.size
    return num / den


def fit_score(xtr, ytr, xte):
    from statsmodels.miscmodels.ordinal_model import OrderedModel
    r = OrderedModel(ytr, xtr, distr="logit").fit(method="bfgs", disp=False, maxiter=2000)
    return xte @ np.asarray(r.params[: xtr.shape[1]])


def cv_scores(X, y4, names, topn, overall, splits):
    sys.path.insert(0, str(MODEL_CODE))
    warnings.filterwarnings("ignore")
    from nested_loocv_radiomics import Preprocessor, select_top_spearman_per_modality
    score = np.full(len(y4), np.nan); chosen = []
    for tr, te in splits:
        prep, xtr = Preprocessor.fit(X[tr], names, 0.90, y4[tr], "spearman")
        xte = prep.transform(X[te])
        idx, _ = select_top_spearman_per_modality(xtr, y4[tr], prep.selected_features, topn, overall)
        score[te] = fit_score(xtr[:, idx], y4[tr], xte[:, idx])
        chosen += [prep.selected_features[i] for i in idx]
    return score, chosen


def boot_ci(y, s, n=NBOOT, seed=SEED):
    rng = np.random.default_rng(seed)
    i1, i0 = np.where(y == 1)[0], np.where(y == 0)[0]
    a = [roc_auc_score(y[np.r_[a1, a0]], s[np.r_[a1, a0]])
         for a1, a0 in ((rng.choice(i1, len(i1)), rng.choice(i0, len(i0))) for _ in range(n))]
    return np.percentile(a, [2.5, 97.5])


def run(X, slope, mods, mode, topn, cuts=CUTS4):
    names = [c for c in X.columns if c.split("__", 1)[0] in mods]
    overall = mode != "per modality"
    Xm = X[names].to_numpy(float)
    y4 = sum((slope < c).astype(int) for c in cuts)
    y3 = (slope < -3).astype(int)
    s, chosen = cv_scores(Xm, y4, names, topn, overall, LeaveOneOut().split(Xm))
    auc = roc_auc_score(y3, s); lo, hi = boot_ci(y3, s)
    kf = []
    for r in range(REPS):
        for tr, te in StratifiedKFold(5, shuffle=True, random_state=SEED + r).split(Xm, y4):
            sk, _ = cv_scores(Xm, y4, names, topn, overall, [(tr, te)])
            kf.append(roc_auc_score(y3[te], sk[te]))
    freq = pd.Series(chosen).value_counts() / len(y4)
    n_model = topn * (len(mods) if mode == "per modality" else 1)
    return dict(feature_set=set_name(mods, mode), modalities=" + ".join(mods), selection=mode, top_n=topn,
                features_in_model=n_model, n_candidates=len(names), LOOCV_AUC_lt3=auc, CI_low=lo, CI_high=hi,
                LOOCV_Cindex_4cat=c_index(y4, s), kfold5x10_mean_AUC_lt3=np.mean(kf),
                kfold5x10_AUC_IQR=f"{np.percentile(kf, 25):.2f}–{np.percentile(kf, 75):.2f}",
                selected=", ".join(f"{f} ({v:.2f})" for f, v in freq.head(max(5, n_model)).items())), s


def heatmap(res, path, scheme_label="4-category"):
    rows = [set_name(m, md) for m, md in SETS]
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 9.5))
    for ax, col, title in zip(axes, ["LOOCV_AUC_lt3", "kfold5x10_mean_AUC_lt3"],
                              ["LOOCV AUC (slope < −3)", "5-fold × 10 mean AUC (slope < −3)"]):
        M = res.pivot_table(index="feature_set", columns="top_n", values=col).reindex(rows)[list(TOPN)].to_numpy()
        im = ax.imshow(M, cmap="RdBu_r", vmin=0.25, vmax=0.75, aspect="auto")
        for i in range(M.shape[0]):
            for j in range(M.shape[1]):
                ax.text(j, i, f"{M[i, j]:.2f}", ha="center", va="center", fontsize=9,
                        color="white" if abs(M[i, j] - 0.5) > 0.17 else INK,
                        fontweight="bold" if M[i, j] >= 0.65 else "normal")
        ax.set_xticks(range(len(TOPN)), [f"top {n}" for n in TOPN]); ax.xaxis.tick_top()
        ax.set_yticks(range(len(rows)), rows if ax is axes[0] else [""] * len(rows), fontsize=9)
        for b in (3.5, 5.5, 7.5, 9.5):
            ax.axhline(b, color="white", lw=2)
        ax.set_title(title, fontsize=11, color=INK, pad=24)
        ax.tick_params(length=0, colors=INK2)
    cb = fig.colorbar(im, ax=axes, fraction=0.03, pad=0.02)
    cb.set_label("held-out AUC for slope < −3 (0.5 = chance); bold ≥ 0.65", color=INK2)
    fig.suptitle(f"Ordinal {scheme_label} (proportional odds) top-N radiomics, scored for eGFR slope < −3\n"
                 "eGFR follow-up ≥ 1.5 yr, kidney-ROI patients (n = 130; < −3: 25)", fontsize=12, color=INK)
    fig.savefig(path.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def roc_singles(res, scores, y3, path, scheme_label="4-category"):
    fig, axes = plt.subplots(1, len(TOPN), figsize=(5 * len(TOPN), 7))
    for ax, n in zip(axes, TOPN):
        for m in ["ESD", "SI", "SS", "HScan"]:
            fpr, tpr, _ = roc_curve(y3, scores[f"{LBL[m]} top{n}"].to_numpy())
            r = res[(res.feature_set == LBL[m]) & (res.top_n == n)].iloc[0]
            ax.step(fpr, tpr, where="post", color=COL[m], lw=2,
                    label=f"{LBL[m]}: {r.LOOCV_AUC_lt3:.2f} ({r.CI_low:.2f}–{r.CI_high:.2f}); 5-fold {r.kfold5x10_mean_AUC_lt3:.2f}")
        ax.plot([0, 1], [0, 1], color="#bdbcb6", lw=1, ls="--")
        ax.set_title(f"Top {n}", fontsize=11, color=INK)
        ax.set_xlabel("1 − specificity", color=INK); ax.set_ylabel("Sensitivity", color=INK)
        ax.set_xlim(0, 1); ax.set_ylim(0, 1.01); ax.set_aspect("equal")
        ax.grid(color="#e6e5e1", lw=0.6); ax.set_axisbelow(True)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        ax.tick_params(colors=INK2)
        ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.14), fontsize=8.5, frameon=False,
                  title="LOOCV AUC (95% CI); 5-fold×10 mean", title_fontsize=8.5)
    fig.suptitle(f"Single-modality ordinal {scheme_label} models, scored for eGFR slope < −3 — eGFR follow-up ≥ 1.5 yr "
                 f"(n = {len(y3)}; < −3: {y3.sum()})", fontsize=12, color=INK)
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    scheme = sys.argv[1] if len(sys.argv) > 1 else "4cat"   # 4cat (default) | 3cat
    cuts, prefix, slabel = SCHEMES[scheme]
    X, slope = load()
    y4 = sum((slope < c).astype(int) for c in cuts)
    y3 = (slope < -3).astype(int)
    print(f"{scheme} counts (low -> high decline):", np.bincount(y4).tolist(), "| < -3:", y3.sum())
    jobs = [(m, md, n) for m, md in SETS for n in TOPN]
    out = Parallel(n_jobs=18, verbose=0)(delayed(run)(X, slope, m, md, n, cuts) for m, md, n in jobs)
    res = pd.DataFrame([o[0] for o in out]).rename(columns={"LOOCV_Cindex_4cat": f"LOOCV_Cindex_{scheme}"})
    scores = pd.DataFrame({f"{set_name(m, md)} top{n}": o[1] for (m, md, n), o in zip(jobs, out)}, index=X.index)
    heatmap(res, RES / f"{prefix}_topN_lt3_min1p5yr_heatmap", slabel)
    roc_singles(res, scores, y3, RES / f"{prefix}_topN_lt3_min1p5yr_roc", slabel)
    with pd.ExcelWriter(RES / f"{prefix}_topN_lt3_min1p5yr.xlsx") as w:
        pd.DataFrame({"Item": ["Method"], "Value": [__doc__]}).to_excel(w, sheet_name="Read_Me", index=False)
        res.to_excel(w, sheet_name="Performance", index=False)
        scores.assign(slope=slope, group=y4, lt_minus3=y3).reset_index().to_excel(w, sheet_name="LOOCV scores", index=False)
    pd.set_option("display.width", 250)
    print(res[["feature_set", "top_n", "features_in_model", "LOOCV_AUC_lt3", "CI_low", "CI_high", f"LOOCV_Cindex_{scheme}",
               "kfold5x10_mean_AUC_lt3", "kfold5x10_AUC_IQR"]].round(3).to_string(index=False))

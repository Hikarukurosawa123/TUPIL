"""LOOCV: single-modality top-3 radiomics (ESD, SI, SS, H-scan) -> binary logistic regression for LMM slope < -3,
in the >= 1 yr and >= 1.5 yr expanded-RF Case 8 cohorts (kidney-ROI patients with radiomics).

Features   kidney_pr ROI (kidney minus pelvis), Original image type, visit mean over the FIRST-visit images that have an
           allowed RF size and a kidney ROI. Source tables (same extractor, extract_radiomics_3roi.py):
           noMinFU allRF table first, then the all-patients table for images not in it.
Pipeline   refitted inside every LOOCV training set (patient held out) with the TUPIL R-score code
           (Results_Rough/02_.../nested_loocv_radiomics.py): median imputation -> greedy |Pearson r| >= 0.9 filter
           (variance priority) -> standardisation -> top-3 |Spearman rho| with the binary outcome (within the one
           modality) -> LogisticRegression(penalty=None, class_weight='balanced').
Outcome    slope < -3 (1) vs >= -3 (0), slope = that cohort's own LMM fit.
Metrics    LOOCV AUC (held-out probabilities), 95% CI by 2000 stratified bootstrap resamples of the held-out
           predictions, sensitivity / specificity at p = 0.5 (balanced weights) and at the Youden cut-off of the
           held-out scores (optimistic), and how often each feature was chosen across the LOOCV folds.
Run with C:/Users/tupil/venvs/radiomics_env/python.exe
"""
import re
import sys
import warnings
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, roc_curve

ROOT = Path(__file__).resolve().parents[2]
RES = ROOT / "Results_Rough" / "case8_min1p5yr_expanded_rf"
sys.path.insert(0, str(ROOT / "Results_Rough" / "02_ci_Rscore_radiomics_top3overall_SI_HScan"))
import nested_loocv_radiomics as M  # noqa: E402

R = Path(r"\\Tupil\tupil data\Hikaru\POCUS_Kidney_Analysis\Analysis_scripts\Redo_analysis_1106_2026")
TABLES = [(R / "v10_figures/patient_groups/noMinFU_radiomics/radiomics_SI_HScan_allRF_noMinFU_jobs.csv",
           R / "v10_figures/patient_groups/noMinFU_radiomics/radiomics_SI_HScan_allRF_noMinFU.csv"),
          (R / "rscore_allpatients_Data1103/radiomics_SI_HScan_allpatients_jobs.csv",
           R / "rscore_allpatients_Data1103/radiomics_SI_HScan_allpatients.csv")]
AUDIT = R / "rscore_allpatients_Data1103/images_audit.csv"
COHORTS = {"≥ 1 yr": RES / "case8_min1yr_expanded_rf_slope_groups.xlsx",
           "≥ 1.5 yr": RES / "case8_min1p5yr_expanded_rf_slope_groups.xlsx"}
ALLOWED = {"2928 x 192", "2944 x 192", "2784 x 176", "2960 x 192", "2960 x 184", "2944 x 176",
           "3488 x 184", "2736 x 192", "2928 x 168", "2912 x 160"}
MODS = ["ESD", "SI", "SS", "HScan"]
TOP, SEED, NBOOT = 3, 2026, 2000
COL = {"ESD": "#2a78d6", "SI": "#eb6834", "SS": "#1baf7a", "HScan": "#6250d6"}  # categorical slots 1, 2, 3, 7
INK, INK2 = "#0b0b0b", "#52514e"

# ---------------- image-level features (first visit, allowed size, kidney ROI) ----------------
aud = pd.read_csv(AUDIT, parse_dates=["Imaging Date"])
aud = aud[aud["RF Size"].isin(ALLOWED) & aud["Kidney ROI"]]
aud["key"] = aud["Visit Folder"] + "/" + aud["RF File"]
feat_cols = None
imgs = []
for jf, rf in TABLES:
    j = pd.read_csv(jf, keep_default_na=False)[["png", "scan_id", "raw"]]
    hdr = pd.read_csv(rf, nrows=0).columns
    want = [c for c in hdr if any(c.startswith(f"kidney_pr__{m}_") for m in MODS)]
    feat_cols = feat_cols or want
    x = pd.read_csv(rf, usecols=["png", "status"] + feat_cols, low_memory=False)
    x = x[x.status == "ok"].merge(j, on="png")
    parts = x.raw.str.replace("/", "\\").str.split("\\")
    x["key"] = parts.str[-3] + "/" + parts.str[-1]
    imgs.append(x.merge(aud[["key", "Patient Case ID", "Imaging Date", "RF Size"]], on="key"))
img = pd.concat(imgs, ignore_index=True).drop_duplicates("key")  # noMinFU wins
print(f"first-visit-eligible images with radiomics: {len(img)} ({img['Patient Case ID'].nunique()} patients)")
rename = {c: c.replace("kidney_pr__", "").replace("_", "__", 1) for c in feat_cols}  # -> 'ESD__firstorder_Mean'


def visit_features(coh):
    m = img.merge(coh[["Patient Case ID", "Imaging Date"]].assign(**{"Imaging Date": pd.to_datetime(coh["Imaging Date"])}),
                  on=["Patient Case ID", "Imaging Date"])
    v = m.groupby("Patient Case ID")[feat_cols].mean().rename(columns=rename)
    v["n_images"] = m.groupby("Patient Case ID").size()
    return v


def bootstrap_ci(y, p, n=NBOOT, seed=SEED):
    rng = np.random.default_rng(seed)
    i1, i0 = np.where(y == 1)[0], np.where(y == 0)[0]
    a = [roc_auc_score(y[np.r_[s1, s0]], p[np.r_[s1, s0]])
         for s1, s0 in ((rng.choice(i1, len(i1)), rng.choice(i0, len(i0))) for _ in range(n))]
    return np.percentile(a, [2.5, 97.5])


def loocv(X, y, names):
    prob = np.empty(len(y))
    chosen = []
    for i in range(len(y)):
        tr = np.r_[0:i, i + 1:len(y)]
        prep, xt = M.Preprocessor.fit(X[tr], names, 0.9, y[tr], "variance")
        xs = prep.transform(X[[i]])
        top, _ = M.select_top_spearman_per_modality(xt, y[tr], prep.selected_features, TOP, True)
        chosen.append([prep.selected_features[k] for k in top])
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            lr = LogisticRegression(penalty=None, class_weight="balanced", max_iter=10000).fit(xt[:, top], y[tr])
        prob[i] = lr.predict_proba(xs[:, top])[0, 1]
    return prob, chosen


rows, featrows, preds, roc = [], [], [], {}
for cname, wb in COHORTS.items():
    coh = pd.read_excel(wb, "Cohort")
    coh = coh[coh["Kidney ROI (first visit, allowed size)"]]
    v = visit_features(coh)
    d = coh.set_index("Patient Case ID").join(v, how="inner")
    y = (d["Slope Group"] == "< -3").astype(int).to_numpy()
    print(f"\n[{cname}] kidney-ROI patients {len(coh)}, with radiomics {len(d)} (< -3: {y.sum()}, >= -3: {(y == 0).sum()}); "
          f"missing radiomics: {sorted(set(coh['Patient Case ID']) - set(d.index))}")
    for mod in MODS:
        names = [c for c in d.columns if c.startswith(f"{mod}__")]
        X = d[names].apply(pd.to_numeric, errors="coerce").to_numpy(float)
        p, chosen = loocv(X, y, names)
        auc = roc_auc_score(y, p)
        lo, hi = bootstrap_ci(y, p)
        fpr, tpr, thr = roc_curve(y, p)
        j = np.argmax(tpr - fpr)
        pred = (p >= 0.5).astype(int)
        rows.append(dict(cohort=cname, modality=mod, n=len(y), n_lt_minus3=int(y.sum()), LOOCV_AUC=auc,
                         CI_low=lo, CI_high=hi,
                         sens_at_0p5=pred[y == 1].mean(), spec_at_0p5=1 - pred[y == 0].mean(),
                         youden_cut=thr[j], sens_youden=tpr[j], spec_youden=1 - fpr[j]))
        cnt = pd.Series([f for c in chosen for f in c]).value_counts()
        for f, k in cnt.items():
            featrows.append(dict(cohort=cname, modality=mod, feature=f, folds_selected=int(k), fraction=k / len(y)))
        preds.append(pd.DataFrame({"cohort": cname, "modality": mod, "Patient Case ID": d.index,
                                   "slope": d["LMM eGFR Slope (mL/min/1.73m2/year)"].to_numpy(), "y_lt_minus3": y,
                                   "LOOCV_prob": p}))
        roc[(cname, mod)] = (fpr, tpr, auc, lo, hi)
        print(f"  {mod:5s} AUC {auc:.3f} (95% CI {lo:.3f}-{hi:.3f}) | top features: "
              + "; ".join(f"{f.split('__', 1)[1]} {k}/{len(y)}" for f, k in cnt.head(3).items()))

res = pd.DataFrame(rows)
fig, axes = plt.subplots(1, 2, figsize=(10.5, 6))
for ax, cname in zip(axes, COHORTS):
    for mod in MODS:
        fpr, tpr, auc, lo, hi = roc[(cname, mod)]
        ax.step(fpr, tpr, where="post", color=COL[mod], lw=2, label=f"{'H-scan' if mod == 'HScan' else mod}: AUC {auc:.2f} ({lo:.2f}–{hi:.2f})")
    ax.plot([0, 1], [0, 1], color="#bdbcb6", lw=1, ls="--")
    r0 = res[res.cohort == cname].iloc[0]
    ax.set_title(f"eGFR follow-up {cname} (n = {r0.n}; slope < −3: {r0.n_lt_minus3})", fontsize=11, color=INK)
    ax.set_xlabel("1 − specificity", color=INK); ax.set_ylabel("Sensitivity", color=INK)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1.01); ax.set_aspect("equal")
    ax.grid(color="#e6e5e1", lw=0.6); ax.set_axisbelow(True)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.tick_params(colors=INK2)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.14), ncol=2, fontsize=9, frameon=False)
fig.suptitle("LOOCV: single-modality top-3 radiomics → logistic regression, predicting eGFR slope < −3", fontsize=12, color=INK)
fig.tight_layout()
fig.savefig(RES / "loocv_top3_per_modality_roc.png", dpi=300)
fig.savefig(RES / "loocv_top3_per_modality_roc.pdf")

with pd.ExcelWriter(RES / "loocv_top3_per_modality.xlsx") as w:
    pd.DataFrame({"Item": ["Method"], "Value": [__doc__]}).to_excel(w, sheet_name="Read_Me", index=False)
    res.to_excel(w, sheet_name="Performance", index=False)
    pd.DataFrame(featrows).to_excel(w, sheet_name="Feature selection", index=False)
    pd.concat(preds).to_excel(w, sheet_name="LOOCV predictions", index=False)
pd.set_option("display.width", 250)
print("\n", res.round(3).to_string(index=False))

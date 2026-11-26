"""R-fibrosis score (SI + H-scan top-3, cumulative L1 ordinal; Results_Rough/02_...) by eGFR-slope group (< -3 vs >= -3)
for the >= 1 yr and >= 1.5 yr expanded-RF Case 8 cohorts, kidney-ROI patients.

Score at each patient's FIRST imaging visit (the slope's time zero):
  * one of the 174 R-score model patients whose model visit is that same date -> held-out CV score (mean of 10 repeats)
  * otherwise -> final model refitted on all 174 (same pipeline/seed as km_rfibrosis_ci_allpatients.py) applied to the
    visit mean of the first-visit kidney-ROI images of an allowed RF size (kidney_pr ROI, Original features; noMinFU
    all-RF table, then the all-patients table; same extractor as the model inputs).
R-score = expected class sum_k k * P(class k), classes ci 0 / 1 / 2 / HD (0-3). Also P(ci > 0).
Test: Mann-Whitney U and AUC of the R-score for slope < -3.
Run with C:/Users/tupil/venvs/radiomics_env/python.exe
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu
from sklearn.metrics import roc_auc_score

ROOT = Path(__file__).resolve().parents[2]
RES = ROOT / "Results_Rough" / "case8_min1p5yr_expanded_rf"
RS = ROOT / "Results_Rough" / "02_ci_Rscore_radiomics_top3overall_SI_HScan"
sys.path.insert(0, str(RS))
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
SEED = 2026
GROUPS = [">= -3", "< -3"]
LABEL = {">= -3": "slope ≥ −3", "< -3": "slope < −3"}
COL = {">= -3": "#2a78d6", "< -3": "#eb6834"}
INK, INK2 = "#0b0b0b", "#52514e"

# ---------------- final R-fibrosis model (identical to km_rfibrosis_ci_allpatients.py) ----------------
meta = pd.read_csv(RS / "parent_context" / "kidney_pr_all_metadata.csv")
feat = pd.read_csv(RS / "parent_context" / "kidney_pr_all_visit_mean_features.csv")
tr = feat.drop(columns=["ci_grade"], errors="ignore").merge(meta[["sample_id", "ci_grade"]], on="sample_id")
cols = [c for c in tr.columns if c.startswith(("SI__", "HScan__"))]
cols = [c for c in cols if pd.to_numeric(tr[c], errors="coerce").notna().any()]
y, _ = M.encode_target(tr.ci_grade, ["0", "1", "2", "3"], None, None)
X = tr[cols].apply(pd.to_numeric, errors="coerce").to_numpy(float)
best_c, _ = M.tune_c(X, y, cols, "ordinal", list(M.DEFAULT_C_GRID), 0.9, 5, SEED, "balanced", 4, None, False, 0.05, 3, True, "variance")
prep, Xt = M.Preprocessor.fit(X, cols, 0.9, y, "variance")
top, _ = M.select_top_spearman_per_modality(Xt, y, prep.selected_features, 3, True)
top_feats = [prep.selected_features[i] for i in top]
model = M.CumulativeL1Ordinal(best_c, SEED, "balanced").fit(Xt[:, top], y, 4)
print(f"final R-fibrosis model: C = {best_c}, features {top_feats}")

# ---------------- CV scores of the 174 model patients ----------------
cv = pd.read_csv(RS / "outer_predictions_per_sample_mean.csv")
coh174 = pd.read_excel(RS / "parent_context" / "cohort_source_firstVisit_174.xlsx")
coh174["sample_id"] = coh174["REDCap Scan ID"].str.strip()
coh174["pid"] = coh174["Participant ID (POCUS)"].str.strip()
cv = cv.merge(coh174[["sample_id", "pid", "Visit date"]], on="sample_id")
cv["Visit date"] = pd.to_datetime(cv["Visit date"])
cv = cv.set_index("pid")

# ---------------- first-visit image features ----------------
aud = pd.read_csv(AUDIT, parse_dates=["Imaging Date"])
aud = aud[aud["RF Size"].isin(ALLOWED) & aud["Kidney ROI"]]
aud["key"] = aud["Visit Folder"] + "/" + aud["RF File"]
need = ["kidney_pr__" + c.replace("__", "_", 1) for c in cols]
imgs = []
for jf, rf in TABLES:
    j = pd.read_csv(jf, keep_default_na=False)[["png", "raw"]]
    x = pd.read_csv(rf, usecols=["png", "status"] + need, low_memory=False)
    x = x[x.status == "ok"].merge(j, on="png")
    parts = x.raw.str.replace("/", "\\").str.split("\\")
    x["key"] = parts.str[-3] + "/" + parts.str[-1]
    imgs.append(x.merge(aud[["key", "Patient Case ID", "Imaging Date"]], on="key"))
img = pd.concat(imgs, ignore_index=True).drop_duplicates("key")


def rscores(coh):
    m = img.merge(coh[["Patient Case ID", "Imaging Date"]].assign(**{"Imaging Date": pd.to_datetime(coh["Imaging Date"])}),
                  on=["Patient Case ID", "Imaging Date"])
    v = m.groupby("Patient Case ID")[need].mean()
    v.columns = cols
    P = model.predict_proba(prep.transform(v.reindex(columns=cols).to_numpy(float))[:, top])
    out = pd.DataFrame({"R_score": P @ np.arange(4), "P_ci_gt0": 1 - P[:, 0], "score_source": "final model"}, index=v.index)
    for pid, r in coh.set_index("Patient Case ID").iterrows():
        if pid in cv.index and cv.at[pid, "Visit date"] == pd.to_datetime(r["Imaging Date"]):
            out.loc[pid, ["R_score", "P_ci_gt0", "score_source"]] = [cv.at[pid, "r_score"], 1 - cv.at[pid, "probability_0"],
                                                                    "CV held-out (model patient, same visit)"]
    return out


fig, axes = plt.subplots(1, 2, figsize=(9, 4.8), sharey=True)
rows, per = [], []
rng = np.random.default_rng(7)
for ax, (name, wb) in zip(axes, COHORTS.items()):
    coh = pd.read_excel(wb, "Cohort")
    coh = coh[coh["Kidney ROI (first visit, allowed size)"]]
    d = coh.set_index("Patient Case ID").join(rscores(coh), how="left")
    per.append(d.assign(cohort=name))
    vals = [d.loc[d["Slope Group"] == g, "R_score"].dropna().to_numpy(float) for g in GROUPS]
    p = mannwhitneyu(*vals).pvalue
    ok = d.R_score.notna()
    auc = roc_auc_score((d.loc[ok, "Slope Group"] == "< -3").astype(int), d.loc[ok, "R_score"].astype(float))
    for i, (g, v) in enumerate(zip(GROUPS, vals)):
        bp = ax.boxplot([v], positions=[i], widths=0.5, showfliers=False, patch_artist=True,
                        medianprops=dict(color=INK, lw=1.6), whiskerprops=dict(color=INK2), capprops=dict(color=INK2),
                        boxprops=dict(edgecolor=INK2))
        bp["boxes"][0].set_facecolor(COL[g]); bp["boxes"][0].set_alpha(0.28)
        ax.scatter(i + rng.uniform(-0.15, 0.15, len(v)), v, s=14, color=COL[g], edgecolor="white", lw=0.5, zorder=3)
    ax.set_xticks(range(2), [f"{LABEL[g]}\n(n = {len(v)})" for g, v in zip(GROUPS, vals)], color=INK)
    ax.set_xlim(-0.6, 1.6); ax.set_ylim(-0.1, 3.1)
    ax.set_title(f"eGFR follow-up {name}\nMann-Whitney p = {p:.3f}; AUC (< −3) = {auc:.2f}", fontsize=10, color=INK)
    ax.grid(axis="y", color="#e6e5e1", lw=0.6); ax.set_axisbelow(True)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.tick_params(colors=INK2)
    rows.append(dict(cohort=name, n=int(ok.sum()), not_scored=int((~ok).sum()),
                     **{f"median R {LABEL[g]}": np.median(v) for g, v in zip(GROUPS, vals)},
                     **{f"IQR R {LABEL[g]}": f"{np.percentile(v, 25):.2f}–{np.percentile(v, 75):.2f}" for g, v in zip(GROUPS, vals)},
                     mann_whitney_p=p, AUC_lt3=auc,
                     CV_scored=int((d.score_source == "CV held-out (model patient, same visit)").sum()),
                     final_model_scored=int((d.score_source == "final model").sum())))
axes[0].set_ylabel("R-fibrosis score (expected ci class, 0–3)", color=INK)
fig.suptitle("R-fibrosis score (SI + H-scan top-3) by eGFR slope group — kidney-ROI patients", fontsize=12, color=INK)
fig.tight_layout()
fig.savefig(RES / "rfibrosis_by_slope_lt3.png", dpi=300, bbox_inches="tight")
fig.savefig(RES / "rfibrosis_by_slope_lt3.pdf", bbox_inches="tight")
st = pd.DataFrame(rows)
with pd.ExcelWriter(RES / "rfibrosis_by_slope_lt3.xlsx") as w:
    pd.DataFrame({"Item": ["Method"], "Value": [__doc__]}).to_excel(w, sheet_name="Read_Me", index=False)
    st.to_excel(w, sheet_name="Stats", index=False)
    pd.concat(per).reset_index().rename(columns={"index": "Patient Case ID"})[
        ["cohort", "Patient Case ID", "Native/Transplant", "Imaging Date", "LMM eGFR Slope (mL/min/1.73m2/year)",
         "Slope Group", "R_score", "P_ci_gt0", "score_source"]].to_excel(w, sheet_name="Per patient", index=False)
pd.set_option("display.width", 250)
print(st.round(3).to_string(index=False))

"""Kaplan-Meier for ALL scanned patients: R-fibrosis predicted ci 0 vs ci > 0.

R-fibrosis model = SI + H-scan, top-3 overall |Spearman|, cumulative L1 ordinal (classes ci 0/1/2/HD),
from Results_Rough/02_ci_Rscore_radiomics_top3overall_SI_HScan (174 first-visit patients).
  * 174 model patients: their held-out (nested CV, mean of 10 repeats) probabilities, at the model visit.
  * every other patient: FINAL model refitted on all 174 with the identical pipeline (inner 5-fold C tuning,
    seed 2026; median impute -> |r|>=0.9 variance-priority filter -> scale -> top-3 |Spearman| -> L1 cumulative
    logistic) and applied to the visit-mean radiomics of their first usable Data_2026_1103 visit
    (rscore_allpatients_Data1103/radiomics_SI_HScan_allpatients.csv; same extractor and kidney_pr ROI).
Group: predicted ci > 0 if P(ci > 0) = 1 - P(ci 0) >= 0.5 (the model's own ci0 | ci>=1 boundary).

Endpoint (time zero = scored scan date)
  baseline  = first eGFR on/after the scan (same day allowed); transplant eGFR restricted to the graft in place
  40% decline (confirmed) = eGFR <= 60% of baseline, with another <= 60% value >= 30 d later and no value
              > 60% in between; event date = first low value
  RRT       = dialysis start (ProspBx_eGFR_RRT_Est), graft failure, or any new transplant after baseline
  censoring = last eGFR date. Patients on RRT at baseline are excluded.
Run with C:/Users/tupil/venvs/radiomics_env/python.exe (scikit-learn 1.6.1, as the model run).
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import chi2

ROOT = Path(__file__).resolve().parents[2]
D = ROOT / "All_csv"
RS = ROOT / "Results_Rough" / "02_ci_Rscore_radiomics_top3overall_SI_HScan"
ALLP = Path(r"\\Tupil\tupil data\Hikaru\POCUS_Kidney_Analysis\Analysis_scripts\Redo_analysis_1106_2026"
            r"\rscore_allpatients_Data1103")
OUT = ROOT / "Results_Rough" / "km_rfibrosis_ci_allpatients"
OUT.mkdir(parents=True, exist_ok=True)
CONFIRM_DAYS = 30
SEED = 2026

sys.path.insert(0, str(RS))
import nested_loocv_radiomics as M  # noqa: E402  (model code used for the R-fibrosis run)

# =============== 1. final model on the 174 ===============
meta = pd.read_csv(RS / "parent_context" / "kidney_pr_all_metadata.csv")
feat = pd.read_csv(RS / "parent_context" / "kidney_pr_all_visit_mean_features.csv")
tr = feat.drop(columns=["ci_grade"], errors="ignore").merge(meta[["sample_id", "ci_grade"]], on="sample_id")
cols = [c for c in tr.columns if c.startswith(("SI__", "HScan__"))]
cols = [c for c in cols if pd.to_numeric(tr[c], errors="coerce").notna().any()]
y, labels = M.encode_target(tr.ci_grade, ["0", "1", "2", "3"], None, None)
X = tr[cols].apply(pd.to_numeric, errors="coerce").to_numpy(float)
settings = dict(task_type="ordinal", c_grid=list(M.DEFAULT_C_GRID), correlation_threshold=0.9,
                inner_folds=5, class_weight="balanced", n_classes=4)
best_c, _ = M.tune_c(X, y, cols, settings["task_type"], settings["c_grid"], settings["correlation_threshold"],
                     settings["inner_folds"], SEED, settings["class_weight"], 4, None, False, 0.05, 3, True, "variance")
prep, Xt = M.Preprocessor.fit(X, cols, 0.9, y, "variance")
top, _ = M.select_top_spearman_per_modality(Xt, y, prep.selected_features, 3, True)
top_feats = [prep.selected_features[i] for i in top]
model = M.CumulativeL1Ordinal(best_c, SEED, "balanced").fit(Xt[:, top], y, 4)
coef = pd.DataFrame(model.coefficients(), columns=top_feats, index=["P(ci>0)", "P(ci>1)", "P(ci>2 / HD)"])
coef["intercept"] = [m.intercept_[0] for m in model.models_]
print(f"Final model: C = {best_c}; features = {top_feats}\n{coef.round(4)}")


def score(frame):
    x = frame.reindex(columns=cols).apply(pd.to_numeric, errors="coerce").to_numpy(float)
    return model.predict_proba(prep.transform(x)[:, top])


# =============== 2. scores: model patients (CV) + all other patients (final model) ===============
cv = pd.read_csv(RS / "outer_predictions_per_sample_mean.csv")
coh = pd.read_excel(RS / "parent_context" / "cohort_source_firstVisit_174.xlsx")
coh["sample_id"] = coh["REDCap Scan ID"].str.strip()
coh["pid"] = coh["Participant ID (POCUS)"].str.strip()
cv = cv.merge(coh[["sample_id", "pid", "Visit date", "Kidney type", "Case type"]], on="sample_id")
s_cv = pd.DataFrame({"pid": cv.pid, "scan_date": pd.to_datetime(cv["Visit date"]), "ktype_model": cv["Kidney type"],
                     "hd_at_scan": cv["Case type"].eq("hemodialysis"), "biopsy_ci": cv.true_index.where(cv.true_index < 3),
                     "P_ci0": cv.probability_0, "r_score": cv.r_score, "score_source": "CV held-out (model patient)"})

jobs = pd.read_csv(ALLP / "radiomics_SI_HScan_allpatients_jobs.csv", keep_default_na=False)[["png", "scan_id"]]
need = ["kidney_pr__" + c.replace("__", "_", 1) for c in cols]
rad = pd.read_csv(ALLP / "radiomics_SI_HScan_allpatients.csv", usecols=["png"] + need, low_memory=False)
v = jobs.merge(rad, on="png", how="left").groupby("scan_id")[need].mean()
v.columns = cols
aud = pd.read_csv(ALLP / "images_audit.csv", parse_dates=["Imaging Date"])
vdate = aud[aud["Qualifying (size + ROI)"]].groupby("Patient Case ID")["Imaging Date"].first()
ok = v[top_feats].notna().all(axis=1)
P = score(v[ok])
s_new = pd.DataFrame({"pid": v.index[ok], "scan_date": vdate.reindex(v.index[ok]).to_numpy(), "ktype_model": None,
                      "hd_at_scan": False, "biopsy_ci": np.nan, "P_ci0": P[:, 0], "r_score": P @ np.arange(4),
                      "score_source": "final model (not a model patient)"})
s_new = s_new[~s_new.pid.isin(s_cv.pid)]
# sanity: final model on the model patients' own visits vs their CV scores
chk = pd.Series(score(tr[cols])[:, 0], index=tr.sample_id).reindex(cv.sample_id)
print(f"final-model P(ci0) vs CV P(ci0) on the 174: r = {np.corrcoef(chk, cv.probability_0)[0, 1]:.3f}; "
      f"group agreement = {((chk.to_numpy() < 0.5) == (cv.probability_0 < 0.5)).mean():.3f}")
sc = pd.concat([s_cv, s_new], ignore_index=True)
sc["scan_date"] = pd.to_datetime(sc.scan_date)
sc["R_group"] = np.where(1 - sc.P_ci0 >= 0.5, "R-fibrosis ci > 0", "R-fibrosis ci 0")
flow = [("radiomics visits in Data_2026_1103 (first usable visit)", len(v)),
        ("  with complete top-3 features", int(ok.sum())),
        ("model patients (CV score)", len(s_cv)), ("other patients scored by final model", len(s_new)),
        ("total scored", len(sc))]

# kidney type of the scored scan (POCUS export; nearest scan date)
pocus = pd.read_csv(D / "PROSPECTIVEBIOPSYFIB-ExportPOCUSScan_DATA_LABELS_2026-11-13_1007.csv", encoding="latin1")
pocus["d"] = pd.to_datetime(pocus["POCUS Scan Date"], dayfirst=True, errors="coerce")
pocus = pocus[pocus["Is this kidney a failed kidney?"] != "Yes"]


def kidney_type(r):
    if isinstance(r.ktype_model, str):
        return r.ktype_model
    p = pocus[pocus["Participant ID:"] == r.pid]
    if p.empty:
        return None
    return p.loc[(p.d - r.scan_date).abs().idxmin(), "Kidney Type"]


sc["ktype"] = sc.apply(kidney_type, axis=1)

# =============== 3. endpoint ===============
bw = pd.read_excel(D / "PROSPBX_DataRequest_10Oct2026.xlsx", "PTX BW")
bw = pd.DataFrame({"pid": bw["StudyIDRequestList::StudyID"], "date": pd.to_datetime(bw["Date"]),
                   "egfr": bw["estGFRckd_epi_2021"], "txdate": pd.to_datetime(bw["StudyIDRequestList::TXDate"]),
                   "src": "Transplant"})
nk = pd.read_csv(D / "PROSPECTIVEBIOPSYFIB-DataEntryNativeKidne_DATA_LABELS_2026-10-22_1519.csv", encoding="latin1")
ecol = [c for c in nk.columns if c.startswith("Predicted eGFR (CKD-EPI")][0]
nk = pd.DataFrame({"pid": nk["Participant ID:"],
                   "date": pd.to_datetime(nk["Date of Lab Results:"], dayfirst=True, errors="coerce"),
                   "egfr": pd.to_numeric(nk[ecol], errors="coerce"), "txdate": pd.NaT, "src": "Native"})
eg = pd.concat([bw, nk]).dropna(subset=["date", "egfr"])
ps = pd.read_excel(D / "ProspBx_eGFR_RRT_Est_26-08-26.xlsx", "Patient_Summary")
ps["rrt_date"] = pd.to_datetime(ps["RRT_initiation_date (dd/mm/yy)"], errors="coerce")
ps = ps.set_index("patient_id")
tx = pd.read_csv(D / "PROSPECTIVEBIOPSYFIB-ExportTransplantHist_DATA_LABELS_2026-11-13_1009.csv", encoding="latin1")
tx["txd"] = pd.to_datetime(tx["Date of Transplant"], dayfirst=True, errors="coerce")
tx["fail"] = pd.to_datetime(tx["Transplant Failed Date"], dayfirst=True, errors="coerce")
txdates = tx.dropna(subset=["txd"]).groupby("Participant ID:").txd.apply(list)
txfail = tx.dropna(subset=["fail"]).groupby("Participant ID:").fail.apply(list)


def endpoint(r):
    if r.hd_at_scan:
        return {"status": "excluded: HD at scan"}
    if r.ktype not in ("Transplant", "Native"):
        return {"status": "excluded: kidney type unknown"}
    e = eg[(eg.pid == r.pid) & (eg.src == r.ktype)]
    if r.ktype == "Transplant":
        e = e[e.txdate == e[e.txdate <= r.scan_date].txdate.max()]
    e = e.groupby("date", as_index=False).egfr.mean().sort_values("date")
    post = e[e.date >= r.scan_date]
    note = ps.Notes.get(r.pid) if r.pid in ps.index else None
    rd = ps.rrt_date.get(r.pid) if r.pid in ps.index else pd.NaT
    if note == "Baseline RRT":
        return {"status": "excluded: RRT at baseline"}
    if len(post) == 0:
        return {"status": "excluded: no eGFR after scan"}
    b = post.iloc[0]
    fu = post.iloc[1:].reset_index(drop=True)
    thr = 0.6 * b.egfr
    if pd.notna(rd) and rd <= b.date:
        return {"status": "excluded: RRT at baseline"}
    cands = ([rd] if pd.notna(rd) else []) + [t for t in txdates.get(r.pid, []) if t > b.date]
    if r.ktype == "Transplant":
        cands += [t for t in txfail.get(r.pid, []) if t > b.date]
    rrt = min(cands) if cands else pd.NaT
    dec = pd.NaT
    for i in fu.index[fu.egfr <= thr]:
        later = fu.loc[i + 1:]
        above = later[later.egfr > thr]
        stop = above.date.iloc[0] if len(above) else pd.Timestamp.max
        okk = later[(later.egfr <= thr) & (later.date >= fu.date[i] + pd.Timedelta(days=CONFIRM_DAYS)) & (later.date < stop)]
        if len(okk):
            dec = fu.date[i]
            break
    if len(fu) == 0 and pd.isna(rrt):
        return {"status": "excluded: no follow-up eGFR"}
    evd = [x for x in (dec, rrt) if pd.notna(x)]
    end = min(evd) if evd else post.date.iloc[-1]
    return {"status": "analysed", "baseline_date": b.date, "baseline_egfr": round(b.egfr, 1),
            "days_scan_to_baseline": (b.date - r.scan_date).days, "last_egfr_date": post.date.iloc[-1],
            "decline40_date": dec, "rrt_date": rrt, "event": int(bool(evd)), "time_days": (end - r.scan_date).days}


res = pd.concat([sc, sc.apply(lambda r: pd.Series(endpoint(r)), axis=1)], axis=1)
res.to_csv(OUT / "km_rfibrosis_ci_allpatients_patient_level.csv", index=False)
a = res[res.status == "analysed"].copy()
a["event"] = a.event.astype(int)
a["time_yr"] = a.time_days.astype(float) / 365.25


# =============== 4. KM + log-rank ===============
def km(t, e):
    t, e = np.asarray(t, float), np.asarray(e, int)
    s, var, T, S, V = 1.0, 0.0, [0.0], [1.0], [0.0]
    for u in np.unique(t[e == 1]):
        n = (t >= u).sum(); d = ((t == u) & (e == 1)).sum()
        s *= 1 - d / n
        var += d / (n * (n - d)) if n > d else 0
        T.append(u); S.append(s); V.append(var)
    if t.max() > T[-1]:
        T.append(t.max()); S.append(s); V.append(var)
    S, V = np.array(S), np.array(V)
    se = S * np.sqrt(V)
    return np.array(T), S, np.clip(S - 1.96 * se, 0, 1), np.clip(S + 1.96 * se, 0, 1)


def logrank(t, e, g):
    t, e, g = np.asarray(t, float), np.asarray(e, int), np.asarray(g)
    g1 = g == np.unique(g)[0]
    O = E = Vv = 0.0
    for u in np.unique(t[e == 1]):
        at = t >= u; n = at.sum(); n1 = (at & g1).sum()
        d = ((t == u) & (e == 1)).sum()
        O += ((t == u) & (e == 1) & g1).sum(); E += d * n1 / n
        if n > 1: Vv += d * (n1 / n) * (1 - n1 / n) * (n - d) / (n - 1)
    stat = (O - E) ** 2 / Vv
    return stat, chi2.sf(stat, 1)


COL = {0: "#2b6cb0", 1: "#c53030"}
LV = ["R-fibrosis ci 0", "R-fibrosis ci > 0"]
panels = [(a, f"All scored patients (n = {len(a)})"),
          (a[a.score_source.str.startswith("final")],
           f"Excluding the 174 model patients (n = {int(a.score_source.str.startswith('final').sum())})")]
tmax = np.ceil(a.time_yr.max() * 2) / 2
ticks = np.arange(0, tmax + 0.01, 0.5)
fig, axes = plt.subplots(2, 2, figsize=(13, 7.2), gridspec_kw={"height_ratios": [4, 1]})
summary = []
for j, (dd, title) in enumerate(panels):
    ax, axt = axes[0, j], axes[1, j]
    for k, lv in enumerate(LV):
        s = dd[dd.R_group == lv]
        T, S, lo, hi = km(s.time_yr, s.event)
        ax.step(T, S, where="post", color=COL[k], lw=2, label=f"{lv} (n={len(s)}, events={s.event.sum()})")
        ax.fill_between(T, lo, hi, step="post", color=COL[k], alpha=0.12)
        cz = s[s.event == 0].time_yr
        ax.plot(cz, [S[np.searchsorted(T, c, side="right") - 1] for c in cz], "|", color=COL[k], ms=7)
        axt.text(-0.06, 0.75 - 0.45 * k, lv.replace("R-fibrosis ", ""), color=COL[k], ha="right", va="center",
                 transform=axt.transAxes, fontsize=9)
        for x in ticks:
            axt.text(x, 0.75 - 0.45 * k, str((s.time_yr >= x).sum()), ha="center", va="center", fontsize=9, color=COL[k])
        for yr in (1, 2, 3):
            i = np.searchsorted(T, yr, side="right") - 1
            summary.append(dict(panel=title, group=lv, n=len(s), events=int(s.event.sum()), year=yr,
                                event_free=S[i], ci_lo=lo[i], ci_hi=hi[i], at_risk=int((s.time_yr >= yr).sum())))
    stat, p = logrank(dd.time_yr, dd.event, dd.R_group)
    summary.append(dict(panel=title, group="log-rank", n=len(dd), events=int(dd.event.sum()), chi2=stat, p=p))
    ax.text(0.03, 0.05, f"Log-rank p = {p:.3g}", transform=ax.transAxes, fontsize=10)
    ax.set_title(title, fontsize=11)
    ax.set_xlim(0, tmax); ax.set_ylim(0, 1.02); ax.set_xticks(ticks)
    ax.set_ylabel("Free of 40% eGFR decline or RRT"); ax.set_xlabel("Years from scan")
    ax.legend(loc="lower left", bbox_to_anchor=(0, 0.1), fontsize=9, frameon=False)
    ax.grid(alpha=0.25)
    axt.set_xlim(0, tmax); axt.set_ylim(0, 1); axt.axis("off")
    axt.set_title("Number at risk", fontsize=9, loc="left")
fig.suptitle("R-fibrosis score (SI + H-scan top-3): predicted ci 0 vs ci > 0 [P(ci>0) ≥ 0.5]\n"
             f"Composite: confirmed 40% eGFR decline (2 values ≥{CONFIRM_DAYS} d apart) or RRT", fontsize=12)
fig.tight_layout()
fig.savefig(OUT / "km_rfibrosis_ci0_vs_gt0_allpatients.png", dpi=300)
fig.savefig(OUT / "km_rfibrosis_ci0_vs_gt0_allpatients.pdf")
pd.DataFrame(summary).to_csv(OUT / "km_summary_allpatients.csv", index=False)
coef.to_csv(OUT / "final_model_coefficients.csv")

# =============== console ===============
for f, n in flow: print(f"{f}: {n}")
print(res.status.value_counts().to_string())
print("analysed by kidney:", a.ktype.value_counts().to_dict(), "| by source:", a.score_source.value_counts().to_dict())
for dd, title in panels:
    print(f"\n[{title}]")
    for lv in LV:
        s = dd[dd.R_group == lv]
        print(f"  {lv:18s} n={len(s):3d} events={s.event.sum():3d} (decline {s.decline40_date.notna().sum()}, "
              f"RRT {s.rrt_date.notna().sum()}) median FU {s.time_yr.median():.2f} y")
print(pd.DataFrame(summary).round(3).to_string())

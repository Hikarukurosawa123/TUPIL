"""Baseline eGFR, ACR and KFRE (2-yr) by LMM eGFR-slope group (< -3 vs >= -3), for the >= 1 yr and >= 1.5 yr
expanded-RF Case 8 cohorts (output of case8_min1p5yr_expanded_rf.py), all patients and kidney-ROI patients.

Baseline rule = kfre_tx.py (Redo_analysis_1106_2026/egfr_slope_case348_all_1p5_2yr), inputs from All_csv:
  value closest to the imaging date within +-90 d (ties -> earlier); transplant patients: if that value is dated before
  transplant + 3 calendar months it is replaced by the first value on/after that cutoff if within 90 d of it, else missing.
  eGFR  transplant: PROSPBX_DataRequest_10Oct2026.xlsx PTX BW estGFRckd_epi_2021
        native:     PROSPECTIVEBIOPSYFIB-DataEntryNativeKidne_DATA_LABELS_2026-10-22_1519.csv  (CKD-EPI 2021)
  ACR   transplant: PTX URINE MACRRandom (if > 0, else CalcACRLevel), mg/mmol
        native:     same native CSV, 'Urine ACR (mg/mmol):'
  Age   visit_age of the VisitReview visit closest to imaging; sex from DemographicBaseline (both 13Oct2026, All_csv).
  KFRE  Tangri 4-variable, 2-yr, North American (S0 = 0.9750), UACR mg/g = mg/mmol x 8.84.
Test: two-sided Mann-Whitney U per panel.
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu

ROOT = Path(__file__).resolve().parents[2]
D = ROOT / "All_csv"
RES = ROOT / "Results_Rough" / "case8_min1p5yr_expanded_rf"
COHORTS = {"≥ 1 yr": RES / "case8_min1yr_expanded_rf_slope_groups.xlsx",
           "≥ 1.5 yr": RES / "case8_min1p5yr_expanded_rf_slope_groups.xlsx"}
WINDOW, S0_2YR = 90, 0.9750
GROUPS = [">= -3", "< -3"]
LABEL = {">= -3": "slope ≥ −3", "< -3": "slope < −3"}
COL = {">= -3": "#2a78d6", "< -3": "#eb6834"}  # categorical slots 1, 2 (validated reference palette)
INK, INK2 = "#0b0b0b", "#52514e"

# ---------------- inputs ----------------
bw = pd.read_excel(D / "PROSPBX_DataRequest_10Oct2026.xlsx", "PTX BW")
ur = pd.read_excel(D / "PROSPBX_DataRequest_10Oct2026.xlsx", "PTX URINE")
nk = pd.read_csv(D / "PROSPECTIVEBIOPSYFIB-DataEntryNativeKidne_DATA_LABELS_2026-10-22_1519.csv", encoding="latin1")
ecol = [c for c in nk.columns if c.startswith("Predicted eGFR (CKD-EPI")][0]
acol = [c for c in nk.columns if c.startswith("Urine ACR")][0]
nk_date = pd.to_datetime(nk["Date of Lab Results:"], dayfirst=True, errors="coerce")
ur["acr"] = ur.MACRRandom.where(ur.MACRRandom > 0, ur.CalcACRLevel)
ur.loc[ur.acr <= 0, "acr"] = np.nan
LABS = {
    ("egfr", "transplant"): pd.DataFrame({"pid": bw["StudyIDRequestList::StudyID"], "d": pd.to_datetime(bw.Date), "v": bw.estGFRckd_epi_2021}),
    ("egfr", "native"): pd.DataFrame({"pid": nk["Participant ID:"], "d": nk_date, "v": pd.to_numeric(nk[ecol], errors="coerce")}),
    ("acr", "transplant"): pd.DataFrame({"pid": ur["StudyIDRequestList::StudyID"], "d": pd.to_datetime(ur.UrineDate), "v": ur.acr}),
    ("acr", "native"): pd.DataFrame({"pid": nk["Participant ID:"], "d": nk_date, "v": pd.to_numeric(nk[acol], errors="coerce")}),
}
LABS = {k: v.dropna() for k, v in LABS.items()}
demo = pd.read_csv(D / "PROSPBX-FIBRO-US_DemographicBaseline_13Oct2026.csv").drop_duplicates("participant_id").set_index("participant_id")
vr = pd.read_csv(D / "PROSPBX-FIBRO-US_VisitReview_13Oct2026.csv").dropna(subset=["visit_age"])
vr["visit_date"] = pd.to_datetime(vr.visit_date, dayfirst=True)
txh = pd.read_csv(D / "PROSPECTIVEBIOPSYFIB-ExportTransplantHist_DATA_LABELS_2026-11-13_1009.csv", encoding="latin1")
txh["txd"] = pd.to_datetime(txh["Date of Transplant"], dayfirst=True, errors="coerce")


def pick(lab, pid, kt, scan, cutoff):
    s = LABS[(lab, kt)]
    s = s[s.pid == pid].assign(gap=lambda x: (x.d - scan).dt.days)
    w = s[s.gap.abs() <= WINDOW].assign(a=lambda x: x.gap.abs()).sort_values(["a", "gap"])
    sel = w.iloc[0] if len(w) else None
    if pd.notna(cutoff) and ((sel is not None and sel.d < cutoff) or (sel is None and scan < cutoff)):
        c = s[(s.d >= cutoff) & ((s.d - cutoff).dt.days <= WINDOW)].sort_values("d")
        sel = c.iloc[0] if len(c) else None
    return (np.nan, np.nan) if sel is None else (float(sel.v), int(sel.gap))


def baseline(coh):
    rows = []
    for _, r in coh.iterrows():
        pid, kt, scan = r["Patient Case ID"], r["Native/Transplant"], pd.to_datetime(r["Imaging Date"])
        cutoff = pd.NaT
        if kt == "transplant":
            t = txh.loc[(txh["Participant ID:"] == pid) & (txh.txd <= scan), "txd"].max()
            cutoff = t + pd.DateOffset(months=3) if pd.notna(t) else pd.NaT
        e, eg_gap = pick("egfr", pid, kt, scan, cutoff)
        a, ac_gap = pick("acr", pid, kt, scan, cutoff)
        v = vr[vr.participant_id == pid]
        age = v.loc[(v.visit_date - scan).abs().idxmin(), "visit_age"] if len(v) else np.nan
        male = {"M": 1, "F": 0}.get(demo.gender.get(pid), np.nan)
        rows.append(dict(pid=pid, egfr=e, egfr_gap=eg_gap, acr=a, acr_gap=ac_gap, age=age, male=male))
    b = pd.DataFrame(rows).set_index("pid")
    L = (-0.2201 * (b.age / 10 - 7.036) + 0.2467 * (b.male - 0.5642) - 0.5567 * (b.egfr / 5 - 7.222)
         + 0.4510 * (np.log(b.acr * 8.84) - 5.137))
    b["kfre2"] = 100 * (1 - S0_2YR ** np.exp(L))
    return b


PANELS = [("egfr", "Baseline eGFR\n(mL/min/1.73 m²)", "linear"), ("acr", "Baseline ACR\n(mg/mmol, log scale)", "log"),
          ("kfre2", "KFRE 2-yr risk\n(%, log scale)", "log")]


def fmt_p(p):
    return "p < 0.001" if p < 0.001 else f"p = {p:.3f}" if p < 0.1 else f"p = {p:.2f}"


def plot(data, title, path):
    fig, axes = plt.subplots(1, 3, figsize=(11, 4.4))
    stats = []
    rng = np.random.default_rng(7)
    for ax, (var, ylab, scale) in zip(axes, PANELS):
        vals = [data.loc[data["Slope Group"] == g, var].dropna().to_numpy() for g in GROUPS]
        for i, (g, v) in enumerate(zip(GROUPS, vals)):
            if not len(v):
                continue
            bp = ax.boxplot([v], positions=[i], widths=0.5, showfliers=False, patch_artist=True,
                            medianprops=dict(color=INK, lw=1.6), whiskerprops=dict(color=INK2, lw=1),
                            capprops=dict(color=INK2, lw=1), boxprops=dict(edgecolor=INK2, lw=1))
            bp["boxes"][0].set_facecolor(COL[g]); bp["boxes"][0].set_alpha(0.28)
            ax.scatter(i + rng.uniform(-0.15, 0.15, len(v)), v, s=14, color=COL[g], edgecolor="white", lw=0.5, zorder=3)
        p = mannwhitneyu(*vals).pvalue if all(len(v) for v in vals) else np.nan
        ax.set_yscale(scale)
        ax.set_xticks(range(2), [f"{LABEL[g]}\n(n = {len(v)})" for g, v in zip(GROUPS, vals)], color=INK)
        ax.set_xlim(-0.6, 1.6)
        ax.set_ylabel(ylab, color=INK)
        ax.set_title(f"Mann-Whitney {fmt_p(p)}", fontsize=10, color=INK2)
        ax.grid(axis="y", color="#e6e5e1", lw=0.6); ax.set_axisbelow(True)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        ax.tick_params(colors=INK2)
        stats.append(dict(panel=title, variable=var, **{f"n {LABEL[g]}": len(v) for g, v in zip(GROUPS, vals)},
                          **{f"median {LABEL[g]}": np.median(v) if len(v) else np.nan for g, v in zip(GROUPS, vals)},
                          **{f"IQR {LABEL[g]}": f"{np.percentile(v, 25):.2f}–{np.percentile(v, 75):.2f}" if len(v) else "" for g, v in zip(GROUPS, vals)},
                          mann_whitney_p=p))
    fig.suptitle(title, fontsize=12, color=INK)
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"), dpi=300)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)
    return stats


all_stats, breakdown, per_patient = [], [], []
for name, wb in COHORTS.items():
    coh = pd.read_excel(wb, "Cohort")
    b = baseline(coh)
    d = coh.set_index("Patient Case ID").join(b)
    d["ROI"] = np.where(d["Kidney ROI (first visit, allowed size)"], "kidney ROI", "no kidney ROI")
    per_patient.append(d.assign(cohort=name))
    bd = pd.crosstab([d["Slope Group"], d["Native/Transplant"]], d["ROI"], margins=True, margins_name="Total")
    breakdown.append((name, bd))
    tag = name.replace("≥ ", "min").replace(" yr", "yr").replace(".", "p")
    for subset, dd in [("all patients", d), ("kidney-ROI patients", d[d.ROI == "kidney ROI"])]:
        t = (f"Case 8 (expanded RF, 2.5 MHz), eGFR follow-up {name} — {subset} (n = {len(dd)})")
        f = RES / f"baseline_egfr_acr_kfre_by_slope_{tag}_{'all' if subset.startswith('all') else 'kidneyROI'}"
        all_stats += plot(dd, t, f)
    print(f"\n===== {name}\n{bd}")
    miss = d[["egfr", "acr", "kfre2"]].isna().sum()
    print("missing baseline (all patients):", miss.to_dict())

st = pd.DataFrame(all_stats)
with pd.ExcelWriter(RES / "baseline_egfr_acr_kfre_by_slope.xlsx") as w:
    pd.DataFrame({"Item": ["Baseline rule", "eGFR", "ACR", "KFRE", "Test", "Script"],
                  "Value": ["closest to imaging within +-90 d (ties earlier); transplant: values before transplant + 3 months "
                            "replaced by first value on/after the cutoff if within 90 d (kfre_tx.py rule)",
                            "PTX BW estGFRckd_epi_2021 (transplant); native CSV CKD-EPI 2021 (native)",
                            "PTX URINE MACRRandom (>0, else CalcACRLevel) mg/mmol (transplant); native CSV Urine ACR mg/mmol (native)",
                            "Tangri 4-variable 2-yr, S0 0.9750; age = VisitReview visit_age closest to imaging; sex = DemographicBaseline",
                            "two-sided Mann-Whitney U", "Scripts/case8_min1p5yr_expanded_rf/baseline_boxplots.py"]}).to_excel(w, sheet_name="Read_Me", index=False)
    st.to_excel(w, sheet_name="Stats", index=False)
    for name, bd in breakdown:
        bd.to_excel(w, sheet_name=f"Breakdown {name.replace('≥ ', 'min ')}")
    pd.concat(per_patient).reset_index().rename(columns={"index": "Patient Case ID"})[
        ["cohort", "Patient Case ID", "Native/Transplant", "Imaging Date", "LMM eGFR Slope (mL/min/1.73m2/year)", "Slope Group",
         "ROI", "egfr", "egfr_gap", "acr", "acr_gap", "age", "male", "kfre2"]].to_excel(w, sheet_name="Per patient", index=False)
pd.set_option("display.width", 250)
print("\n", st.round(3).to_string(index=False))

"""Kaplan-Meier: time from R-score scan to (confirmed 40% eGFR decline or RRT),
stratified by R-fibrosis (SI + H-scan top-3 radiomics) predicted ci 0 vs ci > 0,
with biopsy ci 0 vs ci > 0 as a reference panel.

Endpoint
  baseline  = first eGFR on/after the scan date (same-day allowed)
  40% decline (confirmed) = eGFR <= 60% of baseline, with another value <= 60% at least
              30 days later and no value > 60% in between; event date = first low value
  RRT       = dialysis start (ProspBx_eGFR_RRT_Est), graft failure, or any new transplant after baseline
  censoring = last eGFR date
Patients on RRT at baseline (HD class, 'Baseline RRT' note) are excluded.
"""
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
OUT = ROOT / "Results_Rough" / "km_rfibrosis_ci"
OUT.mkdir(parents=True, exist_ok=True)
CONFIRM_DAYS = 30

# ---------------- cohort: R-score + visit ----------------
rs = pd.read_csv(RS / "outer_predictions_per_sample_mean.csv")
coh = pd.read_excel(RS / "parent_context" / "cohort_source_firstVisit_174.xlsx")
coh = coh.rename(columns={"REDCap Scan ID": "sample_id", "Participant ID (POCUS)": "pid",
                          "Visit date": "scan_date", "Kidney type": "ktype"})
coh["sample_id"] = coh.sample_id.str.strip()  # 7 IDs have a leading space in the source sheet
coh["pid"] = coh.pid.str.strip()
df = rs.merge(coh[["sample_id", "pid", "scan_date", "ktype", "Case type"]], on="sample_id", how="left")
assert df.pid.notna().all(), "unmatched sample_id"
df["scan_date"] = pd.to_datetime(df.scan_date)
df["p_ci_gt0"] = 1 - df.probability_0
df["R_group"] = np.where(df.p_ci_gt0 >= 0.5, "R-fibrosis ci > 0", "R-fibrosis ci 0")
df["ci_group"] = np.where(df.true_index >= 1, "Biopsy ci > 0", "Biopsy ci 0")

# ---------------- eGFR ----------------
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

# ---------------- RRT sources ----------------
ps = pd.read_excel(D / "ProspBx_eGFR_RRT_Est_26-08-26.xlsx", "Patient_Summary")
ps["rrt_date"] = pd.to_datetime(ps["RRT_initiation_date (dd/mm/yy)"], errors="coerce")
ps = ps.set_index("patient_id")
tx = pd.read_csv(D / "PROSPECTIVEBIOPSYFIB-ExportTransplantHist_DATA_LABELS_2026-11-13_1009.csv", encoding="latin1")
tx["txd"] = pd.to_datetime(tx["Date of Transplant"], dayfirst=True, errors="coerce")
tx["fail"] = pd.to_datetime(tx["Transplant Failed Date"], dayfirst=True, errors="coerce")
txdates = tx.dropna(subset=["txd"]).groupby("Participant ID:").txd.apply(list)
txfail = tx.dropna(subset=["fail"]).groupby("Participant ID:").fail.apply(list)


def endpoint(r):
    e = eg[(eg.pid == r.pid) & (eg.src == r.ktype)]
    if r.ktype == "Transplant":  # graft in place at scan only
        g = e[e.txdate <= r.scan_date].txdate.max()
        e = e[e.txdate == g]
    e = e.groupby("date", as_index=False).egfr.mean().sort_values("date")
    post = e[e.date >= r.scan_date]
    out = dict(n_egfr=len(post))
    if r["Case type"] == "hemodialysis":
        return {**out, "status": "excluded: HD at scan"}
    note = ps.Notes.get(r.pid) if r.pid in ps.index else None
    if note == "Baseline RRT":
        return {**out, "status": "excluded: RRT at baseline"}
    if len(post) == 0:
        return {**out, "status": "excluded: no eGFR after scan"}
    b = post.iloc[0]
    fu = post.iloc[1:].reset_index(drop=True)
    thr = 0.6 * b.egfr
    # RRT
    cands = []
    rd = ps.rrt_date.get(r.pid) if r.pid in ps.index else pd.NaT
    if pd.notna(rd):
        if rd <= b.date:
            return {**out, "status": "excluded: RRT at baseline"}
        cands.append(rd)
    cands += [t for t in txdates.get(r.pid, []) if t > b.date]
    if r.ktype == "Transplant":
        cands += [t for t in txfail.get(r.pid, []) if t > b.date]
    rrt = min(cands) if cands else pd.NaT
    # confirmed 40% decline
    dec = pd.NaT
    for i in fu.index[fu.egfr <= thr]:
        later = fu.loc[i + 1:]
        above = later[later.egfr > thr]
        stop = above.date.iloc[0] if len(above) else pd.Timestamp.max
        ok = later[(later.egfr <= thr) & (later.date >= fu.date[i] + pd.Timedelta(days=CONFIRM_DAYS)) & (later.date < stop)]
        if len(ok):
            dec = fu.date[i]
            break
    last = post.date.iloc[-1]
    if len(fu) == 0 and pd.isna(rrt):
        return {**out, "status": "excluded: no follow-up eGFR", "baseline_date": b.date, "baseline_egfr": b.egfr}
    ev_dates = [x for x in [dec, rrt] if pd.notna(x)]
    event = bool(ev_dates)
    end = min(ev_dates) if event else last
    return {**out, "status": "analysed", "baseline_date": b.date, "baseline_egfr": round(b.egfr, 1),
            "last_egfr_date": last, "decline40_date": dec, "rrt_date": rrt,
            "event": int(event), "time_days": (end - r.scan_date).days}


res = pd.concat([df, df.apply(lambda r: pd.Series(endpoint(r)), axis=1)], axis=1)
res.to_csv(OUT / "km_rfibrosis_ci_patient_level.csv", index=False)
a = res[res.status == "analysed"].copy()
a["event"] = a.event.astype(int)
a["time_yr"] = a.time_days.astype(float) / 365.25


# ---------------- KM + log-rank ----------------
def km(t, e):
    t, e = np.asarray(t, float), np.asarray(e, int)
    ts = np.unique(t[e == 1])
    s, S, T = 1.0, [1.0], [0.0]
    var, V = 0.0, [0.0]
    for u in ts:
        n = (t >= u).sum(); d = ((t == u) & (e == 1)).sum()
        s *= 1 - d / n
        var += d / (n * (n - d)) if n > d else 0
        T.append(u); S.append(s); V.append(var)
    if t.max() > T[-1]:  # extend the curve to the last follow-up time
        T.append(t.max()); S.append(s); V.append(var)
    S, V = np.array(S), np.array(V)
    se = S * np.sqrt(V)
    return np.array(T), S, np.clip(S - 1.96 * se, 0, 1), np.clip(S + 1.96 * se, 0, 1)


def logrank(t, e, g):
    t, e, g = np.asarray(t, float), np.asarray(e, int), np.asarray(g)
    lv = np.unique(g); g1 = g == lv[0]
    O = E = Vv = 0.0
    for u in np.unique(t[e == 1]):
        at = t >= u; n = at.sum(); n1 = (at & g1).sum()
        d = ((t == u) & (e == 1)).sum(); d1 = ((t == u) & (e == 1) & g1).sum()
        O += d1; E += d * n1 / n
        if n > 1: Vv += d * (n1 / n) * (1 - n1 / n) * (n - d) / (n - 1)
    stat = (O - E) ** 2 / Vv
    return stat, chi2.sf(stat, 1)


COL = {0: "#2b6cb0", 1: "#c53030"}
panels = [("R_group", ["R-fibrosis ci 0", "R-fibrosis ci > 0"],
           "R-fibrosis score (SI + H-scan, top-3)\npredicted ci 0 vs ci > 0 [P(ci>0) ≥ 0.5]"),
          ("ci_group", ["Biopsy ci 0", "Biopsy ci > 0"], "Reference: biopsy ci 0 vs ci > 0")]
tmax = np.ceil(a.time_yr.max() * 2) / 2
ticks = np.arange(0, tmax + 0.01, 0.5)
fig, axes = plt.subplots(2, 2, figsize=(13, 7.2), gridspec_kw={"height_ratios": [4, 1]})
summary = []
for j, (col, levels, title) in enumerate(panels):
    ax, axt = axes[0, j], axes[1, j]
    for k, lv in enumerate(levels):
        s = a[a[col] == lv]
        T, S, lo, hi = km(s.time_yr, s.event)
        ax.step(T, S, where="post", color=COL[k], lw=2, label=f"{lv} (n={len(s)}, events={s.event.sum()})")
        ax.fill_between(T, lo, hi, step="post", color=COL[k], alpha=0.12)
        cz = s[s.event == 0].time_yr
        ax.plot(cz, [S[np.searchsorted(T, c, side="right") - 1] for c in cz], "|", color=COL[k], ms=7)
        axt.text(-0.06, 0.75 - 0.45 * k, lv.replace("R-fibrosis ", "").replace("Biopsy ", ""), color=COL[k],
                 ha="right", va="center", transform=axt.transAxes, fontsize=9)
        for x in ticks:
            axt.text(x, 0.75 - 0.45 * k, str((s.time_yr >= x).sum()), ha="center", va="center", fontsize=9, color=COL[k])
        for yr in (1, 2):
            i = np.searchsorted(T, yr, side="right") - 1
            summary.append(dict(grouping=col, group=lv, n=len(s), events=int(s.event.sum()),
                                year=yr, event_free=S[i], ci_lo=lo[i], ci_hi=hi[i], at_risk=int((s.time_yr >= yr).sum())))
    stat, p = logrank(a.time_yr, a.event, a[col])
    ax.text(0.03, 0.05, f"Log-rank p = {p:.3g}", transform=ax.transAxes, fontsize=10)
    ax.set_title(title, fontsize=11)
    ax.set_xlim(0, tmax); ax.set_ylim(0, 1.02); ax.set_xticks(ticks)
    ax.set_ylabel("Free of 40% eGFR decline or RRT"); ax.set_xlabel("Years from scan")
    ax.legend(loc="lower left", bbox_to_anchor=(0, 0.1), fontsize=9, frameon=False)
    ax.grid(alpha=0.25)
    axt.set_xlim(0, tmax); axt.set_ylim(0, 1); axt.axis("off")
    axt.set_title("Number at risk", fontsize=9, loc="left")
    summary.append(dict(grouping=col, group="log-rank", n=len(a), events=int(a.event.sum()), year=None,
                        event_free=None, ci_lo=None, ci_hi=None, at_risk=None, chi2=stat, p=p))
fig.suptitle(f"Composite: confirmed 40% eGFR decline (2 values ≥{CONFIRM_DAYS} d apart) or RRT — n = {len(a)}", fontsize=12)
fig.tight_layout()
fig.savefig(OUT / "km_rfibrosis_ci0_vs_gt0.png", dpi=300)
fig.savefig(OUT / "km_rfibrosis_ci0_vs_gt0.pdf")
pd.DataFrame(summary).to_csv(OUT / "km_summary.csv", index=False)

# ---------------- console summary ----------------
print(res.status.value_counts().to_string())
print("\nAnalysed by kidney:", a.ktype.value_counts().to_dict())
print(pd.crosstab(a.ci_group, a.R_group, margins=True))
for col, levels, _ in panels:
    for lv in levels:
        s = a[a[col] == lv]
        print(f"{lv:20s} n={len(s):3d} events={s.event.sum():2d} (decline {s.decline40_date.notna().sum()}, RRT {s.rrt_date.notna().sum()})  median FU {s.time_yr.median():.2f} y")
print(pd.DataFrame(summary).round(3).to_string())

"""Case 8, >= 1.5 yr eGFR follow-up, with the RF-size rule expanded (still 2.5 MHz), LMM slope, slope < -3 vs >= -3.

eGFR values come ONLY from
  All_csv/PROSPBX_DataRequest_10Oct2026.xlsx (PTX BW, estGFRckd_epi_2021)            -> transplant patients
  All_csv/PROSPECTIVEBIOPSYFIB-DataEntryNativeKidne_DATA_LABELS_2026-10-22_1519.csv   -> native patients
and are rebuilt here with the Case 8 eGFR rules:
  on/after first imaging date; transplant: >= transplant date + 3 calendar months; same-day duplicates -> earliest
  (file order); consecutive-day cluster rule = a run of >= 3 consecutive calendar days keeps only its first value
  (runs of 2 are kept); eGFR after RRT initiation removed (ProspBx_eGFR_RRT_Est_26-08-26.xlsx); eGFR on/after a
  follow-up transplant removed (TransplantHist); >= 3 eGFR.
  The rebuild is checked to reproduce case8_no_min_followup_RRT_adjusted_RF2928or2944x192_egfr.xlsx exactly
  (all 374 patients, dates and values) before anything else is done.
Patient-level Case 8 exclusions upstream of the RF rule (no eGFR, < 3 eGFR, early post-transplant scan, Baseline RRT)
are inherited from that workbook: the patients considered are its 374 + the 78 it removed for RF size.
RF rule (expanded): >= 1 image at the first imaging visit (patient_inventory_created_Y2026M05D04, the inventory Case 8
  used) at 2.5 MHz with size in ALLOWED.
Fitting rules: last - first eGFR date >= 1.5 * 365.25 d; >= 3 eGFR >= 30 d apart (greedy from first).
LMM: egfr ~ 1 + t + (1 + t | patient), ML, t = years from scan, fitted once on the final cohort;
     patient slope = fixed + random slope.
Kidney ROI: >= 1 first-visit image of an allowed size with a kidney ROI (rscore_allpatients_Data1103/images_audit.csv).
"""
import warnings
from pathlib import Path
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

ROOT = Path(__file__).resolve().parents[2]
D = ROOT / "All_csv"
CV = Path(r"\\Tupil\tupil data\Hikaru\POCUS_Kidney_Analysis\Analysis_scripts\csv\Clinical Variables")
CASE8 = CV / "case8_no_min_followup_RRT_adjusted_RF2928or2944x192_egfr.xlsx"
CASE8_25 = CV / "case8_no_min_followup_RRT_adjusted_2p5MHz_only_egfr_with_RRT_timing.xlsx"  # patient metadata only
INV = Path(r"\\Tupil\tupil data\Hikaru\Data_2026_1103\#zMASTER_DATA\patient_inventory_created_Y2026M05D04_H08M44S53.xlsx")
AUDIT = Path(r"\\Tupil\tupil data\Hikaru\POCUS_Kidney_Analysis\Analysis_scripts\Redo_analysis_1106_2026"
             r"\rscore_allpatients_Data1103\images_audit.csv")
OUT = ROOT / "Results_Rough" / "case8_min1p5yr_expanded_rf"
OUT.mkdir(parents=True, exist_ok=True)

ALLOWED = {"2928 x 192", "2944 x 192", "2784 x 176", "2960 x 192", "2960 x 184", "2944 x 176",
           "3488 x 184", "2736 x 192", "2928 x 168", "2912 x 160"}
FREQ = 2_500_000
MIN_FU_DAYS = 1.5 * 365.25
MIN_GAP, MIN_SPACED = 30, 3
SLOPE_CUT = -3.0
T = "time_years_from_scan"

# ---------------- source eGFR ----------------
bw = pd.read_excel(D / "PROSPBX_DataRequest_10Oct2026.xlsx", "PTX BW")
nk = pd.read_csv(D / "PROSPECTIVEBIOPSYFIB-DataEntryNativeKidne_DATA_LABELS_2026-10-22_1519.csv", encoding="latin1")
ecol = [c for c in nk.columns if c.startswith("Predicted eGFR (CKD-EPI")][0]
src = pd.concat([
    pd.DataFrame({"pid": bw["StudyIDRequestList::StudyID"], "date": pd.to_datetime(bw["Date"]),
                  "egfr": bw["estGFRckd_epi_2021"], "kt": "transplant", "row": np.arange(len(bw)) + 2,
                  "source": "PROSPBX_DataRequest_10Oct2026.xlsx: PTX BW"}),
    pd.DataFrame({"pid": nk["Participant ID:"],
                  "date": pd.to_datetime(nk["Date of Lab Results:"], dayfirst=True, errors="coerce"),
                  "egfr": pd.to_numeric(nk[ecol], errors="coerce"), "kt": "native", "row": np.arange(len(nk)) + 2,
                  "source": "PROSPECTIVEBIOPSYFIB-DataEntryNativeKidne_DATA_LABELS_2026-10-22_1519.csv"}),
]).dropna(subset=["date", "egfr"])

rrt_est = pd.read_excel(D / "ProspBx_eGFR_RRT_Est_26-08-26.xlsx", "Patient_Summary").set_index("patient_id")
rrt_est["rrt"] = pd.to_datetime(rrt_est["RRT_initiation_date (dd/mm/yy)"], errors="coerce")
txh = pd.read_csv(D / "PROSPECTIVEBIOPSYFIB-ExportTransplantHist_DATA_LABELS_2026-11-13_1009.csv", encoding="latin1")
txh["txd"] = pd.to_datetime(txh["Date of Transplant"], dayfirst=True, errors="coerce")

# ---------------- patients considered: Case 8 374 + 78 RF-size exclusions ----------------
c8_ps = pd.read_excel(CASE8, "Patient_Summary").set_index("patient_id")
c8_eg = pd.read_excel(CASE8, "Filtered_eGFR")
c8_eg["egfr_date"] = pd.to_datetime(c8_eg.egfr_date)
c8_ex = pd.read_excel(CASE8, "Excluded_Patients")
rf_ex = c8_ex.loc[c8_ex.Status.str.startswith("no first-visit RF"), "patient_id"].tolist()
meta25 = pd.read_excel(CASE8_25, "Patient_Summary").set_index("patient_id")
meta = pd.concat([c8_ps[["native_transplant_status", "transplant_date", "first_imaging_scan_date"]],
                  meta25.reindex(rf_ex)[["native_transplant_status", "transplant_date", "first_imaging_scan_date"]]])
meta["source_set"] = ["case8 (374)"] * len(c8_ps) + ["case8 RF-size excluded (78)"] * len(rf_ex)
# patients of the 78 with no row in the 2.5 MHz workbook: scan date + kidney type from the inventory, transplant date from
# TransplantHist (latest transplant on/before the scan)
_inv = pd.read_excel(INV, "Patient Inventory (Flat)")
_inv["patient_id"] = _inv["Patient ID"].map(lambda n: f"REB16236-{int(n):03d}")
_inv["Scan Date"] = pd.to_datetime(_inv["Scan Date"])
_first = _inv.sort_values("Scan Date").groupby("patient_id").first()
for p in meta.index[meta.first_imaging_scan_date.isna()]:
    if p not in _first.index:
        continue
    scan, vt = _first.loc[p, "Scan Date"], str(_first.loc[p, "Visit Type"])
    meta.loc[p, "first_imaging_scan_date"] = scan
    if vt.startswith("PT"):
        txs = txh.loc[(txh["Participant ID:"] == p) & (txh.txd <= scan), "txd"]
        meta.loc[p, ["native_transplant_status", "transplant_date"]] = ["transplant", txs.max()]
    else:
        meta.loc[p, "native_transplant_status"] = "native"
    meta.loc[p, "source_set"] = "case8 RF-size excluded (78); metadata from inventory + TransplantHist"


def build_egfr(p, r):
    s = src[(src.pid == p) & (src.kt == r.native_transplant_status)]
    scan = pd.to_datetime(r.first_imaging_scan_date)
    s = s[s.date >= scan]
    if r.native_transplant_status == "transplant":
        s = s[s.date >= pd.to_datetime(r.transplant_date) + pd.DateOffset(months=3)]
    s = s.sort_values(["date", "row"]).drop_duplicates("date", keep="first")
    run = (s.date.diff().dt.days.fillna(999) > 1).cumsum()
    s = s[(s.groupby(run).date.transform("size") < 3) | (s.groupby(run).cumcount() == 0)]
    rrt = rrt_est.rrt.get(p, pd.NaT)
    if pd.notna(rrt):
        s = s[s.date <= rrt]
    fut = [t for t in txh.loc[txh["Participant ID:"] == p, "txd"] if pd.notna(t) and t > scan]
    fut = [t for t in fut if len(s) and t <= s.date.max()]
    if fut:
        s = s[s.date < min(fut)]
    return s.assign(patient_id=p, native_transplant_status=r.native_transplant_status, scan_date=scan,
                    **{T: (s.date - scan).dt.days / 365.25})


meta_ok = meta.dropna(subset=["first_imaging_scan_date"])
eg = pd.concat([build_egfr(p, r) for p, r in meta_ok.iterrows()], ignore_index=True)

# ---- reproduction check against Case 8 (374)
mism = []
for p in c8_ps.index:
    a = eg[eg.patient_id == p].set_index("date").egfr
    b = c8_eg[c8_eg.patient_id == p].set_index("egfr_date").egfr
    if not (a.index.equals(b.index) and np.allclose(a.values, b.values)):
        mism.append(p)
assert not mism, f"rebuild differs from Case 8 for {mism}"
print(f"Rebuild reproduces Case 8 exactly for all {len(c8_ps)} patients.")

# ---------------- RF rule (expanded sizes, 2.5 MHz, first visit) ----------------
inv = pd.read_excel(INV, "Patient Inventory (Flat)")
inv["patient_id"] = inv["Patient ID"].map(lambda n: f"REB16236-{int(n):03d}")
inv["Scan Date"] = pd.to_datetime(inv["Scan Date"])
inv["size"] = inv["Image Samples Depth"].astype("Int64").astype(str) + " x " + inv["Image RF Lines"].astype("Int64").astype(str)
fv = inv[inv["Scan Date"] == inv.groupby("patient_id")["Scan Date"].transform("min")]
fv = fv.assign(ok=fv["size"].isin(ALLOWED) & (fv["Image Transmit Frequency (Hz)"] == FREQ))
rf_ok = fv.groupby("patient_id").ok.any()
rf_sizes = fv[fv.ok].groupby("patient_id")["size"].agg(lambda s: ", ".join(sorted(set(s))))

# ---------------- kidney ROI ----------------
aud = pd.read_csv(AUDIT, parse_dates=["Imaging Date"])
aud = aud.merge(meta_ok["first_imaging_scan_date"].pipe(pd.to_datetime).rename("scan"), left_on="Patient Case ID",
                right_index=True)
aud = aud[(aud["Imaging Date"] == aud.scan) & aud["RF Size"].isin(ALLOWED)]
roi = aud.groupby("Patient Case ID")["Kidney ROI"].any()
n_roi_img = aud[aud["Kidney ROI"]].groupby("Patient Case ID").size()


# ---------------- per patient ----------------
def spaced(dates):
    n, last = 0, None
    for d in sorted(dates):
        if last is None or (d - last).days >= MIN_GAP:
            n, last = n + 1, d
    return n


per = eg.groupby("patient_id").agg(**{
    "Native/Transplant": ("native_transplant_status", "first"), "Imaging Date": ("scan_date", "first"),
    "eGFR Measurements": ("egfr", "size"), "First eGFR Date": ("date", "min"), "Last eGFR Date": ("date", "max")})
per = per.reindex(meta.index)
per["Source set"] = meta.source_set
per["eGFR Duration (days)"] = (per["Last eGFR Date"] - per["First eGFR Date"]).dt.days
per["Follow-up Years"] = per["eGFR Duration (days)"] / 365.25
per["eGFR >= 30 d apart"] = eg.groupby("patient_id").date.apply(spaced).reindex(per.index)
per["2.5 MHz image of allowed size at first visit"] = rf_ok.reindex(per.index).fillna(False).astype(bool)
per["Qualifying RF sizes"] = rf_sizes.reindex(per.index)
per["Kidney ROI (first visit, allowed size)"] = roi.reindex(per.index).fillna(False).astype(bool)
per["Kidney-ROI images"] = n_roi_img.reindex(per.index).fillna(0).astype(int)


S = "LMM eGFR Slope (mL/min/1.73m2/year)"
COLS = ["Native/Transplant", "Source set", "Imaging Date", S, "Slope Group", "eGFR Measurements", "eGFR >= 30 d apart",
        "Follow-up Years", "First eGFR Date", "Last eGFR Date", "Qualifying RF sizes",
        "Kidney ROI (first visit, allowed size)", "Kidney-ROI images"]
NO_RF = "no 2.5 MHz first-visit image of an allowed size"


def run(min_years):
    """Select the >= min_years cohort, fit the LMM on it once, write the slope-group x ROI workbook."""
    tag = f"{min_years:g}yr".replace(".", "p")
    p = per.copy()

    def status(r):
        if not r["2.5 MHz image of allowed size at first visit"]:
            return NO_RF
        if pd.isna(r["eGFR Measurements"]) or r["eGFR Measurements"] < 3:
            return "fewer than 3 eGFR (rebuild)"
        why = []
        if r["eGFR Duration (days)"] < min_years * 365.25:
            why.append(f"eGFR follow-up < {min_years:g} yr")
        if r["eGFR >= 30 d apart"] < MIN_SPACED:
            why.append("< 3 eGFR spaced >= 30 d")
        return "; ".join(why) or "included"

    p["Status"] = p.apply(status, axis=1)
    keep = p.index[p.Status == "included"]
    fit = eg[eg.patient_id.isin(keep)].copy()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        m = smf.mixedlm(f"egfr ~ {T}", fit, groups=fit.patient_id, re_formula=f"~{T}").fit(reml=False)
    re_ = pd.DataFrame(m.random_effects).T
    p[S] = (m.fe_params[T] + re_[T]).reindex(p.index)
    p["Slope Group"] = np.where(p[S] < SLOPE_CUT, "< -3", ">= -3")
    p.loc[p[S].isna(), "Slope Group"] = None
    coh = p.loc[keep].copy()

    coh["ROI"] = np.where(coh["Kidney ROI (first visit, allowed size)"], "kidney ROI", "no kidney ROI")
    bd = pd.crosstab(coh["Slope Group"], coh["ROI"], margins=True, margins_name="Total")
    bd_kt = pd.crosstab([coh["Slope Group"], coh["Native/Transplant"]], coh["ROI"], margins=True, margins_name="Total")
    rf = p["2.5 MHz image of allowed size at first visit"]
    is374 = p["Source set"] == "case8 (374)"
    flow = pd.DataFrame([
        ("Case 8 patients before the RF rule (374 kept + 78 removed for RF size)", len(p)),
        ("  removed: " + NO_RF, int((p.Status == NO_RF).sum())),
        ("  passing the expanded RF rule", int(rf.sum())),
        ("    of which in the Case 8 374", int((rf & is374).sum())),
        ("    of which newly added by the expanded sizes", int((rf & ~is374).sum())),
        ("  removed: fewer than 3 eGFR (rebuild)", int((p.Status == "fewer than 3 eGFR (rebuild)").sum())),
        (f"  removed: eGFR follow-up < {min_years:g} yr only", int((p.Status == f"eGFR follow-up < {min_years:g} yr").sum())),
        ("  removed: < 3 eGFR >= 30 d apart only", int((p.Status == "< 3 eGFR spaced >= 30 d").sum())),
        ("  removed: both", int(p.Status.str.contains("; ").sum())),
        (f"Final >= {min_years:g} yr cohort", len(coh)),
        ("  transplant / native", f"{(coh['Native/Transplant'] == 'transplant').sum()} / {(coh['Native/Transplant'] == 'native').sum()}"),
        ("  with kidney ROI", int(coh["Kidney ROI (first visit, allowed size)"].sum())),
    ], columns=["Step", "n"])
    model = pd.DataFrame({"Item": ["Formula", "Fitted on (patients / eGFR)", "Converged", "Fixed slope", "Fixed slope SE",
                                   "Fixed slope p", "Fixed intercept"],
                          "Value": ["egfr ~ 1 + t + (1 + t | patient), ML", f"{len(coh)} / {len(fit)}", m.converged,
                                    round(m.fe_params[T], 4), round(m.bse[T], 4), float(m.pvalues[T]),
                                    round(m.fe_params["Intercept"], 3)]})
    readme = pd.DataFrame({"Item": ["eGFR sources", "Rules", "Minimum follow-up", "Allowed RF sizes (2.5 MHz)",
                                    "Slope groups", "Script"],
                           "Value": ["PROSPBX_DataRequest_10Oct2026.xlsx (PTX BW); "
                                     "PROSPECTIVEBIOPSYFIB-DataEntryNativeKidne_DATA_LABELS_2026-10-22_1519.csv",
                                     "Case 8 eGFR rules rebuilt from the two files and verified against the Case 8 "
                                     "workbook (374/374 exact); >= 3 eGFR >= 30 d apart",
                                     f"last - first eGFR date >= {min_years:g} x 365.25 d",
                                     ", ".join(sorted(ALLOWED)),
                                     f"LMM slope (fitted on this cohort) < {SLOPE_CUT} vs >= {SLOPE_CUT}",
                                     "Scripts/case8_min1p5yr_expanded_rf/case8_min1p5yr_expanded_rf.py"]})
    dst = OUT / f"case8_min{tag}_expanded_rf_slope_groups.xlsx"
    with pd.ExcelWriter(dst) as w:
        readme.to_excel(w, sheet_name="Read_Me", index=False)
        flow.to_excel(w, sheet_name="Flow", index=False)
        bd.to_excel(w, sheet_name="Slope x ROI")
        bd_kt.to_excel(w, sheet_name="Slope x type x ROI")
        coh[COLS].reset_index().rename(columns={"index": "Patient Case ID", "patient_id": "Patient Case ID"}) \
            .sort_values("Patient Case ID").to_excel(w, sheet_name="Cohort", index=False)
        p[p.Status != "included"][["Native/Transplant", "Source set", "Status", "eGFR Measurements",
                                   "eGFR >= 30 d apart", "Follow-up Years"]].reset_index() \
            .to_excel(w, sheet_name="Excluded", index=False)
        model.to_excel(w, sheet_name="LMM_Model", index=False)
        fit.rename(columns={"date": "egfr_date"})[["patient_id", "native_transplant_status", "scan_date", "egfr_date",
                                                   "egfr", T, "source", "row"]].to_excel(w, sheet_name="eGFR_Used", index=False)

    print(f"\n==================== >= {min_years:g} yr  ->  {dst.name}")
    print(flow.to_string(index=False))
    print(model.to_string(index=False))
    print("\nSlope group x kidney ROI:\n", bd)
    print("\nby type:\n", bd_kt)
    new = coh[coh["Source set"] != "case8 (374)"]
    print("\nExpanded-size patients in cohort:")
    print(new[["Qualifying RF sizes", "Follow-up Years", S, "Slope Group", "Kidney ROI (first visit, allowed size)"]].round(2).to_string())
    print("No-ROI patients:", sorted(coh.index[~coh["Kidney ROI (first visit, allowed size)"]]))


for yrs in (1.0, 1.5):
    run(yrs)

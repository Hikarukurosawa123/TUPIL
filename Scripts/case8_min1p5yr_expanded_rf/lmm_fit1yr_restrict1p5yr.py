"""LMM slope fitted on the >= 1 yr cohort, then the cohort restricted to >= 1.5 yr follow-up.

Step 1: LMM egfr ~ 1 + t + (1 + t | patient), ML, t = years from scan, fitted on ALL eGFR of the >= 1 yr cohort
        (expanded RF rule, >= 3 eGFR, >= 3 eGFR >= 30 d apart, last - first eGFR >= 1 x 365.25 d).
        Data = eGFR_Used sheet of case8_min1yr_expanded_rf_slope_groups.xlsx (written by
        case8_min1p5yr_expanded_rf.py); the refit is checked against the slopes stored in that workbook.
Step 2: keep patients with last - first eGFR >= 1.5 x 365.25 d; slope = fixed + random slope from the step-1 fit.
Step 3: compare with the slopes from the model fitted on the >= 1.5 yr cohort only
        (case8_min1p5yr_expanded_rf_slope_groups.xlsx).
"""
import warnings
from pathlib import Path
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "Results_Rough" / "case8_min1p5yr_expanded_rf"
W1 = OUT / "case8_min1yr_expanded_rf_slope_groups.xlsx"
W15 = OUT / "case8_min1p5yr_expanded_rf_slope_groups.xlsx"
T = "time_years_from_scan"
S = "LMM eGFR Slope (mL/min/1.73m2/year)"
SLOPE_CUT = -3.0
MIN_DAYS = 1.5 * 365.25

eg = pd.read_excel(W1, "eGFR_Used")
coh1 = pd.read_excel(W1, "Cohort").set_index("Patient Case ID")
coh15_old = pd.read_excel(W15, "Cohort").set_index("Patient Case ID")

# ---- step 1: fit on the >= 1 yr cohort
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    m = smf.mixedlm(f"egfr ~ {T}", eg, groups=eg.patient_id, re_formula=f"~{T}").fit(reml=False)
re_ = pd.DataFrame(m.random_effects).T
slope = m.fe_params[T] + re_[T]
assert set(slope.index) == set(coh1.index)
diff = (slope.reindex(coh1.index) - coh1[S]).abs().max()
assert diff < 1e-3, f"refit differs from the stored >= 1 yr slopes by {diff}"
print(f"Refit on {eg.patient_id.nunique()} patients / {len(eg)} eGFR reproduces the stored >= 1 yr slopes "
      f"(max |diff| = {diff:.2e}).")

# ---- step 2: restrict to >= 1.5 yr
dur = (pd.to_datetime(coh1["Last eGFR Date"]) - pd.to_datetime(coh1["First eGFR Date"])).dt.days
coh = coh1[dur >= MIN_DAYS].copy()
coh[S] = slope.reindex(coh.index)
coh["Slope Group"] = np.where(coh[S] < SLOPE_CUT, "< -3", ">= -3")
assert set(coh.index) == set(coh15_old.index), "restricted cohort differs from the >= 1.5 yr cohort"
print(f"Restricted to >= 1.5 yr: {len(coh)} patients (same IDs as the >= 1.5 yr workbook).")

# ---- step 3: compare with the slopes fitted on the >= 1.5 yr cohort only
cmp_ = pd.DataFrame({
    "Native/Transplant": coh["Native/Transplant"],
    "Follow-up Years": coh["Follow-up Years"],
    "eGFR Measurements": coh["eGFR Measurements"],
    "Slope (fit on >= 1 yr)": coh[S],
    "Group (fit on >= 1 yr)": coh["Slope Group"],
    "Slope (fit on >= 1.5 yr)": coh15_old[S].reindex(coh.index),
    "Group (fit on >= 1.5 yr)": coh15_old["Slope Group"].reindex(coh.index),
})
cmp_["Slope difference (1 yr fit - 1.5 yr fit)"] = cmp_["Slope (fit on >= 1 yr)"] - cmp_["Slope (fit on >= 1.5 yr)"]
cmp_["Group changed"] = cmp_["Group (fit on >= 1 yr)"] != cmp_["Group (fit on >= 1.5 yr)"]
cmp_.index.name = "Patient Case ID"
cross = pd.crosstab(cmp_["Group (fit on >= 1 yr)"], cmp_["Group (fit on >= 1.5 yr)"], margins=True,
                    margins_name="Total")

coh["ROI"] = np.where(coh["Kidney ROI (first visit, allowed size)"], "kidney ROI", "no kidney ROI")
bd = pd.crosstab(coh["Slope Group"], coh["ROI"], margins=True, margins_name="Total")
bd_kt = pd.crosstab([coh["Slope Group"], coh["Native/Transplant"]], coh["ROI"], margins=True, margins_name="Total")

model = pd.DataFrame({"Item": ["Formula", "Fitted on (patients / eGFR)", "Converged", "Fixed slope", "Fixed slope SE",
                               "Fixed slope p", "Fixed intercept", "Slope SD (random)", "Reported on"],
                      "Value": ["egfr ~ 1 + t + (1 + t | patient), ML, t = years from scan",
                                f"{eg.patient_id.nunique()} / {len(eg)} (>= 1 yr cohort)", m.converged,
                                round(m.fe_params[T], 4), round(m.bse[T], 4), float(m.pvalues[T]),
                                round(m.fe_params["Intercept"], 3), round(float(np.sqrt(m.cov_re.loc[T, T])), 4),
                                f"{len(coh)} patients with last - first eGFR >= 1.5 x 365.25 d"]})
readme = pd.DataFrame({"Item": ["Approach", "Data", "Slope groups", "Script"],
                       "Value": ["LMM fitted on the >= 1 yr cohort (all eGFR); per-patient slope = fixed + random "
                                 "slope; cohort then restricted to >= 1.5 yr follow-up",
                                 f"{W1.name}: eGFR_Used / Cohort (refit verified against stored slopes)",
                                 f"< {SLOPE_CUT} vs >= {SLOPE_CUT}",
                                 "Scripts/case8_min1p5yr_expanded_rf/lmm_fit1yr_restrict1p5yr.py"]})

dst = OUT / "case8_min1p5yr_expanded_rf_slope_groups_fit1yr.xlsx"
with pd.ExcelWriter(dst) as w:
    readme.to_excel(w, sheet_name="Read_Me", index=False)
    model.to_excel(w, sheet_name="LMM_Model", index=False)
    bd.to_excel(w, sheet_name="Slope x ROI")
    bd_kt.to_excel(w, sheet_name="Slope x type x ROI")
    coh.drop(columns="ROI").reset_index().sort_values("Patient Case ID").to_excel(w, sheet_name="Cohort", index=False)
    cmp_.reset_index().sort_values("Patient Case ID").to_excel(w, sheet_name="Compare_vs_1p5yr_fit", index=False)
    cross.to_excel(w, sheet_name="Group_crosstab")
    eg[eg.patient_id.isin(coh.index)].to_excel(w, sheet_name="eGFR_reported_patients", index=False)

print(f"\n-> {dst.name}")
print(model.to_string(index=False))
print("\nSlope group x kidney ROI:\n", bd)
print("\nby type:\n", bd_kt)
print("\nGroup (fit on >= 1 yr) vs group (fit on >= 1.5 yr):\n", cross)
d = cmp_["Slope difference (1 yr fit - 1.5 yr fit)"]
print(f"\nSlope difference: mean {d.mean():.3f}, median {d.median():.3f}, min {d.min():.3f}, max {d.max():.3f}; "
      f"r = {cmp_['Slope (fit on >= 1 yr)'].corr(cmp_['Slope (fit on >= 1.5 yr)']):.4f}")
print("\nPatients changing group:")
print(cmp_[cmp_["Group changed"]].round(2).to_string())

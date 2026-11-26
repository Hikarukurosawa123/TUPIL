"""Align cached QUS to the designated Case 8 scan and a single validated kidney episode."""
import json
from pathlib import Path
import pickle

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
STAGE1 = ROOT / "Results_Rough" / "egfr_longtitudinal_all"
OUT = STAGE1 / "qus_case8_2yr"
CASE8 = ROOT / "Results_Rough" / "01_eGFR_slope_model3_min1yr" / "case8_LMM_fit_min1yr_3egfr30d_min2yr_subset_19Nov2026.xlsx"
REDO = ROOT.parent / "POCUS_Kidney_Analysis" / "Analysis_scripts" / "Redo_analysis_1106_2026"
CACHE = REDO / "v10_figures" / "egfr_slope" / "13_multimap_filters_binwidth" / "cache"
EXTRA = REDO / "egfr_slope_redo" / "Heat_map_LMM_fit_all_patients"
MAPS = {"SI": "9", "HScan": "3", "MBF": "2.8", "ESD": "1.9", "EAC": "30", "SS": "1.9"}
FILTERS = {"original", "log-sigma-1-0-mm-3D", "log-sigma-2-0-mm-3D", "log-sigma-3-0-mm-3D",
           "wavelet-LL", "wavelet-LH", "wavelet-HL", "wavelet-HH"}
CLASSES = {"firstorder", "glcm", "glrlm", "glszm", "gldm", "ngtdm"}
YEAR = 365.25
STAGE1_COLUMNS = {c: f"stage1_{c}" for c in (
    "index_date", "baseline_egfr_date", "acr_date", "baseline_egfr_source", "baseline_egfr_row",
    "acr_source", "acr_row", "followup", "n_targets",
)}


def radiomics(eligible):
    image_book = EXTRA / "Eligible_Patients_min2yr_RFimage_KidneyROI_19Nov2026.xlsx"
    images = pd.read_excel(image_book, sheet_name="KidneyROI_Images")
    images = images[images["Patient Case ID"].isin(eligible)].copy()
    records = []
    prefixes = tuple(f"{m}|{b}|" for m, b in MAPS.items())
    cache_sources = []
    for pid, num in zip(images["Patient Case ID"], images["Image Num"]):
        filename = f"{pid}__P{int(pid.split('-')[1])}_Image_{int(num)}.pkl"
        path = CACHE / filename
        if not path.exists():
            path = EXTRA / "extra_cache" / filename
        cache_sources.append(str(path))
        with path.open("rb") as stream:
            data = pickle.load(stream)
        record = {"pid": pid}
        for key, value in data.items():
            if key.startswith(prefixes):
                map_name, _, name = key.split("|", 2)
                filt, family, _ = name.split("_", 2)
                if filt in FILTERS and family in CLASSES:
                    record[f"{map_name}__{name}"] = value
        records.append(record)
    features = pd.DataFrame(records).groupby("pid").mean(numeric_only=True)
    images["cache"] = cache_sources
    for modality in MAPS:
        subset = features.filter(regex=f"^{modality}__")
        assert np.isfinite(subset).any(axis=1).all()
    return features.replace([np.inf, -np.inf], np.nan).sort_index(axis=1), images


def prepare(out):
    out.mkdir(parents=True, exist_ok=True)
    assert not (out / "scan_cohort.csv").exists()
    book = pd.read_excel(CASE8, sheet_name=["LMM_Slopes", "eGFR_Used"])
    case = book["LMM_Slopes"].rename(columns={
        "Patient Case ID": "pid", "Imaging Date": "scan", "Follow-up Years": "workbook_followup",
        "Native/Transplant": "workbook_kind",
        "LMM eGFR Slope (mL/min/1.73m2/year)": "workbook_slope",
    }).set_index("pid")
    assert len(case) == 96 and case.workbook_followup.ge(2).all()
    features, images = radiomics(case.index)
    inventory = pd.read_excel(EXTRA / "Eligible_Patients_min2yr_RFimage_19Nov2026.xlsx",
                              sheet_name="Eligible_Images")
    keys = ["Patient Case ID", "Image Num"]
    images = images.merge(inventory[keys + ["Scan Date"]], on=keys, how="left", validate="1:1")
    assert pd.to_datetime(images["Scan Date"]).eq(images["Patient Case ID"].map(case.scan)).all()
    blood = book["eGFR_Used"]
    assert blood.groupby("patient_id").transplant_date.nunique(dropna=False).le(1).all()
    case["workbook_tx"] = blood.groupby("patient_id").transplant_date.first()
    stage = pd.read_csv(STAGE1 / "clinical_baselines.csv", index_col="pid",
                        parse_dates=["index_date", "tx", "next_tx", "graft_failure", "rrt"])
    cohort = stage.reindex(features.index).join(case[["scan", "workbook_tx", "workbook_slope", "workbook_kind"]])
    assert cohort["index_date"].notna().all()
    assert cohort.workbook_kind.eq(cohort.kind).all()
    assert (cohort.tx.eq(cohort.workbook_tx) | (cohort.tx.isna() & cohort.workbook_tx.isna())).all()
    assert (cohort.tx.isna() | cohort.scan.ge(cohort.tx + pd.DateOffset(months=3))).all()
    cohort["age"] += (cohort.scan - cohort.index_date).dt.days / YEAR
    audit = pd.read_csv(STAGE1 / "lab_audit.csv.gz", low_memory=False, parse_dates=["date"])
    labs = audit[audit.disposition.eq("Eligible")].merge(cohort[["series", "scan"]], on=["pid", "series"],
                                                        validate="m:1")
    labs["lag"] = (labs.scan - labs.date).dt.days
    for label, column, days in [("egfr", "baseline_egfr", 365), ("acr", "acr", 365)]:
        pre = (labs[labs.lab.eq(label) & labs.lag.between(0, days)]
               .sort_values("date").drop_duplicates("pid", keep="last").set_index("pid"))
        cohort[column] = pre.value
        cohort[f"scan_{label}_date"] = pre.date
        cohort[f"scan_{label}_row"] = pre.source_row
        cohort[f"scan_{label}_source"] = pre.source
    cohort["log_acr"] = np.log1p(cohort.acr)
    assert cohort.baseline_egfr.notna().all()
    obs = labs[labs.lab.eq("egfr") & labs.lag.lt(0)].copy()
    obs["t"] = -obs.lag / YEAR
    obs = obs.rename(columns={"value": "egfr"})[["pid", "date", "series", "t", "egfr", "source", "source_row"]]
    obs = obs.sort_values(["pid", "date"]).reset_index(drop=True)
    assert not obs.duplicated(["pid", "date"]).any()
    span = obs.groupby("pid").t.agg(["min", "max"])
    cohort["scan_to_last_years"] = span["max"]
    cohort["postscan_span_years"] = span["max"] - span["min"]
    cohort["observations"] = obs.groupby("pid").size()
    assert cohort.scan_to_last_years.ge(2).all()
    assert cohort.observations.ge(3).all()
    check = obs.join(cohort[["scan", "tx", "next_tx", "graft_failure", "rrt"]], on="pid")
    for stop in ("next_tx", "graft_failure", "rrt"):
        assert (check[stop].isna() | check.date.lt(check[stop])).all()
    cohort.rename(columns=STAGE1_COLUMNS).to_csv(out / "scan_cohort.csv")
    obs.to_csv(out / "postscan_outcomes.csv", index=False)
    features.to_csv(out / "qus_features.csv.gz", compression="gzip")
    images.to_csv(out / "image_provenance.csv", index=False)
    exclusions = case.index.difference(cohort.index)
    disposition = pd.DataFrame({"pid": case.index, "included": ~case.index.isin(exclusions)})
    disposition["reason"] = np.where(disposition.included, "Included", "No qualifying kidney ROI")
    disposition.to_csv(out / "cohort_disposition.csv", index=False)
    assignments = pd.read_csv(STAGE1 / "fold_assignments.csv").query("pid in @cohort.index")
    assert assignments.groupby("repetition").pid.nunique().eq(len(cohort)).all()
    assert assignments.groupby("repetition").fold.nunique().eq(5).all()
    assignments.to_csv(out / "fold_assignments.csv", index=False)
    summary = dict(case8_patients=len(case), patients=len(cohort), images=len(images),
                   observations=len(obs), missing_scan_acr=int(cohort.acr.isna().sum()),
                   counts_by_kind=cohort.kind.value_counts().to_dict(),
                   feature_counts={m: int(sum(features.columns.str.startswith(m + "__"))) for m in MAPS},
                   postscan_span_lt2=int(cohort.postscan_span_years.lt(2).sum()))
    (out / "cohort_summary.json").write_text(json.dumps(summary, indent=2))
    print(summary, flush=True)


if __name__ == "__main__":
    prepare(OUT)



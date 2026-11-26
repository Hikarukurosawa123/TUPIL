"""Extract the 01-analysis radiomics (multimap_filters.extract_one: Original + LoG + Wavelet, first-order + texture,
whole kidney minus pelvis) for first-visit kidney-ROI images of the >= 1.5 yr expanded-RF cohort that are NOT in the
shared cache (v10_figures/egfr_slope/13_multimap_filters_binwidth/cache).

Only the bin widths the 01 models use are computed (SI 9, H-scan 3, ESD 1.9, SS 1.9); keys are identical to the shared
cache ('SI|9|original_firstorder_Mean', ...). The shared cache is only read, never written: output goes to
Results_Rough/case8_min1p5yr_expanded_rf/multimap_cache_extra/ (one pickle per image, resumable).
Jobs (raw / ROI / QUS / H-scan paths) come from the noMinFU all-RF job list.
Run: C:/Users/tupil/venvs/radiomics_env/python.exe extract_multimap_missing.py [--workers 8]
"""
import argparse
import os
import sys
from pathlib import Path

import pandas as pd

RED = Path(r"\\Tupil\tupil data\Hikaru\POCUS_Kidney_Analysis\Analysis_scripts\Redo_analysis_1106_2026")
MMDIR = RED / "v10_figures" / "egfr_slope" / "13_multimap_filters_binwidth"
sys.path.insert(0, str(MMDIR))
import multimap_filters as MF  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
RES = ROOT / "Results_Rough" / "case8_min1p5yr_expanded_rf"
OUT_CACHE = RES / "multimap_cache_extra"
SHARED = MMDIR / "cache"
JOBS = RED / "v10_figures" / "patient_groups" / "noMinFU_radiomics" / "radiomics_SI_HScan_allRF_noMinFU_jobs.csv"
AUDIT = RED / "rscore_allpatients_Data1103" / "images_audit.csv"
COHORT = RES / "case8_min1p5yr_expanded_rf_slope_groups.xlsx"
ALLOWED = {"2928 x 192", "2944 x 192", "2784 x 176", "2960 x 192", "2960 x 184", "2944 x 176",
           "3488 x 184", "2736 x 192", "2928 x 168", "2912 x 160"}
BWS = {"SI": [9], "HScan": [3], "ESD": [1.9], "SS": [1.9], "EAC": [], "MBF": [], "Bmode": []}


def extract_to(job, bws, cache_dir):
    MF.CACHE = Path(cache_dir)  # set inside each worker (spawned processes re-import the module)
    return MF.extract_one(job, bws)


def first_visit_images(cohort_xlsx):
    coh = pd.read_excel(cohort_xlsx, "Cohort")
    coh = coh[coh["Kidney ROI (first visit, allowed size)"]]
    au = pd.read_csv(AUDIT, parse_dates=["Imaging Date"])
    a = au.merge(coh[["Patient Case ID", "Imaging Date"]].assign(**{"Imaging Date": pd.to_datetime(coh["Imaging Date"])}),
                 on=["Patient Case ID", "Imaging Date"])
    a = a[a["RF Size"].isin(ALLOWED) & a["Kidney ROI"]].copy()
    a["num"] = a["RF File"].str.extract(r"Image_(\d+)")[0].astype(int)
    a["stem"] = [f"{p}__P{int(p.split('-')[1])}_Image_{n}.pkl" for p, n in zip(a["Patient Case ID"], a.num)]
    return a


BWS_6MOD = {"SI": [9], "HScan": [3], "ESD": [1.9], "SS": [1.9], "EAC": [30], "MBF": [2.8], "Bmode": []}
OUT_CACHE_6MOD = RES / "multimap_cache_extra_6mod"

if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--six", action="store_true", help="all 6 modalities (adds EAC 30, MBF 2.8) -> multimap_cache_extra_6mod/")
    args = ap.parse_args()
    if args.six:
        BWS, OUT_CACHE = BWS_6MOD, OUT_CACHE_6MOD
    a = first_visit_images(COHORT)
    shared = set(os.listdir(SHARED))
    need = a[~a.stem.isin(shared)]
    j = pd.read_csv(JOBS, keep_default_na=False)
    j["stem"] = j.png.str.replace("|", "__", regex=False) + ".pkl"
    j = j[j.stem.isin(set(need.stem))]
    missing_jobs = sorted(set(need.stem) - set(j.stem))
    print(f"images needed {len(need)} ({need['Patient Case ID'].nunique()} patients); jobs found {len(j)}; "
          f"no job row: {missing_jobs}", flush=True)
    OUT_CACHE.mkdir(parents=True, exist_ok=True)
    L = MF.run_pool(extract_to, j.drop(columns="stem").to_dict("records"), args.workers, (BWS, str(OUT_CACHE)))
    print(L.status.str[:60].value_counts().to_string())

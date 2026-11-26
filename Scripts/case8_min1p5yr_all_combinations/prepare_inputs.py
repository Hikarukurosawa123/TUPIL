#!/usr/bin/env python3
"""Patient-level radiomics + clinical table for the all-combinations grid, >= 1.5 yr expanded-RF Case 8 cohort.

Same ROI / aggregation / feature filtering as Results_Rough/egfr_slope_all_combinations/prepare_inputs.py, applied to:
- Cohort : Results_Rough/case8_min1p5yr_expanded_rf/case8_min1p5yr_expanded_rf_slope_groups.xlsx (sheet Cohort),
           patients with a first-visit kidney ROI on an image of an allowed RF size (n = 130). Slope = that cohort's LMM.
- Images : all first-visit kidney-ROI images of an allowed RF size (images_audit.csv), per patient mean.
- ROI    : whole kidney minus renal pelvis.
- Features: PyRadiomics (multimap_filters.extract_one) SI 9, H-scan 3, MBF 2.8, ESD 1.9, EAC 30, SS 1.9; Original,
           LoG sigma 1/2/3 mm, Wavelet LL/LH/HL/HH; first-order + GLCM/GLRLM/GLSZM/GLDM/NGTDM; kept if >= 90 %
           non-missing and non-constant. Per-image pickles: shared cache (read only), else
           case8_min1p5yr_expanded_rf/multimap_cache_extra_6mod (same extract_one; verified identical on cached images).
- Clinical: baseline eGFR (CKD-EPI 2021), ACR (mg/mmol), KFRE 2-yr (%) from
           case8_min1p5yr_expanded_rf/baseline_egfr_acr_kfre_by_slope.xlsx (sheet Per patient, >= 1.5 yr).
Output: Results_Rough/case8_min1p5yr_all_combinations/inputs/patient_features_6mod.csv (+ input_summary.json).
"""
import json
import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "Scripts" / "case8_min1p5yr_expanded_rf"))
from extract_multimap_missing import first_visit_images, SHARED, OUT_CACHE_6MOD, COHORT  # noqa: E402

RES0 = ROOT / "Results_Rough" / "case8_min1p5yr_expanded_rf"
OUT = ROOT / "Results_Rough" / "case8_min1p5yr_all_combinations" / "inputs"
BIN_WIDTHS = {"SI": "9", "HScan": "3", "MBF": "2.8", "ESD": "1.9", "EAC": "30", "SS": "1.9"}
FILTERS = ("original", "log-sigma-1-0-mm-3D", "log-sigma-2-0-mm-3D", "log-sigma-3-0-mm-3D",
           "wavelet-LL", "wavelet-LH", "wavelet-HL", "wavelet-HH")
CLASSES = ("firstorder", "glcm", "glrlm", "glszm", "gldm", "ngtdm")
S = "LMM eGFR Slope (mL/min/1.73m2/year)"


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    coh = pd.read_excel(COHORT, "Cohort")
    coh = coh[coh["Kidney ROI (first visit, allowed size)"]].set_index("Patient Case ID")
    imgs = first_visit_images(COHORT)
    wanted = tuple(f"{m}|{b}|" for m, b in BIN_WIDTHS.items())
    recs = []
    for pid, stem in zip(imgs["Patient Case ID"], imgs.stem):
        # shared cache first (as the reference); images not in it come from the 6-modality extra cache.
        # 3 shared images (REB16236-134 #2, -138 #2, -160 #4) have no EAC/MBF features (no QUS map for them); those
        # values are missing for that image and the patient mean uses the patient's other images.
        shared = SHARED / stem
        src = "shared" if shared.exists() else "extra_6mod"
        r_sh = pickle.load(open(shared if shared.exists() else OUT_CACHE_6MOD / stem, "rb"))
        recs.append({"patient_id": pid, "src": src, **{k: v for k, v in r_sh.items() if k.startswith(wanted)}})
    I = pd.DataFrame(recs)
    X = I.drop(columns="src").groupby("patient_id").mean(numeric_only=True).reindex(coh.index)
    keep, ren = [], {}
    for c in X.columns:
        mod, bw, name = c.split("|", 2)
        flt, cls, _ = name.split("_", 2)
        if BIN_WIDTHS.get(mod) != bw or flt not in FILTERS or cls not in CLASSES:
            continue
        v = X[c]
        if v.notna().sum() < 0.9 * len(v) or v.std() == 0:
            continue
        keep.append(c); ren[c] = f"{mod}__{name}"
    F = X[keep].rename(columns=ren)
    clin = pd.read_excel(RES0 / "baseline_egfr_acr_kfre_by_slope.xlsx", "Per patient")
    clin = clin[clin.cohort == "≥ 1.5 yr"].set_index("Patient Case ID").reindex(F.index)
    F.insert(0, "kfre2_pct", clin.kfre2.to_numpy(float))
    F.insert(0, "acr_mg_mmol", clin.acr.to_numpy(float))
    F.insert(0, "egfr", clin.egfr.to_numpy(float))
    F.insert(0, "native_transplant", coh["Native/Transplant"].reindex(F.index).to_numpy())
    F.insert(0, "slope", coh[S].reindex(F.index).to_numpy(float))
    F.index.name = "patient_id"
    F.reset_index().to_csv(OUT / "patient_features_6mod.csv", index=False)
    slope = F.slope.to_numpy()
    per_mod = {m: int(sum(c.startswith(f"{m}__") for c in F.columns)) for m in BIN_WIDTHS}
    summary = {"patients": int(len(F)), "images": int(len(I)), "images_from_extra_6mod_cache": int((I.src == "extra_6mod").sum()),
               "features_per_modality": per_mod, "features_total": int(sum(per_mod.values())),
               "groups_4cat (>=-1 | -3..-1 | -5..-3 | <-5)": np.bincount(sum((slope < c).astype(int) for c in (-5, -3, -1))).tolist(),
               "groups_3cat (>=-1 | -3..-1 | <-3)": np.bincount(sum((slope < c).astype(int) for c in (-3, -1))).tolist(),
               "clinical_missing": {k: int(F[k].isna().sum()) for k in ("egfr", "acr_mg_mmol", "kfre2_pct")}}
    (OUT / "input_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()

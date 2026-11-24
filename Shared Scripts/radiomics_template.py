#!/usr/bin/env python3
"""Radiomics template: parametric maps + masks -> features -> cross-validated AUC.

A trimmed-down version of TUPIL_Kidney's compute_radiomics_unified.py that keeps
its method and drops the project-specific plumbing:

  1. extract   PyRadiomics (2D, 'Original' image type) on every image/mask
               pair, with a fixed bin width per modality.
  2. aggregate Average image-level features to one row per visit.
  3. model     Repeated stratified k-fold, grouped by patient. Inside every
               training fold: rank features by univariate AUC, drop near
               duplicates, keep the top N, then fit LR / RF / SVM.
               Reports the mean (and SD) of the per-fold test AUCs.

Inputs
------
manifest CSV (one row per image):
    patient_id, visit_id, modality, image_path, mask_path
    optional: image_key, mask_key  (variable name inside a .mat/.npz file)
    modality names must not contain '_' (they become column prefixes).
labels CSV (one row per visit):
    patient_id, visit_id, label     (binary 0/1)

Images/masks can be .npy, .npz, .mat, or anything SimpleITK reads (.nii, .png...).
Pixels outside the mask or non-finite are ignored.

Usage
-----
    python radiomics_template.py extract --manifest manifest.csv --out results/
    python radiomics_template.py model --features results/features_by_visit.csv \
        --labels labels.csv --out results/ --modalities SI HScan --top-n 3
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import numpy as np
import pandas as pd
import SimpleITK as sitk
from radiomics import featureextractor
from scipy.io import loadmat
from scipy.stats import spearmanr
from sklearn.base import clone
from sklearn.ensemble import RandomForestClassifier
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC

# ---------------------------------------------------------------- settings ---
# Bin widths used in the kidney study (units of each map). Modalities not
# listed fall back to a fixed number of bins.
BIN_WIDTHS = {
    "BMode": 1.5, "MBF": 1.5, "SI": 3.3, "SS": 1.0,
    "ESD": 13.0, "EAC": 10.0, "HScan": 11.0,
}
DEFAULT_BIN_COUNT = 32
PIXEL_SPACING_MM = (1.0, 1.0)   # spacing handed to PyRadiomics
MIN_MASK_PIXELS = 10
# shape2D is left out on purpose: it describes the mask, not the map, so it is
# identical across modalities and mostly measures ROI size.
FEATURE_CLASSES = ["firstorder", "glcm", "glrlm", "glszm", "gldm", "ngtdm"]

N_FOLDS = 5
N_REPEATS = 10
TOP_N = 3
MAX_ABS_RHO = 0.9               # skip a feature this correlated with one already kept
RANDOM_STATE = 42

MODELS = {
    "LogReg": LogisticRegression(C=0.1, max_iter=5000, class_weight="balanced"),
    "RandomForest": RandomForestClassifier(
        n_estimators=300, max_depth=3, min_samples_leaf=2,
        class_weight="balanced", random_state=RANDOM_STATE),
    "SVM": SVC(kernel="rbf", C=1.0, gamma="scale", class_weight="balanced"),
}

ID_COLS = ["patient_id", "visit_id"]
logging.getLogger("radiomics").setLevel(logging.ERROR)


# -------------------------------------------------------------- extraction ---
def load_array(path: str, key: str | None = None) -> np.ndarray:
    path = Path(path)
    suffix = path.suffix.lower()
    if suffix == ".npy":
        arr = np.load(path)
    elif suffix == ".npz":
        data = np.load(path)
        arr = data[key] if key else data[data.files[0]]
    elif suffix == ".mat":
        data = loadmat(path)
        if not key:
            key = next(k for k in data if not k.startswith("__"))
        arr = data[key]
    else:
        arr = sitk.GetArrayFromImage(sitk.ReadImage(str(path)))
    arr = np.squeeze(np.asarray(arr, dtype=float))
    if arr.ndim != 2:
        raise ValueError(f"{path}: expected a 2D array, got shape {arr.shape}")
    return arr


def make_extractor(modality: str) -> featureextractor.RadiomicsFeatureExtractor:
    settings = {"force2D": True}
    if modality in BIN_WIDTHS:
        settings["binWidth"] = BIN_WIDTHS[modality]
    else:
        settings["binCount"] = DEFAULT_BIN_COUNT
    extractor = featureextractor.RadiomicsFeatureExtractor(**settings)
    extractor.disableAllImageTypes()
    extractor.enableImageTypeByName("Original")
    extractor.disableAllFeatures()
    for name in FEATURE_CLASSES:
        extractor.enableFeatureClassByName(name)
    return extractor


def extract_one(extractor, image: np.ndarray, mask: np.ndarray) -> dict[str, float] | None:
    mask = (mask > 0) & np.isfinite(image)
    if mask.sum() < MIN_MASK_PIXELS:
        return None
    image = np.where(np.isfinite(image), image, 0.0)
    img = sitk.GetImageFromArray(image)
    msk = sitk.GetImageFromArray(mask.astype(np.uint8))
    img.SetSpacing(PIXEL_SPACING_MM)
    msk.SetSpacing(PIXEL_SPACING_MM)
    result = extractor.execute(img, msk, label=1)
    return {k.removeprefix("original_"): float(v)
            for k, v in result.items() if k.startswith("original_")}


def cell(row: pd.Series, col: str) -> str | None:
    value = row.get(col)
    return None if value is None or pd.isna(value) or value == "" else str(value)


def run_extract(manifest_csv: Path, out_dir: Path) -> pd.DataFrame:
    manifest = pd.read_csv(manifest_csv)
    extractors: dict[str, featureextractor.RadiomicsFeatureExtractor] = {}
    rows, skipped = [], []
    for i, r in manifest.iterrows():
        modality = str(r["modality"])
        extractor = extractors.setdefault(modality, make_extractor(modality))
        try:
            image = load_array(r["image_path"], cell(r, "image_key"))
            mask = load_array(r["mask_path"], cell(r, "mask_key"))
            if image.shape != mask.shape:
                raise ValueError(f"image {image.shape} vs mask {mask.shape}")
            features = extract_one(extractor, image, mask)
            if features is None:
                raise ValueError(f"mask has < {MIN_MASK_PIXELS} usable pixels")
        except Exception as exc:  # keep going, report at the end
            skipped.append({**r.to_dict(), "error": str(exc)})
            continue
        rows.append({**{c: r[c] for c in ID_COLS}, "modality": modality,
                     "image_path": r["image_path"], **features})
        print(f"[{i + 1}/{len(manifest)}] {modality} {r['image_path']}")

    out_dir.mkdir(parents=True, exist_ok=True)
    per_image = pd.DataFrame(rows)
    per_image.to_csv(out_dir / "features_by_image.csv", index=False)
    if skipped:
        pd.DataFrame(skipped).to_csv(out_dir / "skipped_images.csv", index=False)
        print(f"Skipped {len(skipped)} images -> {out_dir / 'skipped_images.csv'}")

    # one row per visit: mean over images, columns '<modality>_<feature>'
    feature_cols = [c for c in per_image.columns if c not in ID_COLS + ["modality", "image_path"]]
    by_visit = per_image.groupby(ID_COLS + ["modality"])[feature_cols].mean()
    by_visit = by_visit.unstack("modality")
    by_visit.columns = [f"{mod}_{feat}" for feat, mod in by_visit.columns]
    by_visit = by_visit.dropna(axis=1, how="all").reset_index()
    by_visit.to_csv(out_dir / "features_by_visit.csv", index=False)
    print(f"{len(per_image)} images -> {len(by_visit)} visits, {by_visit.shape[1] - 2} features")
    return by_visit


# ---------------------------------------------------------------- modeling ---
def univariate_auc(x: pd.Series, y: np.ndarray) -> float:
    ok = x.notna().to_numpy()
    if ok.sum() < 3 or len(np.unique(y[ok])) < 2:
        return 0.5
    return roc_auc_score(y[ok], x[ok])


def select_features(X: pd.DataFrame, y: np.ndarray, top_n: int) -> list[str]:
    """Top-N by |AUC - 0.5|, skipping features highly correlated with one already kept."""
    strength = X.apply(lambda col: abs(univariate_auc(col, y) - 0.5))
    kept: list[str] = []
    for name in strength.sort_values(ascending=False).index:
        if any(abs(spearmanr(X[name], X[k], nan_policy="omit")[0]) > MAX_ABS_RHO for k in kept):
            continue
        kept.append(name)
        if len(kept) == top_n:
            break
    return kept


def score(model, X: pd.DataFrame) -> np.ndarray:
    # decision_function for SVM: its Platt-scaled predict_proba can rank in reverse
    if hasattr(model, "decision_function"):
        return model.decision_function(X)
    return model.predict_proba(X)[:, 1]


def run_model(features_csv: Path, labels_csv: Path, out_dir: Path,
              modalities: list[str] | None, top_n: int,
              n_folds: int, n_repeats: int) -> pd.DataFrame:
    features = pd.read_csv(features_csv)
    labels = pd.read_csv(labels_csv)[ID_COLS + ["label"]]
    data = features.merge(labels, on=ID_COLS, how="inner").dropna(subset=["label"])
    feature_cols = [c for c in features.columns if c not in ID_COLS]
    if modalities:
        feature_cols = [c for c in feature_cols if c.split("_", 1)[0] in modalities]
    if not feature_cols:
        raise SystemExit("No feature columns left after the modality filter.")

    X = data[feature_cols]
    y = data["label"].astype(int).to_numpy()
    groups = data["patient_id"].to_numpy()
    print(f"{len(data)} visits ({y.sum()} positive) from {len(set(groups))} patients, "
          f"{len(feature_cols)} candidate features")

    fold_rows, picked = [], []
    for rep in range(n_repeats):
        cv = StratifiedGroupKFold(n_splits=n_folds, shuffle=True, random_state=RANDOM_STATE + rep)
        for fold, (tr, te) in enumerate(cv.split(X, y, groups)):
            if len(np.unique(y[te])) < 2:
                continue
            chosen = select_features(X.iloc[tr], y[tr], top_n)
            picked += [{"repeat": rep, "fold": fold, "feature": f} for f in chosen]
            for name, estimator in MODELS.items():
                model = make_pipeline(SimpleImputer(strategy="median"), StandardScaler(),
                                      clone(estimator))
                model.fit(X.iloc[tr][chosen], y[tr])
                auc = roc_auc_score(y[te], score(model, X.iloc[te][chosen]))
                fold_rows.append({"repeat": rep, "fold": fold, "model": name,
                                  "n_test": len(te), "auc": auc})

    folds = pd.DataFrame(fold_rows)
    summary = (folds.groupby("model")["auc"]
               .agg(fold_mean_auc="mean", fold_sd_auc="std", n_folds="count")
               .reset_index().sort_values("fold_mean_auc", ascending=False))
    selection = (pd.DataFrame(picked).groupby("feature").size()
                 .rename("times_selected").sort_values(ascending=False).reset_index())

    out_dir.mkdir(parents=True, exist_ok=True)
    folds.to_csv(out_dir / "fold_aucs.csv", index=False)
    summary.to_csv(out_dir / "auc_summary.csv", index=False)
    selection.to_csv(out_dir / "selected_features.csv", index=False)
    print(summary.to_string(index=False, float_format="%.3f"))
    print("\nMost often selected:\n" + selection.head(10).to_string(index=False))
    return summary


# --------------------------------------------------------------------- CLI ---
def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("extract", help="PyRadiomics features per image, averaged per visit")
    p.add_argument("--manifest", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)

    p = sub.add_parser("model", help="cross-validated classification on visit features")
    p.add_argument("--features", type=Path, required=True, help="features_by_visit.csv")
    p.add_argument("--labels", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--modalities", nargs="*", help="column prefixes to use, e.g. SI HScan")
    p.add_argument("--top-n", type=int, default=TOP_N)
    p.add_argument("--folds", type=int, default=N_FOLDS)
    p.add_argument("--repeats", type=int, default=N_REPEATS)

    args = parser.parse_args()
    if args.command == "extract":
        run_extract(args.manifest, args.out)
    else:
        run_model(args.features, args.labels, args.out, args.modalities,
                  args.top_n, args.folds, args.repeats)


if __name__ == "__main__":
    main()

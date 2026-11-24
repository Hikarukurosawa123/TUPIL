# Scripts

| Script | What it does |
| --- | --- |
| `qus_maps.py` | QUS parametric maps (MBF, SS, SI, Gaussian ESD/EAC) from Clarius RF data, with reference-phantom normalisation. |
| `hscan_maps.py` | H-scan maps of the kidney ROI. Self-contained: its helper modules are embedded in the file. |
| `radiomics_template.py` | Radiomics template: maps + masks → PyRadiomics features → cross-validated AUC. |
| `SMH_LF_BSC.txt` | Reference-phantom backscatter coefficient table that `qus_maps.py` reads. |

Requirements: `python -m pip install -r requirements.txt` from the repo root
(tested with Python 3.13).

## Expected data layout

Both map scripts read from `Data/Raw_data_cleaned/` in this repo by default
(`<case>/<image>/` folders holding `*_rf.raw`, `*_rf.yml`, `*_env.tgc.yml`, and an
`ROIs/` folder with the kidney ROI masks).

## QUS maps

Settings: TGC removed; STFT attenuation correction with 0.50 (kidney) / 0.80
(non-kidney) / 0.439 (phantom) dB/cm/MHz; 6×6 mm windows with 90 % overlap;
−18 dB bandwidth.

```bash
python qus_maps.py --self-test                        # numerical check, no data needed
QUS_COHORT_MODE=all_raw python qus_maps.py            # every mapped RF image
```

Other environment variables:

- `QUS_MAPPING_PATH`: image→phantom CSV. Default: `Data/Raw_data_cleaned/patient_rf_with_phantoms.csv`.
- `QUS_RAW_ROOT`: raw-data folder.
- `QUS_PHANTOM_ROOT`: folder holding the `<Phantom>_rf.raw/.yml` files.
- `QUS_OUTPUT_FOLDER_NAME`: name of the output folder.

Output goes to `Scripts/<output folder>/run_minus18dB_<timestamp>/`: one chunk of
50 cases per worker process, HDF5 maps, and CSV summaries.

## H-scan maps

Settings: TGC removed; attenuation 0.5 (kidney) / 0.8 (non-kidney); Gaussian
taper with σ = 1 MHz; filter bank 1–4 MHz.

```bash
python hscan_maps.py --raw-dir /path/to/Raw_data_cleaned --output-dir /path/to/out
python hscan_maps.py --limit 5        # quick test on five images
```

## Radiomics template

```bash
python radiomics_template.py extract --manifest manifest.csv --out results/
python radiomics_template.py model --features results/features_by_visit.csv \
    --labels labels.csv --out results/ --modalities SI HScan --top-n 3
```

Input files:

- `manifest.csv`: one row per image, with columns `patient_id, visit_id, modality, image_path, mask_path`. Optional `image_key`/`mask_key` columns give the variable name inside a `.mat` or `.npz` file.
- `labels.csv`: columns `patient_id, visit_id, label` (0/1).

Method:

- Features are averaged per visit.
- Cross-validation is stratified k-fold grouped by patient, repeated 10×.
- Feature selection (top N by univariate AUC, skipping features with |ρ| > 0.9) runs inside each training fold.
- The reported score is the mean of the per-fold test AUCs.

Bin widths, feature classes and models are set at the top of the file.

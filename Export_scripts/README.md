# Export_scripts: Results_Main figures

Rebuild everything:

```powershell
powershell -ExecutionPolicy Bypass -File "Export_scripts\run_export.ps1"
```

| Step | Script | Does |
|---|---|---|
| 1 | `export_results_main_data.py` | Reads the Results_Rough tables, recomputes box-plot statistics, writes tidy CSVs to `Results_Main/<section>/` |
| 2 | `make_results_main_figures.m` | Draws every figure from those CSVs: 600 dpi PNG + vector PDF |

Style: `Style/STYLE_GUIDE.md`. The helpers (`tupilStyle`, `tupilBoxScatter`, …) come from
`Shared Scripts/pocus_kidney_main_figures/common`.

## Outputs

| Folder | Figure | Source (Results_Rough) |
|---|---|---|
| `01_Rscore` | `rscore_heldout_by_ci_box`, `rscore_weighted_kappa_confusion` | `02_ci_Rscore_radiomics_top3overall_SI_HScan` |
| `02_Rscore_all_combinations` | `rscore_all_combinations_heatmap_rho_kappa`, `..._auc` | `ci_R_score_all_combinations/results/pivot_*.csv` |
| `03_eGFR_slope` | `egfr_all_combinations_heatmap_cindex_auc1`, `..._auc3_auc5` | `egfr_slope_all_combinations/results/*.xlsx` |
| | `egfr_best_model_box` (H-scan + ESD + SS top 1: radiomics, radiomics + ACR) | `egfr_slope_all_combinations/addon_acr` |
| | `egfr_clinical_markers_{ordinal4, binary_lt1, binary_lt3, binary_lt5}` | `01_eGFR_slope_model3_min1yr/.../baseline_clinical_model3_91.csv` |
| `04_ci_spatial_median_gap` | `auc_by_roi_early_late` | `Final_100/spatial_Final_100.xlsx` |
| | `cross_roi_median_gap` (whole-cohort CV), `cross_time_median_gap` (train early ↔ test late) | `Final_100/cross_analysis_median_gap/spatial_cross_Final_100_median.xlsx` |
| `05_ci_spatial_time_windows` | `auc_by_roi_interval`, `auc_by_roi_cumulative` | `Final_100/ci_spatial_0_1_3_5_months/spatial_time_windows_Final_100.xlsx` |
| | `cross_roi_interval`, `cross_roi_cumulative` | `Final_100/ci_spatial_0_1_3_5_months/cross_analysis/spatial_cross_Final_100_windows.xlsx` |

## Statistics

- **Box plots:** Welch ANOVA + Tukey-Kramer for 3 or more groups; Welch t-test for 2 groups. A bracket is drawn only when p < 0.05.
- **ACR and KFRE:** tested on a log10 scale. ACR is floored at 0.1 mg/mmol, the same as the model's ln-ACR term.
- **Spatial figures (04, 05):** the bar plots come in two versions, ordinal logistic (`auc_by_roi_*`) and frozen model (`auc_by_roi_*_frozen`). The cross-ROI and cross-time heatmaps are ordinal logistic only; the workbooks have no frozen-model cross results, and binary logistic stays in the CSVs.
- **Spatial AUC brackets:** one-sided, unadjusted p-values. DeLong tests from the source workbooks are paired within a window and unpaired across disjoint windows.
  - The early/late, interval and cumulative figures show brackets only for p < 0.05, including All images time comparisons.
  - Pairs per model: 6 for early/late (04, from `four_group_onesided.xlsx`) and 28 for the intervals (05, from the "ROI by time DeLong" sheet).
  - Cumulative All images comparisons use 20,000 patient bootstrap resamples stratified by outcome and entry interval. Shared patients carry the same resampling weights across windows. The tests use saved scores and do not capture model refitting variability. The source omits cumulative cross-window tests for matched and opposite ROIs.
  - The tests are one-sided (H1: first > second). A difference in the opposite direction therefore shows as ns.
  - Bracket p-values are in `*_subroi_pairs.csv`; the older ROI and window tests are in `*_tests.csv`.
- **Combined median-gap bars:** `cross_roi_time_median_gap_pairwise_tests.csv` reports two-sided DeLong tests on the saved ordinal logistic scores. Within each fitting design, tests change one training or testing ROI at a fixed time, or compare early and late test cases at fixed ROIs. For cross-time, the latter reverses the transfer direction. ROI tests are paired and time tests use the disjoint test sets. Holm correction covers all 90 contrasts. The tests use saved predictions and do not account for variability from refitting the models.
- **Heat maps:** the black box marks the best cell in each panel.

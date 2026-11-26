"""Write a compact, reproducible summary of the clinical model results."""
import argparse
import json
from pathlib import Path
from textwrap import dedent

import pandas as pd

from fit import THRESHOLDS


def table(headers, rows):
    return "\n".join(["| " + " | ".join(headers) + " |",
                     "| " + " | ".join(["---"] * len(headers)) + " |",
                     *["| " + " | ".join(map(str, row)) + " |" for row in rows]])


def report(out):
    def read(name):
        return pd.read_csv(out / f"{name}.csv")

    metrics = read("prediction_metrics")
    negative = int(read("out_of_fold_predictions_first_repeat").M3.lt(0).sum())
    pooled = metrics.query("fold == 0 and kidney == 'all' and weighting == 'patient'")
    first = pooled[pooled.repetition.eq(1)]
    gfr = first[first.outcome.eq("egfr")].set_index("model")
    change = first[first.outcome.eq("change")].set_index("model")
    repeated = pooled[pooled.outcome.eq("egfr")].groupby("model")
    intervals = read("bootstrap_ci95_first_repeat")
    m3_interval = intervals.query("scope == 'overall' and outcome == 'egfr' and weighting == 'patient' and model == 'M3' and metric == 'rmse'").iloc[0]
    slope_source = read("slope_metrics").query("repetition == 1 and fold == 0 and kidney == 'all' and reference == 'reference_lmm'")
    s = slope_source.query("minimum_reference_followup == 2").set_index("model")
    decline_source = read("rapid_decline_metrics").query(
        "repetition == 1 and fold == 0 and minimum_reference_followup == 2 and reference == 'reference_lmm' and kidney == 'all'")
    d = decline_source[decline_source.model.eq("M3")].sort_values("cutoff", ascending=False)
    lab = read("lab_filter_counts")
    cohort_summary = read("cohort_summary")
    counts_by_kind = cohort_summary[cohort_summary.minimum_training_followup.eq(0)].set_index("kidney").patients
    cohorts = cohort_summary.query("kidney == 'all'").set_index("minimum_training_followup")
    base = cohorts.loc[0]
    subgroup = metrics.query("repetition == 1 and fold == 0 and weighting == 'patient' and outcome == 'egfr' and kidney != 'all'")
    subset = subgroup.pivot(index="kidney", columns="model", values="rmse")
    horizons = read("horizon_metrics").query("repetition == 1 and weighting == 'patient' and outcome == 'egfr'")
    horizon = horizons[horizons.model.eq("M3")]
    sensitivity = pd.read_csv(out / "all_distinct_dates_sensitivity" / "metrics.csv")
    all_dates = sensitivity.query("evaluation_targets == 'all_distinct_dates' and outcome == 'egfr' and weighting == 'patient'").set_index("model")
    m3 = gfr.loc["M3"]
    m3_slope = s.loc["M3"]
    fixed = read("fixed_effects").query("model == 'M3_full_cohort'").set_index("term")
    ac = fixed.loc["t:log_acr"]

    comparison = []
    for model in ("persistence", "M0", "M1", "M2", "M3"):
        ref = repeated.get_group(model).rmse
        comparison.append((model, f"{gfr.loc[model, 'r2']:.3f}", f"{gfr.loc[model, 'mae']:.2f}",
                           f"{gfr.loc[model, 'rmse']:.2f}", f"{change.loc[model, 'r2']:.3f}",
                           f"{ref.mean():.2f} ({ref.std():.2f})"))

    decline = []
    for r in d.itertuples():
        key = f"reference_lmm_min_2y_lt_{abs(r.cutoff):g}"
        ci = intervals[(intervals.scope == "decline") & (intervals.model == "M3")
                       & (intervals.metric == "roc_auc") & (intervals.outcome == key)].iloc[0]
        decline.append((f"< {r.cutoff:g}", r.positives, r.negatives,
                        f"{r.roc_auc:.3f} ({ci.lower:.3f}-{ci.upper:.3f})", f"{r.pr_auc_ap:.3f}",
                        f"{r.sensitivity:.1%}", f"{r.specificity:.1%}"))

    selection = []
    for t in THRESHOLDS:
        label = "M3" if t == 0 else f"M3_min_{t:g}y"
        r = cohorts.loc[t]
        selection.append((f"{t:g}", int(r.patients), int(r.targets), f"{r.followup_median:.2f}",
                          f"{gfr.loc[label, 'rmse']:.2f}", f"{s.loc[label, 'rmse']:.2f}"))

    ablations = []
    for field in ("baseline_egfr", "age", "female", "log_acr"):
        label = f"without_{field}"
        a = decline_source[decline_source.model.eq(label)].set_index("cutoff")
        b = d.set_index("cutoff")
        ablations.append((field, f"{gfr.loc[label, 'rmse'] - m3.rmse:+.2f}",
                          f"{s.loc[label, 'rmse'] - m3_slope.rmse:+.2f}",
                          f"{a.loc[-3, 'roc_auc'] - b.loc[-3, 'roc_auc']:+.3f}",
                          f"{a.loc[-5, 'roc_auc'] - b.loc[-5, 'roc_auc']:+.3f}"))

    clean_table = lab[(lab.kind == "transplant") & lab.disposition.str.startswith("Before three")]
    clean_counts = clean_table.set_index("lab")
    vg = json.loads((out / "verification.json").read_text())
    episode_conflicts = read("transplant_episode_conflicts").groupby("lab").rows.sum()
    source_names = (
        "PROSPECTIVEBIOPSYFIB-ExportTransplantHist_DATA_LABELS_2026-11-13_1009.csv",
        "PROSPBX_DataRequest_10Oct2026.xlsx",
        "PROSPECTIVEBIOPSYFIB-DataEntryNativeKidne_DATA_LABELS_2026-10-22_1519.csv",
        "PROSPBX-FIBRO-US_DemographicBaseline_13Oct2026.csv",
        "PROSPBX-FIBRO-US_VisitReview_13Oct2026.csv",
    )
    inputs = "\n".join(f"- [{name}](../../All_csv/{name})" for name in source_names)
    text = f"""# Clinical longitudinal eGFR prediction

## Results

The clinical analysis included **{int(base.patients):,} patients**, with **{counts_by_kind["transplant"]} transplant and {counts_by_kind["native"]} native** trajectories, and **{int(base.targets):,} post-baseline eGFR measurements**. Median follow-up was {base.followup_median:.2f} years (IQR {base.followup_p25:.2f}-{base.followup_p75:.2f}). Median measurements per patient, including baseline, were {base.n_targets_median + 1:g} (IQR {base.n_targets_p25 + 1:g}-{base.n_targets_p75 + 1:g}). Each patient contributed one clinical baseline and one kidney episode.

M3 reached patient-weighted RMSE **{m3.rmse:.2f}** (patient bootstrap 95% CI {m3_interval.lower:.2f}-{m3_interval.upper:.2f}) and R² **{m3.r2:.3f}** for absolute future eGFR. The prediction for change from baseline had R² **{change.loc['M3', 'r2']:.3f}**. Absolute and change errors are identical when measured on the same rows because both observed and predicted eGFR have the same baseline subtracted.

{table(["Model", "Absolute R²", "MAE", "RMSE", "Change R²", "RMSE mean (SD), 10 repetitions"], comparison)}

M0 uses time. M1 adds baseline eGFR. M2 adds age and sex. M3 adds ACR. All clinical predictors interact with time. Persistence carries baseline eGFR forward. The table gives pooled out-of-fold metrics in repetition 1, with equal total weight per patient. MAE and RMSE are in mL/min/1.73 m². Observation-weighted metrics and all folds appear in [prediction_metrics.csv](prediction_metrics.csv).

### Prediction by time and kidney type

{table(["Years after baseline", "Patients", "M3 RMSE", "Persistence RMSE"], [(f"{r.horizon:g}", r.patients, f"{r.rmse:.2f}", f"{horizons[(horizons.horizon == r.horizon) & (horizons.model == 'persistence')].iloc[0].rmse:.2f}") for r in horizon.itertuples()])}

The horizon comparison uses the nearest eligible measurement within ±90 days of each target date and the fixed-effect prediction evaluated at the target date. The complete trajectory table includes extrapolations for horizons without eligible observations. Patient-weighted M3 RMSE was {subset.loc['native', 'M3']:.2f} in native patients and {subset.loc['transplant', 'M3']:.2f} in transplant recipients. The corresponding persistence RMSEs were {subset.loc['native', 'persistence']:.2f} and {subset.loc['transplant', 'persistence']:.2f}.

### Slope and rapid decline

Among {int(m3_slope.patients)} patients with at least two years of follow-up and three eGFR dates spaced at least 30 days apart, M3 annual slope RMSE was {m3_slope.rmse:.2f} mL/min/1.73 m²/year, R² was {m3_slope.r2:.3f}, and Spearman correlation was {m3_slope.spearman:.3f}. Predicted slope IQR was {m3_slope.predicted_p25:.2f} to {m3_slope.predicted_p75:.2f}. Reference slope IQR was {m3_slope.actual_p25:.2f} to {m3_slope.actual_p75:.2f}. Bland–Altman limits were {m3_slope.loa_lower:.2f} to {m3_slope.loa_upper:.2f}.

{table(["Reference slope", "Positive", "Negative", "ROC AUC (95% CI)", "Average precision", "Sensitivity", "Specificity"], decline)}

Sensitivity and specificity apply the same numerical cutoff to the predicted annual slope. Reference slopes use a time-only, random-intercept and random-slope model fitted on training patients with at least one year of follow-up and three dates 30 days apart, adapted from Case 8. For held-out reference labels, the fold's fitted population and variance parameters are held fixed while that patient's observed eGFR estimates the random slope. The prediction model uses clinical baseline predictors and fixed effects for each held-out patient. The [slope metrics](slope_metrics.csv) and [rapid decline metrics](rapid_decline_metrics.csv) include 1-, 1.5- and 2-year reference cohorts and patient-specific OLS slopes as a sensitivity analysis.

### Minimum training follow-up

{table(["Minimum years", "Available training patients before CV", "Targets before CV", "Median years", "eGFR RMSE", "Slope RMSE"], selection)}

Each minimum applies only to training patients. The held-out patient assignments and observations are identical across the five M3 fits in repetition 1. Slope RMSE in this table uses the same ≥2-year reference cohort.

The all-distinct-date sensitivity restored **{int(all_dates.iloc[0].observations - base.targets):,}** eligible follow-up dates removed by the Case 8 consecutive-day rule and held the clinical baselines and folds fixed. On the full distinct-date set, M3 trained with all dates reached patient-weighted RMSE **{all_dates.loc['M3_all_dates', 'rmse']:.2f}**, the Case 8 fit reached **{all_dates.loc['M3_case8', 'rmse']:.2f}**, and persistence reached **{all_dates.loc['persistence', 'rmse']:.2f}**. [Sensitivity metrics](all_distinct_dates_sensitivity/metrics.csv)

### Clinical predictors

Removing each predictor changed first-repetition performance by the following amounts relative to M3.

{table(["Predictor removed", "Δ eGFR RMSE", "Δ slope RMSE", "Δ ROC AUC < −3", "Δ ROC AUC < −5"], ablations)}

In the full-cohort fit, the time interaction was {ac.estimate:.2f} mL/min/1.73 m²/year per SD higher log(1 + ACR), with Wald 95% CI {ac.ci_lower:.2f} to {ac.ci_upper:.2f} and p = {ac.p_value:.3g}. See [fixed coefficients](fixed_effects.csv) for main effects and time interactions. These are association estimates from the fitted model.

## Clinical episode and baseline rules

- The transplant lab workbook provides an explicit transplant date per patient. The transplant history had to contain that exact patient-date combination. Retained labs precede any subsequent recorded transplant and recorded graft failure and fall within the workbook observation window. The raw export contained {episode_conflicts["egfr"]} eGFR and {episode_conflicts["acr"]} ACR rows whose dates reached a subsequent recorded transplant despite the older workbook label. These rows were excluded. [Episode conflict audit](transplant_episode_conflicts.csv)
- Native labs precede the first recorded transplant. Native eGFR reflects aggregate native kidney function. The native labs do not identify left or right kidney function.
- Both eGFR and ACR must be on or after three calendar months following transplantation. After earlier exclusion checks, this washout removed {int(clean_counts.loc['egfr', 'rows']):,} eGFR rows from {int(clean_counts.loc['egfr', 'patients'])} patients and {int(clean_counts.loc['acr', 'rows']):,} ACR rows from {int(clean_counts.loc['acr', 'patients'])} patients. Counts are sequential and detailed in [lab_filter_counts.csv](lab_filter_counts.csv).
- The first clinical baseline is the earliest lab day with both a same-episode eGFR and ACR available on or before that day, each within 90 days. Both baseline labs must pass the same kidney and censoring rules as follow-up. Follow-up starts the next day. Eligibility requires two distinct follow-up eGFR dates. Where age was available as visit age, age at the clinical index was approximated by shifting the recorded age by elapsed calendar time. Visit timing was used only for this calculation. The age source provides rounded ages and can postdate the lab index.
- For multiple same-day observations, the earliest timed measurement was used, with file order breaking ties or resolving entirely untimed dates. The primary analysis inherited the Case 8 consecutive-day proxy for intensive inpatient sampling. Runs of at least three consecutive eGFR dates retained only the first day. The source lacks encounter-level inpatient labels.
- Recorded RRT dates came from the sibling Case 8 workbook `case8_no_min_followup_RRT_adjusted_2p5MHz_only_egfr_with_RRT_timing.xlsx`, using its Patient_Summary and Excluded_Patients sheets. Known initiation or graft failure censored measurements on and after the event. Dated RRT prior to a transplant did not censor that graft. Patients identified as receiving RRT with an unavailable initiation date were excluded conservatively. RRT ascertainment is incomplete because the timing workbook covers a subset of patients.
- Transplant ACR used the calculated workbook field with reported ACR as a fallback when it was missing. Native ACR used the reported value. ACR was transformed as log(1 + ACR) so measured zeros were retained. eGFR used the supplied CKD-EPI 2021 field. Continuous predictors were standardized using one row per patient in each training fold. Sex used female = 1 and male = 0.
- The independent verification reconciled {vg['verified_same_episode_baselines']:,} baselines and {vg['verified_same_episode_targets']:,} follow-up measurements to source rows and reproduced {vg['verified_all_repeats_patient_date_predictions']:,} patient-date predictions from the saved fixed effects. [Verification](verification.json)

## Model files and reproduction

Five-fold patient cross-validation used kidney-type stratification. Core comparators ran in 10 independent repetitions. Minimum follow-up and leave-one-predictor-out comparisons ran in repetition 1. Intervals use 1,000 patient-cluster resamples of the first repetition's held-out results. The bootstrap reflects sampling of evaluated patients conditional on the trained fold models. Prediction model fits used maximum likelihood with correlated patient intercepts and slopes. [Mixed-effects model reference](https://www.jstatsoft.org/article/view/v067i01)

- [Fold assignments for future QUS integration](fold_assignments.csv)
- [Held-out predictions](out_of_fold_predictions_first_repeat.csv), [standardized trajectories](standardized_trajectories_first_repeat.csv), and [predicted and reference slopes](predicted_and_reference_slopes_first_repeat.csv)
- [Full-cohort M3 parameters](full_cohort_M3.json), [fixed effects](fixed_effects.csv), and [variance components](random_effects_and_residual_variance.csv)
- [Patient disposition](patient_disposition.csv), [cohort summaries](cohort_summary.csv), and [patient-bootstrap confidence intervals](bootstrap_ci95_first_repeat.csv)
- [Observed versus predicted](figures/observed_vs_predicted.png), [slope agreement](figures/slope_agreement.png), [trajectory examples](figures/example_trajectories.png), [horizon errors](figures/prediction_by_horizon.png), [rapid decline ROC and PR curves](figures/rapid_decline_discrimination.png), [fixed effects](figures/clinical_fixed_effects.png), and [comparisons](figures/model_cohort_and_kidney_comparisons.png)

Every repetition and fold has a `models.json` with preprocessing, fixed coefficients, covariance and residual variance, plus held-out predictions. For QUS augmentation, match the clinical training model and scaling to the intended outer held-out fold. These are internal cross-validation results over enrolled study patients. The longitudinal target reflects patients and dates with available follow-up, and observation and censoring patterns may affect performance. The linear model produced {negative} negative eGFR predictions at observed timestamps in the first repetition. Those predictions were retained in the scoring.

From the repository root, run the following commands with a fresh output directory. The analyzed Python environment used pandas 3.1.1, NumPy 2.5.4, SciPy 1.18.1, statsmodels 0.15.0, scikit-learn 1.9.1, Matplotlib 3.11.2 and openpyxl 3.1.5.

```powershell
$run = 'Results_Rough/egfr_longtitudinal_all_new_run'
python Scripts/egfr_longtitudinal_all/prepare.py $run
python Scripts/egfr_longtitudinal_all/fit.py $run --repetitions 10
python Scripts/egfr_longtitudinal_all/all_dates.py $run
python Scripts/egfr_longtitudinal_all/evaluate.py $run
python Scripts/egfr_longtitudinal_all/verify.py $run
python Scripts/egfr_longtitudinal_all/figures.py $run
python Scripts/egfr_longtitudinal_all/report.py $run
```

### Source exports

{inputs}
"""
    (out / "README.md").write_text(dedent(text), encoding="utf-8")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    report(parser.parse_args().output)



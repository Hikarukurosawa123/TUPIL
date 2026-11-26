"""Score paired held-out predictions and summarize the nested QUS experiment."""
from collections import Counter
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from model import OUT, GRID, MODELS, name
from evaluate import metrics, classifiers

GRID_NAMES = tuple(map(name, GRID))


def score(out, draws=2000):
    cohort = pd.read_csv(out / "scan_cohort.csv", index_col="pid")
    assignments = pd.read_csv(out / "fold_assignments.csv")
    postscan = pd.read_csv(out / "postscan_outcomes.csv")
    variants = MODELS + GRID_NAMES
    entries, frames, slopes, select = [], [], [], []
    for (repeat, fold), group in assignments.groupby(["repetition", "fold"]):
        path = out / "cv" / f"r{repeat:02d}_f{fold}"
        obs = pd.read_csv(path / "oof.csv")
        ref = pd.read_csv(path / "slopes.csv")
        params = json.loads((path / "models.json").read_text())
        ids = set(group.pid)
        assert ids == set(obs.pid) == set(ref.pid)
        assert not obs.duplicated(["pid", "date"]).any()
        assert len(obs) == len(postscan[postscan.pid.isin(ids)])
        assert len(cohort) - len(ids) == params["clinical"]["training_patients"]
        assert len(cohort) - len(ids) == params["reference"]["training_patients"]
        assert params["clinical"]["converged"] and params["reference"]["converged"]
        for v in variants[2:]:
            item = params["qus"][v]
            assert len(item["features"]) == len(set(item["features"]))
            delta = (obs[v] - obs.groupby("pid")[v].transform("first")
                     - (obs.t - obs.groupby("pid").t.transform("first")) * obs.pid.map(ref.set_index("pid")[v]))
            assert delta.abs().max() < 1e-8
            for feature, rho in zip(item["features"], item["rho"]):
                select.append(dict(repetition=repeat, fold=fold, model=v, chosen_grid=item["config"],
                                   penalty=item["penalty"], feature=feature, train_spearman=rho,
                                   intercept_beta=item["qus_intercept_beta"][feature], slope_beta=item["qus_slope_beta"][feature]))
        frames.append(obs.assign(repetition=repeat, fold=fold))
        slopes.append(ref.assign(repetition=repeat, fold=fold))

    raw = pd.concat(frames, ignore_index=True)
    refs = pd.concat(slopes, ignore_index=True)
    assert not raw.duplicated(["repetition", "pid", "date"]).any()
    assert not refs.duplicated(["repetition", "pid"]).any()
    raw.to_csv(out / "oof_predictions.csv.gz", index=False, compression="gzip")
    refs.to_csv(out / "oof_slopes.csv.gz", index=False, compression="gzip")
    select = pd.DataFrame(select)
    select.to_csv(out / "selected_features_by_fold.csv", index=False)
    (select.groupby(["model", "feature"]).agg(
        folds=("fold", "size"), mean_absolute_train_rho=("train_spearman", lambda v: v.abs().mean()),
        mean_slope_beta=("slope_beta", "mean"), mean_intercept_beta=("intercept_beta", "mean"))
     .reset_index().sort_values(["model", "folds"], ascending=[True, False])
     .to_csv(out / "feature_frequency.csv", index=False))
    (select[["model", "repetition", "fold", "chosen_grid", "penalty"]].drop_duplicates()
     .groupby(["model", "chosen_grid", "penalty"]).size().rename("folds").reset_index()
     .to_csv(out / "configuration_frequency.csv", index=False))

    for repeat, r in raw.groupby("repetition"):
        label = refs[refs.repetition.eq(repeat)].set_index("pid")
        assert len(r) == len(postscan) and len(label) == len(cohort)
        near = r.assign(distance=(r.t - 2).abs())
        near = near[near.distance.le(90 / 365.25)].sort_values(["distance", "date"]).drop_duplicates("pid")
        for model in variants:
            for scope, d in [("all_postscan", r), ("through_2_years", r[r.t.le(2)]), ("nearest_2_years", near)]:
                weight = 1 / d.groupby("pid").pid.transform("size")
                valid = dict(repetition=repeat, model=model, patients=d.pid.nunique(), observations=len(d), target="egfr", window=scope)
                entries.append(valid | {"weight": "patient"} | metrics(d.egfr, d[model], weight))
                entries.append(valid | {"weight": "observation"} | metrics(d.egfr, d[model]))
                if scope == "all_postscan":
                    baseline = d.pid.map(cohort.baseline_egfr)
                    entries.append(valid | {"weight": "patient", "target": "change"} |
                                   metrics(d.egfr - baseline, d[model] - baseline, weight))
            for ref_name in ("reference_lmm", "reference_ols"):
                for ref_target in ["slope", "rapid_decline_-3", "rapid_decline_-5"]:
                    valid = dict(repetition=repeat, model=model, patients=len(label), observations=len(label),
                                 target=ref_target, window=ref_name, weight="patient")
                    results = (metrics(label[ref_name], label[model]) if ref_target == "slope" else
                               classifiers(label[ref_name], label[model], -int(ref_target[-1])))
                    entries.append(valid | results)

    perf = pd.DataFrame(entries)
    perf.to_csv(out / "metrics_by_repetition.csv", index=False)
    measure_cols = perf.select_dtypes(include="number").columns.difference(["repetition"])
    agg = perf.groupby(["target", "window", "weight", "model"])[list(measure_cols)].agg(["mean", "std"])
    agg.columns = [f"{metric}_{stat}" for metric, stat in agg.columns]
    agg.reset_index().to_csv(out / "mean_metrics.csv", index=False)
    contrasts = paired_bootstrap(raw, refs, variants, draws)
    contrasts.to_csv(out / "paired_rmse_differences.csv", index=False)
    report(out, agg.reset_index(), contrasts, select)
    figures(out, agg.reset_index())


def paired_bootstrap(raw, refs, variants, draws):
    rng = np.random.default_rng(20261124)
    result = []
    for key, d, actual in [
        ("egfr_all_postscan", raw, "egfr"),
        ("egfr_through_2_years", raw[raw.t.le(2)], "egfr"),
        ("slope_lmm", refs, "reference_lmm"),
        ("slope_ols", refs, "reference_ols"),
    ]:
        squared = d[list(variants)].subtract(d[actual], axis="index")**2
        values = squared.groupby(d.pid).mean()
        n = len(values)
        samples = rng.multinomial(n, np.full(n, 1 / n), size=draws)
        rmse = np.sqrt(samples @ values.to_numpy() / n)
        baseline = list(variants).index("clinical")
        point = np.sqrt(values.mean())
        for j, model in enumerate(variants):
            delta = rmse[:, j] - rmse[:, baseline]
            lo, hi = np.quantile(delta, [.025, .975])
            result.append(dict(target=key, model=model, patients=n, rmse=point[model],
                               clinical_rmse=point["clinical"], delta=point[model] - point["clinical"],
                               low=lo, high=hi, bootstrap_draws=draws))
    return pd.DataFrame(result)


def table(headers, rows):
    return "\n".join(["| " + " | ".join(headers) + " |",
                       "| " + " | ".join(["---"] * len(headers)) + " |",
                       *["| " + " | ".join(map(str, row)) + " |" for row in rows]])


def report(out, agg, ci, select):
    flow = json.loads((out / "cohort_summary.json").read_text())
    a = agg[agg.weight.eq("patient")].set_index(["target", "window", "model"])
    contrast = ci.set_index(["target", "model"])
    names = dict(persistence="Carry baseline forward", clinical="Clinical only",
                 qus_nested_egfr="Clinical + QUS (tuned for eGFR)", qus_nested_slope="Clinical + QUS (tuned for slope)")
    values = []
    for key, title in names.items():
        egfr = a.loc["egfr", "all_postscan", key]
        two = a.loc["egfr", "through_2_years", key]
        slope = a.loc["slope", "reference_lmm", key]
        d3 = a.loc["rapid_decline_-3", "reference_lmm", key]
        d5 = a.loc["rapid_decline_-5", "reference_lmm", key]
        values.append((title, f"{egfr.rmse_mean:.2f}", f"{egfr.r2_mean:.3f}", f"{two.rmse_mean:.2f}",
                       f"{slope.rmse_mean:.2f}", f"{slope.r2_mean:.3f}", f"{d3.roc_auc_mean:.3f}", f"{d5.roc_auc_mean:.3f}"))
    differences = []
    for endpoint in ("egfr_all_postscan", "egfr_through_2_years", "slope_lmm", "slope_ols"):
        for key in ("qus_nested_egfr", "qus_nested_slope"):
            d = contrast.loc[endpoint, key]
            differences.append((endpoint, names[key], f"{d.delta:+.2f}", f"{d.low:+.2f} to {d.high:+.2f}"))
    grid = []
    for target, scope, top in GRID:
        if scope not in ("overall", "per_modality"):
            continue
        key = name((target, scope, top))
        e = a.loc["egfr", "all_postscan", key]
        s = a.loc["slope", "reference_lmm", key]
        grid.append((target, scope.replace("_", " "), top, top * (6 if scope == "per_modality" else 1),
                     f"{e.rmse_mean:.2f}", f"{s.rmse_mean:.2f}"))
    best_single = []
    for scope in ("SI", "HScan", "MBF", "ESD", "EAC", "SS"):
        keys = [name(c) for c in GRID if c[1] == scope]
        best_e = min(keys, key=lambda k: a.loc["egfr", "all_postscan", k].rmse_mean)
        best_s = min(keys, key=lambda k: a.loc["slope", "reference_lmm", k].rmse_mean)
        best_single.append((scope, best_e.replace("__", " "), f"{a.loc['egfr', 'all_postscan', best_e].rmse_mean:.2f}",
                            best_s.replace("__", " "), f"{a.loc['slope', 'reference_lmm', best_s].rmse_mean:.2f}"))
    frequency = select[select.model.isin(["qus_nested_egfr", "qus_nested_slope"])]
    top_features = []
    for key, features in frequency.groupby("model"):
        freq = features.feature.value_counts().head(6)
        top_features.extend((names[key], feature, n) for feature, n in freq.items())
    transfer = pd.read_csv(out / "stage1_transfer_sensitivity" / "paired_rmse_differences.csv")
    transfer_table = []
    for row in transfer[transfer.target.isin(["egfr_all_postscan", "slope_lmm"])].itertuples():
        transfer_table.append((row.target, row.model, f"{row.rmse:.2f}",
                               f"{row.delta:+.2f} ({row.low:+.2f} to {row.high:+.2f})"))
    text = f"""# QUS-augmented longitudinal eGFR prediction

## Held-out performance

The designated Case 8 workbook supplied {flow['case8_patients']} patients with at least two years of follow-up. The analysis included **{flow['patients']}** with a qualifying QUS kidney ROI ({flow['counts_by_kind']['transplant']} transplant, {flow['counts_by_kind']['native']} native), across **{flow['images']}** first-scan images and **{flow['observations']:,}** same-episode post-scan eGFR measurements.

{table(["Model", "eGFR RMSE", "eGFR R²", "First 2 years RMSE", "Slope RMSE", "Slope R²", "AUC below −3", "AUC below −5"], values)}

Values are means over 10 repetitions of five-fold patient-level validation. eGFR measures use equal total weight per patient. Slope refers to annual eGFR slope from a fold-specific time-only mixed model. The rapid-decline cutoffs are in mL/min/1.73 m²/year. The two QUS rows use separate inner-loop selection criteria and share the same outer held-out patients as clinical-only.

### Paired differences from clinical-only

{table(["Outcome", "Model", "RMSE change", "95% patient-bootstrap interval"], differences)}

Negative values favor adding QUS. The interval resamples patients with all their observations and repeats together. Each contrast pools the squared errors across repeated held-out predictions before taking the root. This differs slightly from averaging per-repeat RMSE above. The interval reflects evaluation-sample variability conditional on the fitted models.

### Requested QUS grid

{table(["Training ranking target", "Selection", "Top N", "QUS predictors", "eGFR RMSE", "Slope RMSE"], grid)}

Each grid row tuned ridge strength on inner-validation eGFR error. The two nested results above also tune the ranking target, feature pool, and N. The full exploratory grid tested each modality on its own. The most favorable single-modality rows below were identified after comparing outer results.

{table(["Modality", "Lowest eGFR RMSE configuration", "eGFR RMSE", "Lowest slope RMSE configuration", "Slope RMSE"], best_single)}

Most frequently selected features in the nested models appear below. These counts describe selection stability across folds.

{table(["Model", "Feature", "Folds selected"], top_features)}

### Using the larger Stage 1 clinical fit

{table(["Outcome", "Model", "Pooled RMSE", "Change from recalibrated Stage 1 (95% interval)"], transfer_table)}

The sensitivity analysis applies the original Stage 1 fold's clinical coefficients to scan-time clinical covariates. The original Stage 1 training data exclude each scanned held-out patient. Both variants recalibrate the intercept and time slope using scanned training patients. The QUS variants add ridge-penalized QUS level and slope corrections using the features and penalty chosen inside the primary experiment. This tests whether starting from the larger clinical training set changes the incremental QUS result. [Transfer outputs](stage1_transfer_sensitivity/paired_rmse_differences.csv)

## Methods

- The prediction index is the qualifying scan date. Clinical-only and QUS models used the same scan-index cohorts, observations, predictors, and folds. The clinical model preserves the Stage 1 predictor and random-effects structure and was refit on the scanned training patients. The original Stage 1 clinical baseline preceded some scans by years.
- Baseline eGFR and ACR are the latest eligible same-episode measurements within 365 days on or before the scan. The clinical data are day-level. The analysis includes {flow['missing_scan_acr']} patient with missing prescan ACR, filled by the training-fold median of log(1 + ACR). Age is the Stage 1 age approximation shifted to the scan date.
- All outcomes are strictly post-scan and pass Stage 1 kidney episode, graft failure, RRT, three-month transplant, and Case 8 concentrated-sampling filters. The supplied Case 8 follow-up rule determines inclusion. In the strictly post-scan audited labs, {flow['postscan_span_lt2']} patients have less than two years between the first and last retained outcome, even though all reach at least two years after scan.
- The six modalities are SI, MBF, ESD, EAC, SS, and H-scan. First-order intensity and GLCM, GLRLM, GLSZM, GLDM, and NGTDM textures were averaged across the qualifying kidney ROI images for each patient. The cache supplies original images, LoG at 1–3 mm, and four first-level wavelet bands. Bin widths were SI 9, MBF 2.8, ESD 1.9, EAC 30, SS 1.9, and H-scan 3. The source workflow generally used the kidney mask minus the pelvis and used available whole-kidney masks for images without a pelvis mask.
- Each training split applied 90% feature completeness, median imputation, variance filtering, standardization, absolute Spearman ranking, and greedy absolute Pearson redundancy filtering at 0.90. The category ranking used four ordered strata with boundaries −1, −3, and −5. Per-modality selection takes N per modality and pools the six sets. Overall selection ranks all modalities together and takes N total. Training slope and category labels were derived entirely within each split from post-scan observations. The workbook's published slopes were retained for provenance.
- Clinical-only prediction uses a maximum-likelihood random-intercept and random-slope model. The augmented model estimates QUS main effects and time interactions jointly with the clinical coefficients by generalized least squares, holding the training clinical model's covariance fixed. It applies ridge shrinkage to the QUS coefficients. Held-out prediction uses fixed effects. Held-out reference slopes use the fold-trained population parameters with each held-out patient's measured trajectory. OLS reference slopes are a sensitivity analysis.
- Outer folds reuse Stage 1 assignments. Inner three-fold patient validation selects ridge strength and the reported nested configurations. See the [nested validation reference](https://pubmed.ncbi.nlm.nih.gov/16504092/).
- Results describe the small, mostly transplant scan cohort. Five workbook patients lacked qualifying QUS features. The analysis evaluates performance only where QUS is available.

## Results files

- [Full metrics](mean_metrics.csv), [metrics by repetition](metrics_by_repetition.csv), and [paired bootstrap contrasts](paired_rmse_differences.csv)
- [Selected features by fold](selected_features_by_fold.csv) and [feature frequencies](feature_frequency.csv)
- [Held-out eGFR predictions](oof_predictions.csv.gz), [held-out slopes](oof_slopes.csv.gz), and [fold assignments](fold_assignments.csv)
- [Predictor grid and comparison figure](qus_comparison.png) and [cohort summary](cohort_summary.json)
- Per-fold models, coefficients, scaling, inner-loop scores, and predictions in the `cv` directory

Rerun `prepare.py`, `model.py`, `transfer.py`, and `results.py` in order from `Scripts/egfr_longtitudinal_qus` with a fresh output directory. Preparation reads the original trusted radiomics cache and the saved Stage 1 files.
"""
    (out / "README.md").write_text(text, encoding="utf-8")


def figures(out, agg):
    a = agg[agg.weight.eq("patient")].set_index(["target", "window", "model"])
    plt.rcParams.update({"font.family": "Arial", "font.size": 8, "axes.spines.top": False, "axes.spines.right": False})
    blue, rose, gray = "#4C9AC6", "#E84A75", "#555555"
    fig, axes = plt.subplots(1, 3, figsize=(9, 3.35), layout="constrained")
    for ax, target, window in zip(axes[:2], ("egfr", "slope"), ("all_postscan", "reference_lmm")):
        values = [a.loc[target, window, v].rmse_mean for v in MODELS]
        ax.barh(["Persistence", "Clinical", "QUS tuned for eGFR", "QUS tuned for slope"], values,
                color=[gray, blue, rose, rose])
        ax.invert_yaxis()
        ax.set_xlabel("eGFR RMSE" if target == "egfr" else "Annual slope RMSE")
    ax = axes[2]
    reference = a.loc["egfr", "all_postscan", "clinical"].rmse_mean
    for target, linestyle in [("slope", "-"), ("category", "--")]:
        for scope, color in [("overall", blue), ("per_modality", rose)]:
            values = [a.loc["egfr", "all_postscan", f"{target}__{scope}__top{n}"].rmse_mean for n in (1, 3, 5)]
            ax.plot([1, 3, 5], np.array(values) - reference, marker="o", linestyle=linestyle,
                    color=color, label=f"{target.capitalize()}, {scope.replace('_', ' ')}")
    ax.axhline(0, ls=":", color=gray, lw=.8)
    ax.set(xlabel="Top N", ylabel="eGFR RMSE change from clinical", xticks=[1, 3, 5])
    ax.legend(frameon=False, fontsize=6.5)
    fig.savefig(out / "qus_comparison.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    score(OUT)


"""Evaluate held-out predictions, including patient-cluster confidence intervals."""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score

from fit import HORIZONS, MODELS, THRESHOLDS
from prepare import YEAR

METADATA = {"pid", "series", "date", "t", "egfr", "source", "source_row",
            "kind", "baseline_egfr", "repetition", "fold"}


def metrics(actual, predicted, weights=None):
    y, p = np.asarray(actual), np.asarray(predicted)
    w = np.ones(len(y)) if weights is None else np.asarray(weights)
    if not len(y):
        return {}
    w = w / w.sum()

    def moments(a, b):
        da, db = a - w @ a, b - w @ b
        return w @ (da * da), w @ (db * db), w @ (da * db)

    def weighted_rank(a):
        _, index = np.unique(a, return_inverse=True)
        counts = np.bincount(index, weights=w)
        return (np.cumsum(counts) - counts / 2)[index]

    def correlation(va, vb, cov):
        return cov / np.sqrt(va * vb) if va > 1e-20 and vb > 1e-20 else np.nan

    vy, vp, cov = moments(y, p)
    rank = moments(weighted_rank(y), weighted_rank(p))
    err = p - y
    mse = w @ (err * err)
    bias = w @ err
    denom = vy + vp + bias * bias
    return dict(r2=1 - mse / vy if vy > 0 else np.nan, mae=w @ np.abs(err), rmse=np.sqrt(mse),
                pearson=correlation(vy, vp, cov), spearman=correlation(*rank),
                bias=bias, ccc=2 * cov / denom if denom > 0 else np.nan)


def classifiers(actual, predicted, cutoff):
    y = np.asarray(actual) < cutoff
    score = -np.asarray(predicted)
    call = score > -cutoff
    tp, tn = sum(y & call), sum(~y & ~call)
    positives, n = sum(y), len(y)
    negatives, fp, fn = n - positives, sum(~y & call), sum(y & ~call)

    def div(a, b):
        return a / b if b else np.nan

    se, sp = div(tp, positives), div(tn, negatives)
    return dict(positives=positives, negatives=negatives,
                roc_auc=roc_auc_score(y, score) if positives and negatives else np.nan,
                pr_auc_ap=average_precision_score(y, score) if positives and negatives else np.nan,
                balanced_accuracy=(se + sp) / 2, sensitivity=se, specificity=sp,
                ppv=div(tp, tp + fp), npv=div(tn, tn + fn), f1=div(2 * tp, 2 * tp + fp + fn))


def read_cv(out, table, repetition):
    root = out / "cv" / f"rep_{repetition:02d}"
    paths = sorted(root.glob(f"fold_*/{table}.csv"))
    assert len(paths) == 5
    return pd.concat([pd.read_csv(path).assign(fold=int(path.parent.name.split("_")[1]),
                                              repetition=repetition) for path in paths], ignore_index=True)


def trajectory_metrics(obs, identity, models, options=("egfr", "change")):
    if obs.empty:
        return
    by_patient = 1 / obs.groupby("pid").pid.transform("size").to_numpy()
    for model in models:
        for outcome in options:
            baseline = 0 if outcome == "egfr" else obs.baseline_egfr
            for scheme, weight in [("observation", None), ("patient", by_patient)]:
                scores = metrics(obs.egfr - baseline, obs[model] - baseline, weight)
                yield identity | dict(model=model, outcome=outcome, weighting=scheme,
                                      patients=obs.pid.nunique(), observations=len(obs)) | scores


def slope_metrics(slopes, identity, models):
    for minimum in (1, 1.5, 2):
        d = slopes[slopes.followup.ge(minimum) & slopes.spaced_dates.ge(3)]
        if d.empty:
            continue
        for reference in ("reference_lmm", "reference_ols"):
            actual = d[reference]
            for model in models:
                pred = d[model]
                diff = pred - actual
                row = identity | dict(model=model, minimum_reference_followup=minimum,
                                       reference=reference, patients=len(d))
                yield (row | metrics(actual, pred) |
                       dict(loa_lower=diff.mean() - 1.96 * diff.std(), loa_upper=diff.mean() + 1.96 * diff.std(),
                            actual_p25=actual.quantile(.25), actual_median=actual.median(), actual_p75=actual.quantile(.75),
                            predicted_p25=pred.quantile(.25), predicted_median=pred.median(),
                            predicted_p75=pred.quantile(.75)))


def decline_metrics(slopes, identity, models):
    for minimum in (1, 1.5, 2):
        d = slopes[slopes.followup.ge(minimum) & slopes.spaced_dates.ge(3)]
        if d.empty:
            continue
        for reference in ("reference_lmm", "reference_ols"):
            for model in models:
                for cutoff in (-3, -5):
                    yield identity | dict(model=model, minimum_reference_followup=minimum, reference=reference,
                                          cutoff=cutoff, patients=len(d)) | classifiers(d[reference], d[model], cutoff)


def bootstrap_trajectories(obs, models, draws, rng):
    ys = obs.egfr.to_numpy()
    base = obs.baseline_egfr.to_numpy()
    pid_index, patients = pd.factorize(obs.pid, sort=True)
    n = len(patients)
    sizes = np.bincount(pid_index)
    resamples = rng.multinomial(n, np.full(n, 1 / n), size=draws)
    for scheme in ("patient", "observation"):
        row_weight = 1 / sizes[pid_index] if scheme == "patient" else np.ones(len(obs))

        def per_patient(array):
            return np.bincount(pid_index, weights=array * row_weight, minlength=n)

        N = resamples @ per_patient(np.ones(len(obs)))
        for target in ("egfr", "change"):
            y = ys if target == "egfr" else ys - base
            TSS = resamples @ per_patient(y**2) - (resamples @ per_patient(y))**2 / N
            for model in models:
                err = obs[model].to_numpy() - ys
                mse = resamples @ per_patient(err**2) / N
                vals = dict(mae=(resamples @ per_patient(np.abs(err))) / N,
                            rmse=np.sqrt(mse), r2=1 - mse * N / TSS,
                            bias=(resamples @ per_patient(err)) / N)
                for metric, vector in vals.items():
                    lo, hi = np.quantile(vector, [.025, .975])
                    yield dict(scope="overall", model=model, outcome=target, weighting=scheme, metric=metric,
                               lower=lo, upper=hi, bootstrap_draws=draws)


def bootstrap_decline(slopes, draws, rng):
    d = slopes[slopes.followup.ge(2) & slopes.spaced_dates.ge(3)]
    n = len(d)
    selections = rng.integers(0, n, size=(draws, n))
    actual = d.reference_lmm.to_numpy()
    for model in ("M1", "M2", "M3"):
        pred = d[model].to_numpy()
        vals = [metrics(actual[sel], pred[sel]) for sel in selections]
        for metric in vals[0]:
            lo, hi = np.nanquantile([v[metric] for v in vals], [.025, .975])
            yield dict(scope="slope", model=model, outcome="reference_lmm_min_2y",
                       weighting="patient", metric=metric, lower=lo, upper=hi, bootstrap_draws=draws)
        for cutoff in (-3, -5):
            vals = [classifiers(actual[sel], pred[sel], cutoff) for sel in selections]
            for metric in ("roc_auc", "pr_auc_ap"):
                lo, hi = np.nanquantile([v[metric] for v in vals], [.025, .975])
                yield dict(scope="decline", model=model, outcome=f"reference_lmm_min_2y_lt_{abs(cutoff)}",
                           weighting="patient", metric=metric, lower=lo, upper=hi, bootstrap_draws=draws)


def cohort_summaries(cohort):
    for minimum in THRESHOLDS:
        selection = cohort[cohort.followup.ge(minimum)]
        for label, d in [("all", selection), *selection.groupby("kind")]:
            row = dict(minimum_training_followup=minimum, kidney=label, patients=len(d),
                       targets=sum(d.n_targets), egfr_including_baseline=sum(d.n_targets) + len(d),
                       female=sum(d.female))
            for field in ("age", "baseline_egfr", "acr", "followup", "n_targets"):
                for name, q in [("min", 0), ("p25", .25), ("median", .5), ("p75", .75), ("max", 1)]:
                    row[f"{field}_{name}"] = d[field].quantile(q)
            yield row


def evaluate(out, draws):
    cohort = pd.read_csv(out / "clinical_baselines.csv", index_col="pid")
    assignments = pd.read_csv(out / "fold_assignments.csv")
    eligible = pd.read_csv(out / "eligible_followup.csv")
    regression, slope, decline, horizons, cis = [], [], [], [], []
    for repetition, fold_map in assignments.groupby("repetition"):
        assert len(fold_map) == len(cohort)
        assert not fold_map.pid.duplicated().any() and set(fold_map.pid) == set(cohort.index)
        pred = read_cv(out, "predictions", repetition).join(cohort[["baseline_egfr", "kind"]], on="pid")
        slopes = read_cv(out, "slopes", repetition).join(cohort[["followup", "kind"]], on="pid")
        trajectory = read_cv(out, "standardized_trajectories", repetition)
        assert not pred.duplicated(["pid", "date"]).any()
        assert len(pred) == len(eligible)
        assert pred.merge(fold_map[["pid", "fold"]], on=["pid", "fold"], how="inner").shape[0] == len(pred)
        assert not slopes.pid.duplicated().any()
        assert len(trajectory) == len(cohort) * len(HORIZONS)
        names = [c for c in pred if c not in METADATA]
        pred_and_slopes = pred.merge(slopes[["pid", "M3"]], on="pid", suffixes=("", "_slope"))
        assert np.allclose(pred_and_slopes["M3"] - pred_and_slopes["M3_slope"] * pred_and_slopes.t,
                           pred_and_slopes.groupby("pid").M3.transform("first")
                           - pred_and_slopes["M3_slope"] * pred_and_slopes.groupby("pid").t.transform("first"))
        tags = dict(repetition=repetition, fold=0, kidney="all")
        for fold, d in [(0, pred), *pred.groupby("fold")]:
            ident = tags | dict(fold=fold)
            regression.extend(trajectory_metrics(d, ident, names))
            s = slopes if not fold else slopes[slopes.fold.eq(fold)]
            slope.extend(slope_metrics(s, ident, names))
            decline.extend(decline_metrics(s, ident, names))
        if repetition == 1:
            for kind, d in pred.groupby("kind"):
                regression.extend(trajectory_metrics(d, tags | dict(kidney=kind), names))
                s = slopes[slopes.kind.eq(kind)]
                slope.extend(slope_metrics(s, tags | dict(kidney=kind), names))
                decline.extend(decline_metrics(s, tags | dict(kidney=kind), names))
            pred.to_csv(out / "out_of_fold_predictions_first_repeat.csv", index=False)
            slopes.to_csv(out / "predicted_and_reference_slopes_first_repeat.csv", index=False)
            trajectory.to_csv(out / "standardized_trajectories_first_repeat.csv", index=False)
            rng = np.random.default_rng(737)
            cis.extend(bootstrap_trajectories(pred, ["persistence", *MODELS], draws, rng))
            cis.extend(bootstrap_decline(slopes, draws, rng))
        for horizon in HORIZONS:
            delta = (pred.t - horizon).abs()
            candidates = pred.loc[delta.le(90 / YEAR), ["pid", "date", "t", "egfr", "baseline_egfr"]].copy()
            candidates["horizon_distance"] = delta[candidates.index]
            nearest = candidates.sort_values(["horizon_distance", "date"]).drop_duplicates("pid")
            exact = trajectory[trajectory.t.eq(horizon)].drop(columns=["t", "fold", "repetition"])
            merged = nearest.merge(exact, on="pid", validate="1:1")
            horizons.extend(trajectory_metrics(merged, tags | dict(horizon=horizon), names))
    outputs = {"prediction_metrics": regression, "slope_metrics": slope, "rapid_decline_metrics": decline,
               "horizon_metrics": horizons, "bootstrap_ci95_first_repeat": cis,
               "cohort_summary": list(cohort_summaries(cohort))}
    for name, rows in outputs.items():
        pd.DataFrame(rows).to_csv(out / f"{name}.csv", index=False)
    for name, id_cols in [
        ("prediction_metrics", ["model", "outcome", "weighting"]),
        ("slope_metrics", ["model", "reference", "minimum_reference_followup"]),
        ("rapid_decline_metrics", ["model", "reference", "minimum_reference_followup", "cutoff"]),
        ("horizon_metrics", ["model", "outcome", "weighting", "horizon"]),
    ]:
        d = pd.DataFrame(outputs[name])
        d = d[d.kidney.eq("all")]
        for label, selected in [("fold", d[d.fold.gt(0)]), ("repetition", d[d.fold.eq(0)])]:
            if len(selected):
                value_columns = list(set(selected.select_dtypes("number").columns)
                                     - {"fold", "repetition", *id_cols})
                g = selected.groupby(id_cols)[sorted(value_columns)].agg(["mean", "std", "count"])
                g.columns = [f"{metric}_{stat}" for metric, stat in g.columns]
                g.to_csv(out / f"{name}_across_{label}s.csv")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    parser.add_argument("--bootstrap", type=int, default=1000)
    args = parser.parse_args()
    evaluate(args.output, args.bootstrap)

"""Test QUS corrections on the broader Stage 1 clinical fit at scan time."""
import json

import numpy as np
import pandas as pd
from scipy.linalg import cholesky, solve, solve_triangular

from model import FEATURES, OUT
from prepare import STAGE1
from evaluate import metrics

NAMES = ("stage1_scan_calibrated", "stage1_plus_qus_egfr", "stage1_plus_qus_slope")


def base_at_scan(cohort, train, params):
    z = cohort[FEATURES].copy()
    z["log_acr"] = z.log_acr.fillna(cohort.loc[train, "log_acr"].median())
    scales = params["preprocessing"]
    for feature, mean in scales["center"].items():
        z[feature] = (z[feature] - mean) / scales["scale"][feature]
    beta = params["fixed_effects"]
    level = beta["Intercept"] + sum(beta[c] * z[c] for c in FEATURES)
    slope = beta["t"] + sum(beta[f"t:{c}"] * z[c] for c in FEATURES)
    return level, slope


def corrections(train, test, measurements, covariance, residual, level, slope, x, params):
    white_z, white_y, train_owners = [], [], []
    for i, (pid, d) in enumerate(measurements[measurements.pid.isin(train)].groupby("pid", sort=True)):
        assert pid == train[i]
        z = np.column_stack((np.ones(len(d)), d.t))
        lower = cholesky(z @ covariance @ z.T + residual * np.eye(len(d)), lower=True)
        white_z.append(solve_triangular(lower, z, lower=True))
        error = d.egfr - level[pid] - d.t * slope[pid]
        white_y.append(solve_triangular(lower, error, lower=True))
        train_owners.extend([i] * len(d))
    z, y = np.concatenate(white_z), np.concatenate(white_y)
    train_owners = np.array(train_owners)
    predictions = {}
    b = solve(z.T @ z, z.T @ y, assume_a="pos")
    predictions[NAMES[0]] = (level.loc[test] + b[0], slope.loc[test] + b[1])
    fit = {NAMES[0]: {"calibration_intercept": b[0], "calibration_slope": b[1]}}
    for key, model_name in zip(("qus_nested_egfr", "qus_nested_slope"), NAMES[1:]):
        model = params[key]
        names = model["features"]
        q = ((x[names].fillna(pd.Series(model["feature_median"]))
              - pd.Series(model["feature_center"])) / pd.Series(model["feature_scale"]))
        q_train = q.loc[train].to_numpy()[train_owners]
        design = np.column_stack((z, q_train * z[:, :1], q_train * z[:, 1:]))
        n = len(names)
        penalty = np.r_[0, 0, np.repeat(1 / np.maximum(np.diag(covariance), 1e-8), n)] * model["penalty"]
        beta = solve(design.T @ design + np.diag(penalty), design.T @ y, assume_a="pos")
        intercept_update = beta[0] + q.loc[test] @ beta[2:2 + n]
        slope_update = beta[1] + q.loc[test] @ beta[2 + n:]
        predictions[model_name] = (level.loc[test] + intercept_update, slope.loc[test] + slope_update)
        fit[model_name] = {"chosen_grid": model["config"], "penalty": model["penalty"],
                           "calibration_intercept": beta[0], "calibration_slope": beta[1],
                           "feature_median": model["feature_median"], "feature_center": model["feature_center"],
                           "feature_scale": model["feature_scale"],
                           "qus_intercept": dict(zip(names, beta[2:2 + n])),
                           "qus_slope": dict(zip(names, beta[2 + n:]))}
    return predictions, fit


def run(out):
    cohort = pd.read_csv(out / "scan_cohort.csv", index_col="pid")
    measurements = pd.read_csv(out / "postscan_outcomes.csv")
    feature = pd.read_csv(out / "qus_features.csv.gz", index_col="pid")
    folds = pd.read_csv(out / "fold_assignments.csv")
    all_folds = pd.read_csv(STAGE1 / "fold_assignments.csv")
    labs, slopes, configs = [], [], []
    for (rep, fold), group in folds.groupby(["repetition", "fold"]):
        key = out / "cv" / f"r{rep:02d}_f{fold}"
        test = pd.Index(sorted(group.pid))
        train = cohort.index.difference(test).sort_values()
        stage1_test = all_folds[all_folds.repetition.eq(rep) & all_folds.fold.eq(fold)]
        assert set(test).issubset(stage1_test.pid)
        base_path = STAGE1 / "cv" / f"rep_{rep:02d}" / f"fold_{fold:02d}" / "models.json"
        base = json.loads(base_path.read_text())["models"]["M3"]
        stage2_models = json.loads((key / "models.json").read_text())
        level, slope = base_at_scan(cohort, train, base)
        params = stage2_models["clinical"]
        covariance = np.array(params["random_effect_covariance"])
        residual = params["residual_variance"]
        predictions, fitted = corrections(train, test, measurements, covariance, residual, level, slope,
                                           feature, stage2_models["qus"])
        configs.append(dict(repetition=int(rep), fold=int(fold), stage1_model=str(base_path), models=fitted))
        r = measurements[measurements.pid.isin(test)].copy()
        s = pd.read_csv(key / "slopes.csv", usecols=["pid", "reference_lmm", "reference_ols"])
        for label, (level, slope) in predictions.items():
            r[label] = r.pid.map(level) + r.t * r.pid.map(slope)
            s[label] = s.pid.map(slope)
        labs.append(r.assign(repetition=rep, fold=fold))
        slopes.append(s.assign(repetition=rep, fold=fold))
    labs = pd.concat(labs, ignore_index=True)
    slopes = pd.concat(slopes, ignore_index=True)
    combined = out / "stage1_transfer_sensitivity"
    combined.mkdir(exist_ok=True)
    labs.to_csv(combined / "oof_predictions.csv.gz", compression="gzip", index=False)
    slopes.to_csv(combined / "oof_slopes.csv.gz", compression="gzip", index=False)
    (combined / "models.json").write_text(json.dumps(configs, allow_nan=False))
    results = []
    for rep, a in labs.groupby("repetition"):
        b = slopes[slopes.repetition.eq(rep)]
        assert len(a) == len(measurements) and len(b) == len(cohort)
        for target, frame, ref in [
            ("egfr_all_postscan", a, "egfr"), ("egfr_through_2_years", a[a.t.le(2)], "egfr"),
            ("slope_lmm", b, "reference_lmm"), ("slope_ols", b, "reference_ols"),
        ]:
            weights = 1 / frame.groupby("pid").pid.transform("size").to_numpy()
            for name in NAMES:
                results.append(dict(repetition=rep, model=name, target=target) |
                               metrics(frame[ref], frame[name], weights))
    pd.DataFrame(results).to_csv(combined / "metrics_by_repetition.csv", index=False)
    rows = []
    rng = np.random.default_rng(23121)
    for target, frame, ref in [
        ("egfr_all_postscan", labs, "egfr"), ("egfr_through_2_years", labs[labs.t.le(2)], "egfr"),
        ("slope_lmm", slopes, "reference_lmm"), ("slope_ols", slopes, "reference_ols"),
    ]:
        errors = frame[list(NAMES)].subtract(frame[ref], axis="index")**2
        by_patient = errors.groupby(frame.pid).mean()
        n = len(by_patient)
        means = np.sqrt(by_patient.mean()).to_numpy()
        samples = rng.multinomial(n, np.full(n, 1 / n), size=2000)
        bootstrap = np.sqrt(samples @ by_patient.to_numpy() / n)
        for j, model in enumerate(NAMES):
            lo, hi = np.quantile(bootstrap[:, j] - bootstrap[:, 0], [.025, .975])
            rows.append(dict(model=model, target=target, patients=n, rmse=means[j], delta=means[j]-means[0],
                             low=lo, high=hi))
    pd.DataFrame(rows).to_csv(combined / "paired_rmse_differences.csv", index=False)


if __name__ == "__main__":
    run(OUT)

"""Fit clinical random-intercept and random-slope models and save held-out predictions."""
import os
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
import argparse
import json
from pathlib import Path
import warnings

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy.linalg import solve
from sklearn.model_selection import StratifiedKFold

from prepare import YEAR

FEATURES = ["baseline_egfr", "age", "female", "log_acr"]
MODELS = {"M0": [], "M1": FEATURES[:1], "M2": FEATURES[:3], "M3": FEATURES}
MODELS.update({f"without_{feature}": [f for f in FEATURES if f != feature] for feature in FEATURES})
HORIZONS = np.array([.5, 1, 1.5, 2, 3, 4, 5])
THRESHOLDS = [0, .5, 1, 1.5, 2]


def fit(data, formula):
    model = smf.mixedlm(formula, data, groups="pid", re_formula="~t")
    history = []
    for method in ("lbfgs", "bfgs", "powell"):
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            try:
                result = model.fit(reml=False, method=method, maxiter=350)
            except (ValueError, np.linalg.LinAlgError) as error:
                history.append(f"{method}: {error}")
                continue
        history.extend(f"{method}: {item.message}" for item in caught)
        if result.converged and np.isfinite(result.llf) and np.isfinite(result.fe_params).all():
            return result, method, history
        history.append(f"{method}: did not converge")
    raise RuntimeError(f"Mixed model failed for {formula}: {history}")


def transform(train, test, features):
    continuous = [feature for feature in features if feature != "female"]
    center, scale = train[continuous].mean(), train[continuous].std(ddof=0).replace(0, 1)
    train = train.assign(**((train[continuous] - center) / scale))
    test = test.assign(**((test[continuous] - center) / scale))
    return train, test, {"center": center.to_dict(), "scale": scale.to_dict()}


def parameters(result, method, history, formula, scaling, n_patients):
    fixed = result.fe_params
    ci = result.conf_int().loc[fixed.index]
    return {
        "formula": formula, "estimation": "ML", "optimizer": method, "converged": bool(result.converged),
        "warnings": history, "training_patients": n_patients, "training_rows": int(result.nobs),
        "log_likelihood": result.llf, "preprocessing": scaling,
        "fixed_effects": fixed.to_dict(), "fixed_effects_se": result.bse_fe.to_dict(),
        "fixed_effects_ci95": {key: ci.loc[key].tolist() for key in fixed.index},
        "fixed_effects_p": result.pvalues.loc[fixed.index].to_dict(),
        "random_effect_labels": ["intercept", "time_years"],
        "random_effect_covariance": result.cov_re.to_numpy().tolist(), "residual_variance": result.scale,
    }


def run_model(train, test, obs, features):
    train, test, scaling = transform(train, test, features)
    train_obs = obs[obs.pid.isin(train.index)].join(train[features], on="pid")
    formula = "egfr ~ t" + (" * (" + " + ".join(features) + ")" if features else "")
    result, method, history = fit(train_obs, formula)
    saved = parameters(result, method, history, formula, scaling, len(train))
    coeff = result.fe_params
    intercept = np.full(len(test), coeff.Intercept)
    slope = np.full(len(test), coeff.t)
    for feature in features:
        intercept += test[feature].to_numpy() * coeff[feature]
        slope += test[feature].to_numpy() * coeff[f"t:{feature}"]
    return saved, pd.DataFrame({"intercept": intercept, "slope": slope}, index=test.index)


def reference_data(cohort, obs):
    baseline = cohort.reset_index().rename(columns={"baseline_egfr": "egfr", "baseline_egfr_date": "date"})
    baseline["t"] = (pd.to_datetime(baseline.date) - pd.to_datetime(baseline.index_date)).dt.days / YEAR
    ref = pd.concat([baseline[["pid", "t", "egfr"]], obs[["pid", "t", "egfr"]]], ignore_index=True)
    ref_groups = {pid: group.sort_values("t") for pid, group in ref.groupby("pid")}

    def spaced(group):
        last, n = -np.inf, 0
        for t in group.t:
            if t - last >= (30 / YEAR - 1e-10):
                last, n = t, n + 1
        return n

    counts = pd.Series({pid: spaced(group) for pid, group in ref_groups.items()})
    return ref, ref_groups, counts


def reference_slopes(train, test, ref, groups, counts):
    subset = train.index[train.followup.ge(1) & counts.reindex(train.index).ge(3)]
    formula = "egfr ~ t"
    fit_result, method, history = fit(ref[ref.pid.isin(subset)], formula)
    params = parameters(fit_result, method, history, formula, {}, len(subset))
    beta, covariance, noise = fit_result.fe_params.to_numpy(), fit_result.cov_re.to_numpy(), fit_result.scale
    labels = []
    for pid in test.index:
        data = groups[pid]
        z = np.column_stack([np.ones(len(data)), data.t])
        y = data.egfr.to_numpy()
        # Posterior random effects using the marginal covariance identity.
        random = covariance @ solve(noise * np.eye(2) + (z.T @ z) @ covariance, z.T @ (y - z @ beta))
        ols = np.linalg.lstsq(z, y, rcond=None)[0]
        labels.append((pid, beta[1] + random[1], ols[1], counts[pid]))
    return params, pd.DataFrame(labels, columns=["pid", "reference_lmm", "reference_ols", "spaced_dates"]).set_index("pid")


def grid(repetition):
    if repetition == 1:
        yield 0, list(MODELS)
        for threshold in THRESHOLDS[1:]:
            yield threshold, ["M3"]
    else:
        yield 0, ["M0", "M1", "M2", "M3"]


def run(out, repetitions):
    cohort = pd.read_csv(out / "clinical_baselines.csv", index_col="pid")
    obs = pd.read_csv(out / "eligible_followup.csv")
    ref, groups, counts = reference_data(cohort, obs)
    assignments, records = [], []
    checkpoint = out / "cv"
    checkpoint.mkdir()
    for repetition in range(1, repetitions + 1):
        repeat_dir = checkpoint / f"rep_{repetition:02d}"
        repeat_dir.mkdir()
        splits = StratifiedKFold(n_splits=5, shuffle=True, random_state=2026 + repetition - 1)
        for fold, (training, testing) in enumerate(splits.split(cohort, cohort.kind), 1):
            print(f"Repetition {repetition}/{repetitions}, fold {fold}/5", flush=True)
            base_train, test = cohort.iloc[training], cohort.iloc[testing]
            assert not base_train.index.intersection(test.index).size
            labels_file = repeat_dir / f"fold_{fold:02d}"
            labels_file.mkdir()
            test_obs = obs[obs.pid.isin(test.index)].copy()
            ref_params, slopes = reference_slopes(base_train, test, ref, groups, counts)
            assignments.append(test[["kind"]].reset_index().assign(repetition=repetition, fold=fold))
            saved = {"reference": ref_params, "models": {}}
            slopes["persistence"] = 0.0
            test_obs["persistence"] = test_obs.pid.map(test.baseline_egfr)
            trajectory = pd.MultiIndex.from_product([test.index, HORIZONS], names=["pid", "t"]).to_frame(index=False)
            trajectory["persistence"] = trajectory.pid.map(test.baseline_egfr)
            for minimum, variants in grid(repetition):
                train = base_train[base_train.followup.ge(minimum)]
                for name in variants:
                    label = name if minimum == 0 else f"{name}_min_{minimum:g}y"
                    params, pred = run_model(train, test, obs, MODELS[name])
                    saved["models"][label] = params
                    slopes[label] = pred.slope
                    for frame in (test_obs, trajectory):
                        frame[label] = frame.pid.map(pred.intercept) + frame.t * frame.pid.map(pred.slope)
            for tab, frame, index in [("slopes", slopes, True), ("predictions", test_obs, False),
                                       ("standardized_trajectories", trajectory, False)]:
                frame.to_csv(labels_file / f"{tab}.csv", index=index)
            (labels_file / "models.json").write_text(json.dumps(saved, indent=2, allow_nan=False))
            records.extend(flatten(saved["models"] | {"reference": ref_params}, repetition, fold))
        pd.concat(assignments, ignore_index=True).to_csv(out / "fold_assignments.csv", index=False)
    final, _ = run_model(cohort, cohort, obs, FEATURES)
    (out / "full_cohort_M3.json").write_text(json.dumps(final, indent=2, allow_nan=False))
    records.extend(flatten({"M3_full_cohort": final}, 0, 0))
    flat = pd.DataFrame(records)
    fixed = flat[flat.term.notna()].drop(columns=["var_intercept", "var_slope", "cov_intercept_slope", "residual_variance"])
    variance = flat[flat.term.isna()].drop(columns=["term", "estimate", "se", "ci_lower", "ci_upper", "p_value"])
    fixed.to_csv(out / "fixed_effects.csv", index=False)
    variance.to_csv(out / "random_effects_and_residual_variance.csv", index=False)


def flatten(models, repetition, fold):
    for name, params in models.items():
        keys = dict(model=name, repetition=repetition, fold=fold)
        yield keys | dict(var_intercept=params["random_effect_covariance"][0][0],
                          cov_intercept_slope=params["random_effect_covariance"][0][1],
                          var_slope=params["random_effect_covariance"][1][1],
                          residual_variance=params["residual_variance"])
        for term, estimate in params["fixed_effects"].items():
            lo, hi = params["fixed_effects_ci95"][term]
            yield keys | dict(term=term, estimate=estimate, se=params["fixed_effects_se"][term],
                              ci_lower=lo, ci_upper=hi, p_value=params["fixed_effects_p"][term])


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    parser.add_argument("--repetitions", type=int, default=10)
    args = parser.parse_args()
    run(args.output, args.repetitions)

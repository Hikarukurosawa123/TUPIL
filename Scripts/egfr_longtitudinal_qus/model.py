"""Nested validation of scan-time clinical models with longitudinal QUS effects."""
import os
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"

import argparse
from itertools import product
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy.linalg import cholesky, solve, solve_triangular
from scipy.stats import rankdata
from sklearn.model_selection import KFold

from prepare import MAPS, OUT

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "egfr_longtitudinal_all"))
from fit import fit, parameters, FEATURES

COMBINED = ("overall", "per_modality")
TARGETS = ("slope", "category")
SIZES = (1, 3, 5)
PENALTIES = (1., 10., 100.)
GRID = tuple(product(TARGETS, COMBINED + tuple(MAPS), SIZES))
PREDICTORS = ["intercept", "t", *FEATURES, *(f"t:{name}" for name in FEATURES)]
MODELS = ("persistence", "clinical", "qus_nested_egfr", "qus_nested_slope")


def name(config):
    target, scope, count = config
    return f"{target}__{scope}__top{count}"


def x_design(z, t):
    t = np.asarray(t, dtype=float)
    return np.column_stack((np.ones(len(t)), t, z, z * t[:, None]))


def patient_error(y, pred, person):
    error = np.bincount(person, weights=(y - pred)**2) / np.bincount(person)
    return error.sum()


def clinical_data(cohort, train, test):
    continuous = [c for c in FEATURES if c != "female"]
    train_data = cohort.loc[train, FEATURES].copy()
    test_data = cohort.loc[test, FEATURES].copy()
    median = train_data[continuous].median()
    train_data[continuous] = train_data[continuous].fillna(median)
    test_data[continuous] = test_data[continuous].fillna(median)
    center = train_data[continuous].mean()
    scale = train_data[continuous].std(ddof=0).replace(0, 1)
    for d in (train_data, test_data):
        d[continuous] = (d[continuous] - center) / scale
    return train_data, test_data, {"median": median.to_dict(), "center": center.to_dict(), "scale": scale.to_dict()}


def ref_slopes(obs, model, ids):
    beta = model.fe_params.to_numpy()
    covariance, residual = model.cov_re.to_numpy(), model.scale
    values = []
    for pid in ids:
        d = obs[pid]
        z = np.column_stack([np.ones(len(d)), d.t])
        y = d.egfr.to_numpy()
        random = covariance @ solve(residual * np.eye(2) + (z.T @ z) @ covariance, z.T @ (y - z @ beta))
        ols = np.linalg.lstsq(z, y, rcond=None)[0][1]
        values.append((beta[1] + random[1], ols))
    return np.array(values)


class Split:
    def __init__(self, cohort, obs, feature, train, test):
        assert len(set(train) & set(test)) == 0
        self.train, self.test = train, test
        self.train_clinical, self.test_clinical, self.scaling = clinical_data(cohort, train, test)
        self.train_obs = pd.concat([obs[pid] for pid in train], ignore_index=True).join(self.train_clinical, on="pid")
        self.test_obs = pd.concat([obs[pid] for pid in test], ignore_index=True)
        self.test_y = self.test_obs.egfr.to_numpy()
        self.test_person = test.get_indexer(self.test_obs.pid)
        self.test_t = self.test_obs.t.to_numpy()

        formula = "egfr ~ t * (" + " + ".join(FEATURES) + ")"
        model, method, warnings = fit(self.train_obs, formula)
        self.covariance, self.residual = model.cov_re.to_numpy(), model.scale
        self.clinical_params = parameters(model, method, warnings, formula, self.scaling, len(train))
        self.test_x = x_design(self.test_clinical.to_numpy()[self.test_person], self.test_t)
        self.test_slope_x = self.slope_design(self.test_clinical.to_numpy())
        self.base_beta = model.fe_params.reindex(["Intercept", *PREDICTORS[1:]]).to_numpy()
        reference, method, warnings = fit(self.train_obs, "egfr ~ t")
        self.reference_params = parameters(reference, method, warnings, "egfr ~ t", {}, len(train))
        self.train_ref = ref_slopes(obs, reference, train)[:, 0]
        self.test_ref = ref_slopes(obs, reference, test)
        self.radiomics(feature)
        self.gls_setup()
        assert np.allclose(solve(self.XtX, self.Xty, assume_a="pos"), self.base_beta, atol=1e-6)

    @staticmethod
    def slope_design(z):
        return np.column_stack((np.zeros(len(z)), np.ones(len(z)), np.zeros_like(z), z))

    def radiomics(self, feature):
        candidates = feature.loc[self.train]
        enough = candidates.notna().sum().ge(.9 * len(candidates))
        candidates = candidates.loc[:, enough]
        median = candidates.median()
        filled = candidates.fillna(median)
        scale = filled.std(ddof=0)
        keep = scale.gt(0) & np.isfinite(scale)
        filled = filled.loc[:, keep]
        self.feature_names = filled.columns
        self.feature_median = median.loc[self.feature_names]
        self.feature_center = filled.mean()
        self.feature_scale = scale.loc[self.feature_names]
        self.q_train = ((filled - self.feature_center) / self.feature_scale).to_numpy()
        self.q_test = ((feature.loc[self.test, self.feature_names].fillna(self.feature_median)
                        - self.feature_center) / self.feature_scale).to_numpy()
        self.map = self.feature_names.str.split("__").str[0].to_numpy()
        np.testing.assert_allclose(np.mean(self.q_train, axis=0), 0, atol=1e-8)
        ranks = rankdata(self.q_train, axis=0)
        ranks -= ranks.mean(axis=0)
        ranks /= np.linalg.norm(ranks, axis=0)
        category = sum((self.train_ref < edge).astype(int) for edge in (-1, -3, -5))
        self.ranking = {}
        self.rho = {}
        for target, values in [("slope", self.train_ref), ("category", category)]:
            yrank = rankdata(values)
            yrank -= yrank.mean()
            r = ranks.T @ yrank / np.linalg.norm(yrank)
            self.rho[target] = r
            self.ranking[target] = np.argsort(-np.abs(r), kind="stable")

    def gls_setup(self):
        x_rows, y_rows, z_rows, pid_rows = [], [], [], []
        for i, (pid, d) in enumerate(self.train_obs.groupby("pid", sort=False)):
            assert pid == self.train[i]
            t = d.t.to_numpy()
            z = np.column_stack([np.ones(len(t)), t])
            marginal = z @ self.covariance @ z.T + self.residual * np.eye(len(t))
            lower = cholesky(marginal, lower=True)
            x_rows.append(solve_triangular(lower, x_design(d[FEATURES].to_numpy(), t), lower=True))
            z_rows.append(solve_triangular(lower, z, lower=True))
            y_rows.append(solve_triangular(lower, d.egfr.to_numpy(), lower=True))
            pid_rows.extend([i] * len(t))
        self.white_x, self.white_z, self.white_y = map(np.concatenate, (x_rows, z_rows, y_rows))
        self.train_person = np.array(pid_rows)
        self.XtX = self.white_x.T @ self.white_x
        self.Xty = self.white_x.T @ self.white_y

    def pick(self, target, group, top):
        order = self.ranking[target]
        if group != "overall":
            order = order[self.map[order] == group]
        out = []
        for i in order:
            if not out or np.max(np.abs(self.q_train[:, out].T @ self.q_train[:, i])) < .9 * len(self.train):
                out.append(i)
                if len(out) == top:
                    break
        assert len(out) == top
        return out

    def select(self, config):
        target, scope, top = config
        if scope == "per_modality":
            return [i for m in MAPS for i in self.pick(target, m, top)]
        return self.pick(target, scope, top)

    def model_grid(self, config, penalty_values):
        selected = self.select(config)
        train_q = self.q_train[self.train_person][:, selected]
        wq = np.column_stack((train_q * self.white_z[:, :1], train_q * self.white_z[:, 1:]))
        full = np.column_stack((self.white_x, wq))
        gram, rhs = full.T @ full, full.T @ self.white_y
        sizes = np.maximum(np.diag(self.covariance), 1e-8)
        diag = np.r_[np.zeros(len(PREDICTORS)), np.repeat(1 / sizes, len(selected))]
        n_terms = len(PREDICTORS)
        test_q = self.q_test[:, selected]
        for penalty in penalty_values:
            beta = solve(gram + np.diag(penalty * diag), rhs, assume_a="pos")
            intercept = self.test_clinical.to_numpy() @ beta[2:2 + len(FEATURES)] + beta[0]
            intercept += test_q @ beta[n_terms:n_terms + len(selected)]
            slope = self.test_slope_x @ beta[:n_terms] + test_q @ beta[n_terms + len(selected):]
            pred = intercept[self.test_person] + self.test_t * slope[self.test_person]
            yield penalty, beta, pred, slope, selected

    def baseline(self):
        return self.test_x @ self.base_beta, self.test_slope_x @ self.base_beta

    def describe(self, config, penalty, beta, indices):
        names = self.feature_names[indices]
        target, _, _ = config
        return dict(penalty=penalty, features=names.tolist(), rho=self.rho[target][indices].tolist(),
                    feature_median=self.feature_median.loc[names].to_dict(),
                    feature_center=self.feature_center.loc[names].to_dict(),
                    feature_scale=self.feature_scale.loc[names].to_dict(),
                    clinical_beta=dict(zip(PREDICTORS, beta[:len(PREDICTORS)])),
                    qus_intercept_beta=dict(zip(names, beta[len(PREDICTORS):len(PREDICTORS) + len(names)])),
                    qus_slope_beta=dict(zip(names, beta[len(PREDICTORS) + len(names):])))


def tune(cohort, obs, x, train, seed):
    scores = {(config, penalty): np.zeros(2) for config in GRID for penalty in PENALTIES}
    for a, b in KFold(n_splits=3, shuffle=True, random_state=seed).split(train):
        split = Split(cohort, obs, x, train[a], train[b])
        for config in GRID:
            for penalty, _, pred, slope, _ in split.model_grid(config, PENALTIES):
                scores[config, penalty] += [
                    patient_error(split.test_y, pred, split.test_person),
                    np.sum((split.test_ref[:, 0] - slope)**2),
                ]
    return scores


def run(out, first=1, last=10):
    cohort = pd.read_csv(out / "scan_cohort.csv", index_col="pid")
    measurements = pd.read_csv(out / "postscan_outcomes.csv")
    obs = {pid: d for pid, d in measurements.groupby("pid")}
    features = pd.read_csv(out / "qus_features.csv.gz", index_col="pid")
    assignments = pd.read_csv(out / "fold_assignments.csv")
    root = out / "cv"
    root.mkdir(exist_ok=True)
    assignments = assignments[assignments.repetition.between(first, last)]
    for (repetition, fold), holdout in assignments.groupby(["repetition", "fold"]):
        key = root / f"r{repetition:02d}_f{fold}"
        if (key / "inner_scores.csv").exists():
            continue
        print(f"Repetition {repetition}, fold {fold}", flush=True)
        test = pd.Index(sorted(holdout.pid))
        train = cohort.index.difference(test).sort_values()
        scores = tune(cohort, obs, features, train, 30000 + repetition * 100 + fold)
        combined = [key for key in scores if key[0][1] in COMBINED]
        nested = {f"qus_nested_{label}": min(combined, key=lambda key: scores[key][i])
                  for i, label in enumerate(("egfr", "slope"))}
        penalties = {config: min(PENALTIES, key=lambda p: scores[config, p][0]) for config in GRID}

        split = Split(cohort, obs, features, train, test)
        pred_data = {"persistence": split.test_obs.pid.map(cohort.baseline_egfr).to_numpy()}
        slope_data = {"persistence": np.zeros(len(test)),
                      "reference_lmm": split.test_ref[:, 0], "reference_ols": split.test_ref[:, 1]}
        pred_data["clinical"], slope_data["clinical"] = split.baseline()
        fitted = {}
        for config, penalty in penalties.items():
            keys = {penalty: [name(config)]}
            for label, chosen in nested.items():
                if chosen[0] == config:
                    keys.setdefault(chosen[1], []).append(label)
            for p, beta, pred, slope, selected in split.model_grid(config, keys):
                params = split.describe(config, p, beta, selected)
                for label in keys[p]:
                    pred_data[label], slope_data[label] = pred, slope
                    fitted[label] = params | {"config": name(config)}
        key.mkdir(exist_ok=True)
        (key / "models.json").write_text(json.dumps(
            {"clinical": split.clinical_params, "reference": split.reference_params, "qus": fitted}, allow_nan=False))
        split.test_obs.join(pd.DataFrame(pred_data, index=split.test_obs.index)).to_csv(key / "oof.csv", index=False)
        pd.DataFrame(slope_data, index=test.rename("pid")).to_csv(key / "slopes.csv")
        pd.DataFrame([{"configuration": name(config), "penalty": p, "inner_patient_egfr_rmse": np.sqrt(v[0] / len(train)),
                       "inner_slope_rmse": np.sqrt(v[1] / len(train))}
                      for (config, p), v in scores.items()]).to_csv(key / "inner_scores.csv", index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--first", type=int, default=1)
    parser.add_argument("--last", type=int, default=10)
    args = parser.parse_args()
    run(OUT, args.first, args.last)


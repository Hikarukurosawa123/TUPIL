"""Independently reconcile retained measurements and saved fold predictions."""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from prepare import MATCH_DAYS, date, read_csv


def verify(out):
    cohort = pd.read_csv(out / "clinical_baselines.csv", index_col="pid",
                         parse_dates=["index_date", "baseline_egfr_date", "acr_date",
                                      "tx", "next_tx", "graft_failure", "rrt"])
    obs = pd.read_csv(out / "eligible_followup.csv", parse_dates=["date"])
    audit = pd.read_csv(out / "lab_audit.csv.gz", low_memory=False, parse_dates=["date", "tx"])
    assert audit.groupby(["source", "source_row", "lab"]).size().eq(1).all()
    eligible = audit[audit.disposition.eq("Eligible")].set_index(["source", "source_row", "lab"])
    assigned = [("baseline_egfr", "egfr", "baseline_egfr_date"), ("acr", "acr", "acr_date")]
    for field, lab, date_field in assigned:
        keys = pd.MultiIndex.from_arrays([cohort[f"{field}_source"], cohort[f"{field}_row"],
                                         [lab] * len(cohort)])
        d = eligible.reindex(keys)
        assert d.pid.tolist() == list(cohort.index)
        assert d.series.tolist() == cohort.series.tolist()
        assert np.allclose(d.value, cohort[field])
        assert np.array_equal(d.date, cohort[date_field])
        assert (cohort.index_date - cohort[date_field]).dt.days.between(0, MATCH_DAYS).all()
    point_keys = pd.MultiIndex.from_arrays([obs.source, obs.source_row, ["egfr"] * len(obs)])
    provenance = eligible.reindex(point_keys)
    assert np.array_equal(provenance.pid, obs.pid)
    assert np.array_equal(provenance.series, obs.series)
    assert np.array_equal(provenance.date, obs.date)
    assert np.allclose(provenance.value, obs.egfr)
    assert len(obs) == cohort.n_targets.sum()
    assert obs.groupby("pid").size().sort_index().equals(cohort.n_targets.sort_index())
    assert set(cohort.index) == set(obs.pid)
    complement = obs.join(cohort.drop(columns=["series"]), on="pid")
    assert complement.date.gt(complement.index_date).all()
    for cutoff in ("next_tx", "graft_failure", "rrt"):
        assert (complement[cutoff].isna() | complement.date.lt(complement[cutoff])).all()
    transplant = complement[complement.kind.eq("transplant")]
    assert transplant.date.ge(transplant.tx + pd.DateOffset(months=3)).all()
    assert cohort.loc[cohort.tx.notna(), "acr_date"].ge(
        cohort.tx.dropna() + pd.DateOffset(months=3)).all()

    history = read_csv("*ExportTransplantHist*").iloc[:, [0, 3]].copy()
    history.columns = ["pid", "transplant_date"]
    history["transplant_date"] = date(history.transplant_date)
    transplanted_labs = audit[audit.kind.eq("transplant") & audit.pid.notna()].copy()
    source_tx = transplanted_labs[["pid", "tx"]].dropna().rename(columns={"tx": "transplant_date"})
    earliest = pd.concat([history, source_tx]).groupby("pid").transplant_date.min()
    native = complement[complement.kind.eq("native")]
    first_tx = native.pid.map(earliest)
    assert (first_tx.isna() | native.date.lt(first_tx)).all()
    assert (first_tx.isna() | native.baseline_egfr_date.lt(first_tx)).all()
    joined = transplanted_labs.merge(history, on="pid", how="inner")
    conflicted = joined[joined.tx.lt(joined.transplant_date) & joined.date.ge(joined.transplant_date)]
    conflicted = conflicted.drop_duplicates(["source", "source_row"])
    assert not conflicted.disposition.eq("Eligible").any()
    conflicts = (conflicted.groupby(["lab", "disposition"])
                 .agg(rows=("source_row", "size"), patients=("pid", "nunique")).reset_index())
    conflicts.to_csv(out / "transplant_episode_conflicts.csv", index=False)
    assignments = pd.read_csv(out / "fold_assignments.csv")
    n_predictions = n_fits = 0
    status = []
    for (repetition, fold), test in assignments.groupby(["repetition", "fold"]):
        test = cohort.loc[test.pid]
        root = out / "cv" / f"rep_{repetition:02d}" / f"fold_{fold:02d}"
        predictions = pd.read_csv(root / "predictions.csv")
        training = cohort.loc[~cohort.index.isin(test.index)]
        assert set(predictions.pid) == set(test.index)
        assert len(predictions) == test.n_targets.sum()
        assert np.allclose(predictions.persistence, predictions.pid.map(test.baseline_egfr))
        n_predictions += len(predictions)
        saved = json.loads((root / "models.json").read_text())
        record = {"repetition": int(repetition), "fold": int(fold),
                  "models": len(saved["models"]), "n_held_out": len(test)}
        status.append(record)
        for name, model in saved["models"].items():
            assert model["converged"]
            minimum = float(name.split("_min_")[-1][:-1]) if "_min_" in name else 0
            train = training[training.followup.ge(minimum)]
            assert model["training_patients"] == len(train)
            assert model["training_rows"] == train.n_targets.sum()
            scaling = model["preprocessing"]
            beta = model["fixed_effects"]
            intercept = pd.Series(beta["Intercept"], index=test.index)
            slope = pd.Series(beta["t"], index=test.index)
            for feature, center in scaling["center"].items():
                scale = scaling["scale"][feature]
                assert np.isclose(train[feature].mean(), center)
                assert np.isclose(train[feature].std(ddof=0), scale)
                z = (test[feature] - center) / scale
                intercept += beta[feature] * z
                slope += beta[f"t:{feature}"] * z
            if "female" in beta:
                intercept += beta["female"] * test.female
                slope += beta["t:female"] * test.female
            reproduced = predictions.pid.map(intercept) + predictions.t * predictions.pid.map(slope)
            assert np.allclose(reproduced, predictions[name], rtol=0, atol=1e-8)
            covariance = np.array(model["random_effect_covariance"])
            assert np.linalg.eigvalsh(covariance).min() > 0
            assert model["residual_variance"] > 0
            n_fits += 1
    report = {
        "verified_same_episode_baselines": len(cohort),
        "verified_same_episode_targets": len(obs),
        "verified_transplant_targets": len(transplant),
        "verified_all_repeats_patient_date_predictions": n_predictions,
        "verified_prediction_fits": n_fits,
        "folds": status,
    }
    (out / "verification.json").write_text(json.dumps(report, indent=2))
    print(json.dumps({k: v for k, v in report.items() if k != "folds"}), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    verify(parser.parse_args().output)


"""Refit M3 retaining all distinct eligible dates after the same clinical baseline."""
import argparse
import json
from pathlib import Path

import pandas as pd

from evaluate import read_cv, trajectory_metrics
from fit import FEATURES, run_model
from prepare import YEAR


def run(out):
    target = out / "all_distinct_dates_sensitivity"
    target.mkdir(exist_ok=True)
    assert not (target / 'models.json').exists()
    cohort = pd.read_csv(out / "clinical_baselines.csv", index_col="pid", parse_dates=["index_date"])
    primary = pd.read_csv(out / "eligible_followup.csv", parse_dates=["date"])
    assignments = pd.read_csv(out / "fold_assignments.csv").query("repetition == 1")
    audit = pd.read_csv(out / "lab_audit.csv.gz", low_memory=False, parse_dates=["date"])
    keep = audit.disposition.eq("Later day of consecutive-day cluster") & audit.lab.eq("egfr")
    extra = audit.loc[keep].merge(cohort[["series", "index_date"]], on=["pid", "series"], validate="m:1")
    extra = extra[extra.date.gt(extra.index_date)].copy()
    extra["t"] = (extra.date - extra.index_date).dt.days / YEAR
    extra = extra.rename(columns={"value": "egfr"})[primary.columns]
    full = pd.concat([primary.assign(in_case8_series=True), extra.assign(in_case8_series=False)], ignore_index=True)
    assert not full.duplicated(["pid", "date"]).any()
    baseline = cohort.reindex(full.pid).reset_index()
    assert full.series.eq(baseline.series).all()
    for column in ("next_tx", "graft_failure", "rrt"):
        dates = pd.to_datetime(baseline[column])
        assert (dates.isna() | full.date.lt(dates)).all()
    held, fitted = [], {}
    for fold in sorted(assignments.fold.unique()):
        testing = assignments.loc[assignments.fold.eq(fold), "pid"]
        train, test = cohort[~cohort.index.isin(testing)], cohort.loc[testing]
        params, estimates = run_model(train, test, full, FEATURES)
        fitted[str(fold)] = params
        primary_params = json.loads((out / "cv" / "rep_01" / f"fold_{fold:02d}" / "models.json").read_text())["models"]["M3"]
        d = full[full.pid.isin(testing)].join(test[["baseline_egfr", "kind"]], on="pid").copy()
        d["M3_all_dates"] = d.pid.map(estimates.intercept) + d.t * d.pid.map(estimates.slope)
        d["persistence"] = d.baseline_egfr
        scaling = primary_params["preprocessing"]
        z = (test[list(scaling["center"])] - pd.Series(scaling["center"])) / pd.Series(scaling["scale"])
        z["female"] = test.female
        beta = primary_params["fixed_effects"]
        intercept = beta["Intercept"] + sum(z[c] * beta[c] for c in FEATURES)
        slope = beta["t"] + sum(z[c] * beta[f"t:{c}"] for c in FEATURES)
        d["M3_case8"] = d.pid.map(intercept) + d.t * d.pid.map(slope)
        held.append(d.assign(fold=fold))
    predictions = pd.concat(held, ignore_index=True)
    match = predictions[predictions.in_case8_series].assign(date=lambda d: d.date.dt.strftime("%Y-%m-%d"))
    reference = match.merge(read_cv(out, "predictions", 1)[["pid", "date", "M3"]], on=["pid", "date"],
                            validate="1:1")
    assert len(reference) == len(primary)
    assert (reference.M3_case8 - reference.M3).abs().max() < 1e-8
    tables = []
    for label, d in [("all_distinct_dates", predictions),
                     ("case8_dates_only", predictions[predictions.in_case8_series])]:
        tables.extend(trajectory_metrics(d, {"evaluation_targets": label},
                                         ["persistence", "M3_case8", "M3_all_dates"]))
    pd.DataFrame(tables).to_csv(target / "metrics.csv", index=False)
    predictions.to_csv(target / "out_of_fold_predictions.csv", index=False)
    (target / "models.json").write_text(json.dumps(fitted, indent=2, allow_nan=False))
    print(f"All-distinct-date fit added {len(extra)} outcomes for {extra.pid.nunique()} patients.", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    run(parser.parse_args().output)



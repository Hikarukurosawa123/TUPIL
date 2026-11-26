"""Build one clinical baseline and one kidney trajectory per patient."""
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
INPUT = ROOT / "All_csv"
CASE8 = ROOT.parent / "POCUS_Kidney_Analysis" / "Analysis_scripts" / "csv" / "Clinical Variables"
MATCH_DAYS = 90
YEAR = 365.25


def read_csv(pattern):
    return pd.read_csv(next(INPUT.glob(pattern)), low_memory=False)


def date(values):
    return pd.to_datetime(values, format="%d/%m/%Y", errors="raise")


def as_time(values):
    value = values.astype(str)
    value = value.where(value.str.count(":") != 1, value + ":00")
    return pd.to_timedelta(value, errors="coerce")


def source_data():
    hist = read_csv("*ExportTransplantHist*").iloc[:, [0, 3, 10, 11]].copy()
    hist.columns = ["pid", "tx", "failed", "failure"]
    hist["tx"], hist["failure"] = date(hist.tx), date(hist.failure)
    assert hist[["pid", "tx"]].notna().all().all()
    assert not hist.duplicated(["pid", "tx"]).any()

    native = read_csv("*DataEntryNativeKidne*")
    native.columns = ["pid", "event", "repeat", "date", "time", "creatinine",
                      "egfr", "egfr_krt", "albumin", "urine_creatinine", "acr"]
    native["date"] = date(native.date)
    native["time"] = as_time(native.time)
    native["source_row"] = native.index + 2
    native["tx"] = pd.NaT
    native["export_end"] = pd.NaT

    book = pd.read_excel(next(INPUT.glob("PROSPBX_DataRequest*.xlsx")), sheet_name=None)
    bw, urine = book["PTX BW"].copy(), book["PTX URINE"].copy()
    bw.columns = ["pid", "tx", "export_end", "date", "time", "creatinine", "egfr", "egfr_krt"]
    urine.columns = ["pid", "tx", "export_end", "date", "time", "albumin",
                     "urine_creatinine", "acr_reported", "acr_calculated"]
    urine["acr"] = urine.acr_calculated.combine_first(urine.acr_reported)
    sets = {}
    for lab, transplant_labs, sheet in [("egfr", bw, "PTX BW"), ("acr", urine, "PTX URINE")]:
        transplant_labs["source_row"] = transplant_labs.index + 2
        cols = ["pid", "date", "time", "tx", "export_end", "source_row", lab]
        labs = pd.concat([transplant_labs[cols].assign(kind="transplant", source=sheet),
                          native[cols].assign(kind="native", source="Native REDCap")], ignore_index=True)
        labs[lab] = pd.to_numeric(labs[lab], errors="coerce")
        labs["lab"] = lab
        sets[lab] = labs.rename(columns={lab: "value"})
    return hist, sets


def rrt_data(source=None):
    if source is None:
        source = next(CASE8.glob("case8_no_min_followup_RRT_adjusted_2p5MHz_only_egfr_with_RRT_timing.xlsx"))
    sheets = pd.read_excel(source, sheet_name=["Patient_Summary", "Excluded_Patients"])
    cols = ["patient_id", "RRT Status", "RRT Initiation Date"]
    records = pd.concat([frame[cols] for frame in sheets.values()]).rename(
        columns=dict(zip(cols, ["pid", "status", "rrt"])))
    records["rrt"] = pd.to_datetime(records.rrt)
    dates = {pid: list(group.rrt.dropna().unique()) for pid, group in records.groupby("pid")}
    undated = set(records.loc[records.status.eq("Baseline RRT") & records.rrt.isna(), "pid"])
    undated -= {pid for pid, value in dates.items() if value}
    return dates, undated


def demographics():
    sex = read_csv("PROSPBX*DemographicBaseline*").set_index("participant_id").gender
    age = read_csv("PROSPBX*VisitReview*").dropna(subset=["visit_date", "visit_age"])
    age["visit_date"] = date(age.visit_date)
    age = age.sort_values("visit_date").drop_duplicates("participant_id").set_index("participant_id")
    return sex, age


def clean(labs, history, rrt_dates, undated):
    episodes = history.sort_values(["pid", "tx"]).copy()
    episodes["next_tx"] = episodes.groupby("pid").tx.shift(-1)
    episodes["first_tx"] = episodes.groupby("pid").tx.transform("min")
    native_cutoff = episodes.drop_duplicates("pid").set_index("pid").first_tx

    def episode_rrt(row):
        return min((pd.Timestamp(d) for d in rrt_dates.get(row.pid, [])
                    if pd.isna(row.tx) or d >= row.tx), default=pd.NaT)

    episodes["rrt"] = episodes.apply(episode_rrt, axis=1)
    labs = labs.merge(episodes[["pid", "tx", "next_tx", "failure", "rrt", "failed"]],
                      how="left", on=["pid", "tx"], validate="m:1")
    nat = labs.kind.eq("native")
    labs.loc[nat, "next_tx"] = labs.loc[nat, "pid"].map(native_cutoff)
    first_rrt = {pid: min(dates) for pid, dates in rrt_dates.items() if dates}
    labs.loc[nat, "rrt"] = labs.loc[nat, "pid"].map(first_rrt)
    labs["washout_end"] = labs.tx + pd.DateOffset(months=3)
    labs["series"] = labs.tx.dt.strftime("%Y-%m-%d").where(~nat, "native")
    minimum = labs.lab.eq("egfr").astype(int)
    missing = (labs.pid.isna() | labs.date.isna() | ~np.isfinite(labs.value)
               | ((labs.value <= 0) & minimum.eq(1)) | (labs.value < 0))
    checks = [
        (missing, "Missing or invalid lab"),
        (~nat & labs.failed.isna(), "Workbook episode absent from transplant history"),
        (labs.date.ge(labs.next_tx), "On or after subsequent transplant"),
        (~nat & labs.date.gt(labs.export_end), "After workbook observation end"),
        (~nat & labs.date.ge(labs.failure), "On or after graft failure"),
        (labs.pid.isin(undated), "Known RRT with unavailable initiation date"),
        (labs.date.ge(labs.rrt), "On or after recorded RRT initiation"),
        (~nat & labs.date.lt(labs.washout_end), "Before three calendar months post transplant"),
    ]
    labs["disposition"] = np.select([mask.fillna(False) for mask, _ in checks],
                                    [reason for _, reason in checks], default="Eligible")
    labs = labs.sort_values(["pid", "series", "date", "time", "source_row"], kind="stable")
    keys = ["pid", "series", "date"]
    keep = labs.disposition.eq("Eligible")
    dupe = labs.loc[keep].duplicated(keys)
    labs.loc[dupe.index[dupe], "disposition"] = "Later same-day measurement"
    return labs


def prepare(out):
    hist, raw = source_data()
    rrt, undated = rrt_data()
    sex, age = demographics()
    egfr, acr = [clean(raw[lab], hist, rrt, undated) for lab in ("egfr", "acr")]

    # Case 8 proxy for concentrated lab sampling during inpatient stays.
    daily = egfr.loc[egfr.disposition.eq("Eligible")].copy()
    episode_keys = ["pid", "series"]
    daily["run"] = daily.groupby(episode_keys).date.diff().dt.days.ne(1).cumsum()
    drop = daily.run.map(daily.run.value_counts()).ge(3) & daily.run.duplicated()
    egfr.loc[daily.index[drop], "disposition"] = "Later day of consecutive-day cluster"

    acr_groups = dict(tuple(acr[acr.disposition.eq("Eligible")].groupby(episode_keys)))
    baselines, outcomes, cohort_flow = [], [], []
    eligible = egfr[egfr.disposition.eq("Eligible")]
    for pid, all_egfr in eligible.groupby("pid"):
        if pid not in sex or sex[pid] not in ("F", "M") or pid not in age.index:
            cohort_flow.append((pid, "Missing age or sex"))
            continue
        choices = []
        for series, blood in all_egfr.groupby("series"):
            urine = acr_groups.get((pid, series))
            if urine is None:
                continue
            bd, ud = pd.DatetimeIndex(blood.date), pd.DatetimeIndex(urine.date)
            days = bd.union(ud)
            bidx, uidx = bd.searchsorted(days, side="right") - 1, ud.searchsorted(days, side="right") - 1
            good = ((bidx >= 0) & (uidx >= 0) & ((days - bd[bidx]).days <= MATCH_DAYS)
                    & ((days - ud[uidx]).days <= MATCH_DAYS))
            options = np.flatnonzero(good)
            if len(options):
                i = options[0]
                choices.append((days[i], series, blood.iloc[bidx[i]], urine.iloc[uidx[i]], blood))
        if not choices:
            cohort_flow.append((pid, "No eligible same-episode ACR/eGFR pair"))
            continue
        index, series, baseline, urine, blood = min(choices, key=lambda x: (x[0], x[1]))
        follow = blood[blood.date.gt(index)].copy()
        if len(follow) < 2:
            cohort_flow.append((pid, "First clinical baseline has fewer than two follow-up dates"))
            continue
        age_anchor = age.loc[pid]
        baseline_age = age_anchor.visit_age + (index - age_anchor.visit_date).days / YEAR
        if not 18 <= baseline_age <= 120:
            cohort_flow.append((pid, "Age at index outside adult analysis range"))
            continue
        span = (follow.date.max() - index).days / YEAR
        baselines.append(dict(pid=pid, kind=baseline.kind, series=series, tx=baseline.tx,
                              next_tx=baseline.next_tx, graft_failure=baseline.failure,
                              rrt=baseline.rrt, index_date=index, baseline_egfr_date=baseline.date,
                              acr_date=urine.date, age=baseline_age, female=int(sex[pid] == "F"),
                              acr=urine.value, log_acr=np.log1p(urine.value),
                              baseline_egfr=baseline.value, followup=span, n_targets=len(follow),
                              baseline_egfr_source=baseline.source, baseline_egfr_row=baseline.source_row,
                              acr_source=urine.source, acr_row=urine.source_row))
        cohort_flow.append((pid, "Included"))
        follow["t"] = (follow.date - index).dt.days / YEAR
        outcomes.append(follow[["pid", "series", "date", "t", "value", "source", "source_row"]])
    cohort = pd.DataFrame(baselines).set_index("pid")
    obs = pd.concat(outcomes, ignore_index=True).rename(columns={"value": "egfr"})
    flow = pd.DataFrame(cohort_flow, columns=["pid", "status"])
    no_egfr = sorted(set(raw["egfr"].loc[raw["egfr"].value.gt(0), "pid"]) - set(flow.pid))
    flow = pd.concat([flow, pd.DataFrame({"pid": no_egfr, "status": "No eligible eGFR after filtering"})])
    labs = pd.concat([egfr, acr], ignore_index=True)
    counts = (labs.groupby(["lab", "kind", "disposition"], dropna=False)
              .agg(rows=("disposition", "size"), patients=("pid", "nunique")).reset_index())
    provenance = cohort.reindex(obs.pid).reset_index()
    assert np.array_equal(obs.series, provenance.series)
    assert obs.date.gt(provenance.index_date).all()
    assert (provenance.next_tx.isna() | obs.date.lt(provenance.next_tx)).all()
    assert (provenance.rrt.isna() | obs.date.lt(provenance.rrt)).all()
    assert (provenance.graft_failure.isna() | obs.date.lt(provenance.graft_failure)).all()
    assert (provenance.tx.isna() | obs.date.ge(provenance.tx + pd.DateOffset(months=3))).all()
    assert cohort.baseline_egfr_date.le(cohort.index_date).all()
    assert cohort.acr_date.le(cohort.index_date).all()
    assert not obs.duplicated(["pid", "date"]).any()
    assert cohort.n_targets.sum() == len(obs)
    cohort.to_csv(out / "clinical_baselines.csv")
    obs.to_csv(out / "eligible_followup.csv", index=False)
    flow.to_csv(out / "patient_disposition.csv", index=False)
    counts.to_csv(out / "lab_filter_counts.csv", index=False)
    labs.drop(columns=["time"]).to_csv(out / "lab_audit.csv.gz", index=False, compression="gzip")
    return cohort, obs


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True)
    cohort, obs = prepare(args.output)
    print(cohort.groupby("kind").agg(patients=("series", "size"), targets=("n_targets", "sum")))

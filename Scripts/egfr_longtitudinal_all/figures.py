"""Generate review figures and the results summary from completed metrics."""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import precision_recall_curve, roc_curve

from fit import THRESHOLDS

BLUE, ROSE, GRAY = "#4C9AC6", "#E84A75", "#555555"
plt.rcParams.update({"font.family": "Arial", "font.size": 8, "axes.spines.top": False,
                     "axes.spines.right": False, "axes.linewidth": .6, "axes.labelsize": 8,
                     "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7})
COMPARATORS = ["persistence", "M0", "M1", "M2", "M3"]


def render(out):
    def read(table):
        return pd.read_csv(out / f"{table}.csv")

    figures = out / "figures"
    figures.mkdir(exist_ok=True)

    def export(fig, name):
        fig.savefig(figures / f"{name}.png", dpi=300, bbox_inches="tight")
        plt.close(fig)

    pred = read("out_of_fold_predictions_first_repeat")
    cohort = read("clinical_baselines").set_index("pid")
    sl = read("predicted_and_reference_slopes_first_repeat").query("followup >= 2 and spaced_dates >= 3")
    overall = read("prediction_metrics")
    first = overall.query("repetition == 1 and fold == 0 and weighting == 'patient'")
    classification = read("rapid_decline_metrics").query(
        "repetition == 1 and fold == 0 and minimum_reference_followup == 2 and reference == 'reference_lmm'")
    hz = read("horizon_metrics").query("repetition == 1 and model in ['persistence', 'M3'] and weighting == 'patient'")

    fig, axes = plt.subplots(1, 3, figsize=(8.8, 3.05), layout="constrained")
    comparisons = [
        (pred.egfr, pred.M3, "eGFR"),
        (pred.egfr - pred.baseline_egfr, pred.M3 - pred.baseline_egfr, "eGFR change"),
        (sl.reference_lmm, sl.M3, "annual eGFR slope"),
    ]
    for ax, (y, p, label) in zip(axes, comparisons):
        if len(y) > 1000:
            ax.hexbin(y, p, gridsize=50, cmap="Blues", bins="log", mincnt=1)
        else:
            ax.scatter(y, p, color=BLUE, s=9, alpha=.62, linewidths=0)
        lim = [min(y.min(), p.min()), max(y.max(), p.max())]
        ax.plot(lim, lim, "--", color=GRAY, lw=.8)
        ax.set(xlabel=f"Observed {label}", ylabel=f"Predicted {label}", xlim=lim, ylim=lim)
    export(fig, "observed_vs_predicted")

    fig, ax = plt.subplots(figsize=(4, 3), layout="constrained")
    diff = sl.M3 - sl.reference_lmm
    ax.scatter((sl.M3 + sl.reference_lmm) / 2, diff, color=BLUE, alpha=.65, s=11)
    ax.axhline(diff.mean(), color=GRAY)
    for bound in [diff.mean() - 1.96 * diff.std(), diff.mean() + 1.96 * diff.std()]:
        ax.axhline(bound, linestyle="--", color=GRAY)
    ax.set(xlabel="Mean predicted and reference slopes", ylabel="Predicted minus reference slope")
    export(fig, "slope_agreement")

    fig, axes = plt.subplots(2, 3, figsize=(8, 4.6), sharey=True, layout="constrained")
    for i, kind in enumerate(("native", "transplant")):
        d = cohort[cohort.kind.eq(kind)].sort_values("followup")
        for j, fraction in enumerate((.25, .5, .75)):
            ax = axes[i, j]
            pid = d.index[int(fraction * len(d))]
            base = cohort.loc[pid]
            person = pred[pred.pid.eq(pid)].sort_values("t")
            line = person.M3 - sl.M3[sl.pid.eq(pid)].iloc[0] * person.t if pid in set(sl.pid) else None
            ax.scatter(person.t, person.egfr, color=BLUE, s=9, label="Follow-up eGFR")
            if line is None:
                intercept = np.polyfit(person.t, person.M3, 1)[1]
            else:
                intercept = line.iloc[0]
            ax.plot(np.r_[0, person.t], np.r_[intercept, person.M3], color=ROSE, label="Held-out prediction")
            base_t = (pd.Timestamp(base.baseline_egfr_date) - pd.Timestamp(base.index_date)).days / 365.25
            ax.scatter(base_t, base.baseline_egfr, marker="D", color=GRAY, s=14, label="Baseline eGFR")
            ax.set_title(f"{kind.capitalize()}  /  follow-up percentile {int(fraction * 100)}", fontsize=8)
            if i == 1:
                ax.set_xlabel("Years from baseline")
            if j == 0:
                ax.set_ylabel("eGFR")
    axes[0, 0].legend(frameon=False, loc="best")
    export(fig, "example_trajectories")

    fig, axes = plt.subplots(1, 2, figsize=(7, 2.7), layout="constrained")
    for ax, mode in zip(axes, ("egfr", "change")):
        for m, color in [("M3", BLUE), ("persistence", GRAY)]:
            d = hz[hz.outcome.eq(mode) & hz.model.eq(m)]
            ax.plot(d.horizon, d.rmse, "o-", color=color, ms=3, label=m.capitalize())
        d = hz[hz.outcome.eq(mode) & hz.model.eq("M3")]
        ax.set(xticks=d.horizon, xticklabels=[f"{x:g}\n{n}" for x, n in zip(d.horizon, d.patients)],
               xlabel="Years from baseline and patient count",
               ylabel="eGFR RMSE" if mode == "egfr" else "eGFR change RMSE")
    axes[0].legend(frameon=False)
    export(fig, "prediction_by_horizon")

    fig, axes = plt.subplots(2, 2, figsize=(6, 5.25), layout="constrained")
    for i, cutoff in enumerate((-3, -5)):
        outcome = sl.reference_lmm.lt(cutoff)
        observed_rate = outcome.mean()
        for j, (model, color) in enumerate(zip(["M1", "M2", "M3"], [GRAY, BLUE, ROSE])):
            score = -sl[model]
            fpr, tpr, _ = roc_curve(outcome, score)
            precision, recall, _ = precision_recall_curve(outcome, score)
            d = classification[(classification.cutoff == cutoff) & classification.model.eq(model)
                               & classification.kidney.eq("all")].iloc[0]
            axes[0, i].plot(fpr, tpr, color=color, label=f"{model}  {d.roc_auc:.2f}")
            axes[1, i].plot(recall, precision, color=color, label=f"{model}  {d.pr_auc_ap:.2f}")
        axes[0, i].plot([0, 1], [0, 1], "--", color=GRAY, lw=.6)
        axes[1, i].axhline(observed_rate, linestyle="--", color=GRAY, lw=.6)
        for ax in axes[:, i]:
            ax.legend(frameon=False)
            ax.set_xlim(0, 1)
            ax.set_ylim(0, 1.02)
        axes[0, i].set_title(f"Reference slope < {cutoff}")
        axes[0, i].set(xlabel="False positive rate", ylabel="Sensitivity")
        axes[1, i].set(xlabel="Recall", ylabel="Precision")
    export(fig, "rapid_decline_discrimination")

    fig, axes = plt.subplots(1, 2, figsize=(7, 2.6), layout="constrained")
    terms = read("fixed_effects").query("model == 'M3_full_cohort'").set_index("term")
    pretty = {"baseline_egfr": "Baseline eGFR", "age": "Age", "female": "Female",
              "log_acr": "Log(1 + ACR)"}
    for ax, prefix in zip(axes, ("", "t:")):
        d = terms.loc[[prefix + c for c in pretty]]
        ax.errorbar(d.estimate, range(len(d)),
                    xerr=[d.estimate - d.ci_lower, d.ci_upper - d.estimate], fmt="o", ms=4, color=BLUE)
        ax.axvline(0, ls="--", color=GRAY, lw=.7)
        ax.set(yticks=range(len(d)), yticklabels=list(pretty.values()),
               xlabel="Main-effect coefficient" if not prefix else "Time-interaction coefficient")
        ax.invert_yaxis()
    export(fig, "clinical_fixed_effects")

    fig, axes = plt.subplots(1, 3, figsize=(8.3, 2.8), layout="constrained")
    values = first.query("kidney == 'all' and outcome == 'egfr'").set_index("model")
    axes[0].barh(COMPARATORS, values.loc[COMPARATORS].rmse, color=[GRAY, BLUE, BLUE, BLUE, ROSE])
    axes[0].invert_yaxis()
    axes[0].set_xlabel("Patient-weighted eGFR RMSE")
    minimum_models = ["M3" if t == 0 else f"M3_min_{t:g}y" for t in THRESHOLDS]
    axes[1].plot(THRESHOLDS, values.loc[minimum_models].rmse, "o-", color=ROSE, ms=4)
    axes[1].set(xlabel="Minimum training follow-up (years)", ylabel="Patient-weighted eGFR RMSE", xticks=THRESHOLDS)
    sub = first.query("kidney != 'all' and outcome == 'egfr'").pivot(index="kidney", columns="model", values="rmse")
    x = np.arange(len(sub))
    axes[2].bar(x - .18, sub.persistence, .36, color=GRAY, label="Persistence")
    axes[2].bar(x + .18, sub.M3, .36, color=ROSE, label="M3")
    axes[2].set(xticks=x, xticklabels=[s.capitalize() for s in sub.index], ylabel="Patient-weighted eGFR RMSE")
    axes[2].legend(frameon=False)
    export(fig, "model_cohort_and_kidney_comparisons")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    render(parser.parse_args().output)

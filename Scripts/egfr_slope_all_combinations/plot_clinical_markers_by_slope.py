"""Plot baseline ACR, eGFR and 2-year KFRE in the Model 3 radiomics cohort."""

from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import FuncFormatter


ROUGH = Path(__file__).resolve().parents[2] / "Results_Rough"
OUTPUT = ROUGH / "egfr_slope_all_combinations"
SLOPES = OUTPUT / "inputs" / "patient_features_6mod.csv"
CLINICAL = (
    ROUGH / "01_eGFR_slope_model3_min1yr" / "Heat_map" / "model_cv"
    / "clinical_addon" / "baseline_clinical_model3_91.csv"
)
MARKERS = (
    ("baseline_ACR_mg_per_mmol", "Baseline urine ACR (mg/mmol)"),
    ("baseline_eGFR_ckdepi2021", "Baseline eGFR (mL/min/1.73 m²)"),
    ("KFRE_2yr_risk", "KFRE 2-year risk (%)"),
)
BINARY_COLORS = ("#4C9AC6", "#E84A75")
ORDINAL_COLORS = ("#4C9AC6", "#F1B06D", "#F08078", "#E84A75")
INK = "#242424"


def load_data():
    slopes = pd.read_csv(SLOPES, usecols=["patient_id", "slope"])
    clinical = pd.read_csv(CLINICAL, usecols=["patient_id", *(col for col, _ in MARKERS)])
    data = slopes.merge(clinical, on="patient_id", how="left", validate="one_to_one", indicator=True)
    if not data["_merge"].eq("both").all() or not np.isfinite(data["slope"]).all():
        raise ValueError("Some cohort patients lack a matched clinical record or a finite slope.")
    data["KFRE_2yr_risk"] *= 100
    return data


def draw_figure(data, categories, labels, colors, heading):
    fig, axes = plt.subplots(1, 3, figsize=(9.2, 3.45), layout="constrained")
    fig.suptitle(f"{heading}  |  Model 3 radiomics cohort (N = {len(data)})", fontsize=10)
    fig.supxlabel("LMM eGFR slope (mL/min/1.73 m²/year)", fontsize=8)
    rng = np.random.default_rng(7)

    for ax, (column, ylabel) in zip(axes, MARKERS):
        values = data[column].where(np.isfinite(data[column]))
        if column == "baseline_ACR_mg_per_mmol":
            values = values.where(values.ge(0))
            ax.set_yscale("symlog", linthresh=0.1)
            ax.set_yticks([0, 0.1, 1, 10, 100])
            ylabel += "\n(linear to 0.1, log above)"
        elif column == "KFRE_2yr_risk":
            values = values.where(values.gt(0))
            ax.set_yscale("log")
            ylabel += "\n(log scale)"
        ax.yaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value:g}"))

        groups = [values[categories == i].dropna().to_numpy() for i in range(len(labels))]
        for i, (group, color) in enumerate(zip(groups, colors)):
            if not group.size:
                continue
            box = ax.boxplot(
                [group], positions=[i], widths=0.54, whis=1.5, showfliers=False, patch_artist=True,
                boxprops={"edgecolor": INK, "linewidth": 0.7},
                medianprops={"color": INK, "linewidth": 1.2},
                whiskerprops={"color": INK, "linewidth": 0.7},
                capprops={"color": INK, "linewidth": 0.7},
            )
            box["boxes"][0].set_facecolor(matplotlib.colors.to_rgba(color, 0.24))
            ax.scatter(
                i + rng.uniform(-0.19, 0.19, group.size), group, s=13,
                color=color, alpha=0.85, edgecolors="white", linewidths=0.35, zorder=3,
            )
            box["medians"][0].set_zorder(4)

        ax.set_xticks(range(len(labels)), [f"{label}\n(n = {len(group)})" for label, group in zip(labels, groups)])
        ax.set_xlim(-0.55, len(labels) - 0.45)
        ax.set_ylabel(ylabel)
        ax.spines[["top", "right"]].set_visible(False)
        ax.spines[["left", "bottom"]].set_linewidth(0.6)
        ax.tick_params(direction="out", width=0.6, length=3)
        ax.margins(y=0.07)

    return fig


def main():
    plt.rcParams.update({
        "font.family": ["Arial", "DejaVu Sans"],
        "font.size": 7,
        "axes.labelsize": 8,
        "pdf.fonttype": 42,
    })
    data = load_data()
    slope = data["slope"]
    schemes = [
        (
            f"binary_lt{cut}", (slope < -cut).astype(int),
            (f"≥ −{cut}", f"< −{cut}"), BINARY_COLORS, f"Binary cutoff < −{cut}",
        )
        for cut in (1, 3, 5)
    ]
    schemes.append((
        "ordinal4", sum((slope < cut).astype(int) for cut in (-1, -3, -5)),
        ("≥ −1", "[−3, −1)", "[−5, −3)", "< −5"), ORDINAL_COLORS, "Four ordinal groups",
    ))

    OUTPUT.mkdir(parents=True, exist_ok=True)
    stem = "acr_egfr_kfre_by_slope"
    with PdfPages(OUTPUT / f"{stem}.pdf") as pdf:
        for name, categories, labels, colors, heading in schemes:
            fig = draw_figure(data, categories, labels, colors, heading)
            fig.savefig(OUTPUT / f"{stem}_{name}.png", dpi=600)
            pdf.savefig(fig)
            plt.close(fig)
            print(name, categories.value_counts().sort_index().to_list())


if __name__ == "__main__":
    main()

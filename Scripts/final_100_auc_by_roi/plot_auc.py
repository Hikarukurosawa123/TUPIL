"""Plot the existing Final_100 AUC and DeLong results by ROI and time window."""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch

RESULTS = Path(__file__).resolve().parents[2] / "Results_Rough" / "Final_100"
MODELS = {
    "frozen": "Frozen source model",
    "polr": "Ordinal logistic retrained within each half",
}
ROIS = {
    "all_images": "All images",
    "matched": "Matched third",
    "opposite": "Opposite third",
}
HALVES = {
    "early": ("Early (0–43 d, n = 50)", "#488FB4", -0.18),
    "late": ("Late (45–144 d, n = 50)", "#B4D5E5", 0.18),
}
ROI_COMPARISONS = (
    "All images vs Matched third",
    "All images vs Opposite third",

)


def p_label(p):
    return f"{p:.3f}" + ("*" if p < 0.05 else "")


def main(
    workbook=RESULTS / "spatial_Final_100.xlsx",
    output_dir=RESULTS / "auc_by_roi_and_interval",
    cv_label=None,
    include_binary=False,
):
    with pd.ExcelFile(workbook) as source:
        auc = pd.read_excel(source, "AUC").set_index(["model", "half", "spatial"]).sort_index()
        timing = pd.read_excel(source, "Early vs late DeLong").set_index(["model", "spatial"])
        paired = pd.read_excel(source, "ROI comparisons DeLong").set_index(
            ["model", "half", "comparison"]
        )

    plt.rcParams.update({"font.family": "DejaVu Sans", "pdf.fonttype": 42})
    models = MODELS.copy()
    if include_binary:
        models["logistic"] = "Binary logistic retrained within each half"
    panel_count = len(models)
    fig = plt.figure(figsize=(6.2 * panel_count, 7.8), facecolor="white")
    grid = fig.add_gridspec(
        2, panel_count, left=0.072, right=0.975, bottom=0.13, top=0.78,
        height_ratios=[3.4, 1], wspace=0.15, hspace=0.31,
    )
    fig.suptitle("AUC by ROI criterion and scan-to-biopsy interval", y=0.962,
                 fontsize=19, fontweight="bold")
    subtitle = "100 cases · Kidney minus pelvis · CI 0 vs CI 1–2"
    if cv_label:
        subtitle += f" · {cv_label}"
    fig.text(0.5, 0.905, subtitle,
             ha="center", fontsize=11, color="#394754")
    fig.legend(
        handles=[Patch(facecolor=color, edgecolor="#375063", label=label)
                 for label, color, _ in HALVES.values()],
        loc="upper center", bbox_to_anchor=(0.5, 0.876), ncol=2,
        frameon=False, fontsize=10, columnspacing=2.6,
    )
    positions = np.arange(len(ROIS))
    for col, (model, title) in enumerate(models.items()):
        ax = fig.add_subplot(grid[0, col])
        ax.set_title(title, fontsize=12, pad=12)
        ax.spines[["top", "right"]].set_visible(False)
        ax.spines[["bottom", "left"]].set_color("#65717C")
        ax.set_axisbelow(True)
        ax.axhline(0.5, color="#65717C", linestyle=":", linewidth=1)
        ax.set(xlim=(-0.59, 2.59), ylim=(0, 1.40), yticks=np.arange(0, 1.01, 0.2))
        ax.set_xticks(positions, list(ROIS.values()), fontsize=10)
        ax.tick_params(axis="x", length=0, pad=8)
        ax.tick_params(axis="y", labelsize=9, colors="#394754")
        if col == 0:
            ax.set_ylabel("AUC with DeLong 95% CI", fontsize=10)

        for half, (_, color, offset) in HALVES.items():
            rows = auc.loc[(model, half)].loc[list(ROIS)]
            values = rows["AUC"].to_numpy()
            low, high = rows["CI_low"].to_numpy(), rows["CI_high"].to_numpy()
            ax.bar(
                positions + offset, values, width=0.31, color=color,
                edgecolor="#375063", linewidth=0.75,
                yerr=np.array([values - low, high - values]), capsize=3.5,
                error_kw={"ecolor": "#28343D", "elinewidth": 1, "capthick": 1},
            )
            for x, value, upper in zip(positions + offset, values, high):
                ax.text(x, upper + 0.014, f"{value:.2f}", ha="center",
                        va="bottom", fontsize=9)

        for x, roi in zip(positions, ROIS):
            test = timing.loc[(model, roi)]
            np.testing.assert_allclose(
                [test["AUC_early"], test["AUC_late"]],
                [auc.loc[(model, half, roi), "AUC"] for half in HALVES],
            )
            left, right = (x + offset for _, _, offset in HALVES.values())
            ax.plot([left, left, right, right], [1.055, 1.077, 1.077, 1.055],
                    color="#28343D", linewidth=0.85)
            ax.text(x, 1.09, f"p = {p_label(test['p_DeLong_1sided'])}",
                    ha="center", va="bottom", fontsize=9.2)

        for (half, (_, _, offset)), height, color in zip(
            HALVES.items(), (1.18, 1.29), ("#23556F", "#415662"),
        ):
            p = paired.loc[
                (model, half, "Matched third vs Opposite third"), "p_DeLong_1sided"
            ]
            left, right = positions[1:] + offset
            linestyle = "-" if half == "early" else "--"
            ax.plot([left, left, right, right],
                    [height - 0.025, height, height, height - 0.025],
                    color=color, linewidth=1, linestyle=linestyle)
            ax.text((left + right) / 2, height + 0.012,
                    f"{half.capitalize()} matched > opposite, p = {p_label(p)}",
                    ha="center", va="bottom", color=color, fontsize=8.8)

        table_ax = fig.add_subplot(grid[1, col])
        table_ax.set_axis_off()
        table_ax.text(0, 1.01, "All-images comparisons within each time window",
                      fontsize=10, va="bottom", transform=table_ax.transAxes)
        table_rows = [
            [comparison, *[p_label(paired.loc[(model, half, comparison), "p_DeLong_1sided"])
                           for half in HALVES]]
            for comparison in ROI_COMPARISONS
        ]
        table = table_ax.table(
            cellText=table_rows, colLabels=["Comparison", "Early p", "Late p"],
            colWidths=[0.60, 0.20, 0.20], cellLoc="center", bbox=(0, 0, 1, 0.94),
        )
        table.auto_set_font_size(False)
        table.set_fontsize(9)
        for (row, column), cell in table.get_celld().items():
            cell.set_facecolor("#EFF3F6" if row == 0 else "white")
            cell.set_edgecolor("#D4DCE3")
            cell.set_linewidth(0.55)
            if column == 0:
                cell.set_text_props(ha="left")
                cell.PAD = 0.035

    fig.text(
        0.072, 0.079,
        "Early > late within each ROI uses unpaired DeLong tests. Matched > opposite and the table comparisons use paired DeLong tests.",
        fontsize=9, color="#394754",
    )
    fig.text(
        0.072, 0.051,
        "Table p-values test the first ROI > the second. "
        "All p-values are one-sided (* p < 0.05). No multiple-comparison correction was applied.",
        fontsize=9, color="#394754",
    )


    output_dir.mkdir(parents=True, exist_ok=True)
    for extension in ("png", "pdf"):
        fig.savefig(output_dir / f"AUC_by_ROI_and_interval.{extension}", dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    main()

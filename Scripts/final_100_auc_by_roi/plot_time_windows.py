"""Render the reference AUC layout for the Final_100 month-based windows."""

from itertools import combinations

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch

from plot_auc import RESULTS, ROIS, ROI_COMPARISONS, p_label


COHORT = RESULTS / "ci_spatial_0_1_3_5_months"
WORKBOOK = COHORT / "spatial_time_windows_Final_100.xlsx"
OUTPUT = COHORT / "auc_by_roi_and_interval"
MODELS = {
    "frozen": "Frozen source model",
    "polr": "Ordinal logistic retrained within each window",
}
WINDOWS = {
    "interval": {
        "day0": ("Same day", "#213E56"),
        ">0-1m": (">0-1 mo", "#488FB4"),
        ">1-3m": (">1-3 mo", "#82B5CF"),
        ">3-5m": (">3-5 mo", "#C2DCE8"),
    },
    "cumulative": {
        "0-1m": ("0-1 mo", "#213E56"),
        "0-3m": ("0-3 mo", "#488FB4"),
        "0-5m": ("0-5 mo", "#C2DCE8"),
    },
}
PAIRED_ROIS = (
    ("all_images", "matched"),
    ("all_images", "opposite"),
    ("matched", "opposite"),
)
COMPARISONS = (*ROI_COMPARISONS, "Matched third vs Opposite third")
P_VALUE = "p_DeLong_1sided"


def draw_table(ax, title, headers, rows, widths):
    ax.set_axis_off()
    ax.set_title(title, loc="left", fontsize=10.5, pad=4)
    table = ax.table(
        cellText=rows, colLabels=headers, colWidths=widths,
        cellLoc="center", bbox=(0, 0, 1, 1),
    )
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    for (row, col), cell in table.get_celld().items():
        cell.set_facecolor("#EFF3F6" if row == 0 else "white")
        cell.set_edgecolor("#D4DCE3")
        cell.set_linewidth(0.55)
        if col == 0:
            cell.set_text_props(ha="left")
            cell.PAD = 0.035


def draw_bars(ax, title, windows, aucs):
    positions = np.arange(len(ROIS))
    step = 0.82 / len(windows)
    offsets = (np.arange(len(windows)) - (len(windows) - 1) / 2) * step
    ax.set_title(title, fontsize=12, pad=10)
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["bottom", "left"]].set_color("#65717C")
    ax.set_axisbelow(True)
    ax.axhline(0.5, color="#65717C", linestyle=":", linewidth=1)
    ax.set(xlim=(-0.55, 2.55), ylim=(0, 1.10), yticks=np.arange(0, 1.01, 0.2))
    ax.set_xticks(positions, list(ROIS.values()), fontsize=10)
    ax.tick_params(axis="x", length=0, pad=8)
    ax.tick_params(axis="y", labelsize=9, colors="#394754")

    for offset, (window, (_, color)) in zip(offsets, windows.items()):
        rows = aucs.loc[window].loc[list(ROIS)]
        values = rows["AUC"].to_numpy()
        low = rows["CI_low"].to_numpy()
        high = rows["CI_high"].to_numpy()
        ax.bar(
            positions + offset, values, width=step * 0.87, color=color,
            edgecolor="#375063", linewidth=0.75,
            yerr=[values - low, high - values], capsize=3.5,
            error_kw={"ecolor": "#28343D", "elinewidth": 1, "capthick": 1},
        )
        for x, value, upper in zip(positions + offset, values, high):
            ax.text(x, upper + 0.015, f"{value:.2f}", ha="center",
                    va="bottom", fontsize=8.5)


def paired_table(paired, aucs, windows):
    rows = []
    for comparison, rois in zip(COMPARISONS, PAIRED_ROIS):
        row = [comparison]
        for window in windows:
            test = paired.loc[(window, comparison)]
            np.testing.assert_allclose(
                test[["AUC_1", "AUC_2"]].to_numpy(float),
                aucs.loc[window].loc[list(rois), "AUC"],
                atol=1e-12, rtol=0,
            )
            row.append(p_label(test[P_VALUE]))
        rows.append(row)
    return rows


def interval_table(timing, aucs, windows):
    rows = []
    for earlier, later in combinations(windows, 2):
        row = [f"{windows[earlier][0]} vs {windows[later][0]}"]
        for roi in ROIS:
            test = timing.loc[(earlier, later, roi)]
            np.testing.assert_allclose(
                test[["AUC_earlier", "AUC_later"]].to_numpy(float),
                [aucs.loc[(earlier, roi), "AUC"], aucs.loc[(later, roi), "AUC"]],
                atol=1e-12, rtol=0,
            )
            row.append(p_label(test[P_VALUE]))
        rows.append(row)
    return rows


def plot(window_set, metadata, aucs, paired, timing, output):
    windows = WINDOWS[window_set]
    disjoint = window_set == "interval"
    figure = plt.figure(figsize=(13.6, 10.2 if disjoint else 7.6), facecolor="white")
    grid = figure.add_gridspec(
        3 if disjoint else 2, len(MODELS), left=0.065, right=0.975,
        bottom=0.13 if disjoint else 0.16, top=0.805 if disjoint else 0.765,
        height_ratios=[3.5, 1.0, 1.6] if disjoint else [3.7, 1.05],
        wspace=0.15, hspace=0.4 if disjoint else 0.38,
    )
    title_window = "interval" if disjoint else "cumulative window"
    figure.suptitle(f"AUC by ROI criterion and scan-to-biopsy {title_window}",
                    y=0.959, fontsize=19, fontweight="bold")
    figure.text(0.5, 0.901 if disjoint else 0.890,
                "Final_100 · Kidney minus pelvis · CI 0 vs CI 1-2 · 1 month = 30 days",
                ha="center", fontsize=11, color="#394754")
    figure.legend(
        handles=[
            Patch(
                facecolor=color, edgecolor="#375063",
                label=f"{label}\n(n = {metadata.loc[key, 'n']}, "
                      f"CI 1-2 = {metadata.loc[key, 'n_ci12']})",
            )
            for key, (label, color) in windows.items()
        ],
        loc="upper center", bbox_to_anchor=(0.5, 0.885), ncol=len(windows),
        frameon=False, fontsize=10, columnspacing=2.6,
    )

    for col, (model, title) in enumerate(MODELS.items()):
        model_aucs = aucs.loc[model]
        ax = figure.add_subplot(grid[0, col])
        draw_bars(ax, title, windows, model_aucs)
        if col == 0:
            ax.set_ylabel("AUC with DeLong 95% CI", fontsize=10)
        draw_table(
            figure.add_subplot(grid[1, col]),
            "ROI comparisons within each window (paired DeLong p)",
            ["Comparison", *[label for label, _ in windows.values()]],
            paired_table(paired.loc[model], model_aucs, windows),
            [0.43, *[0.57 / len(windows)] * len(windows)],
        )
        if disjoint:
            draw_table(
                figure.add_subplot(grid[2, col]),
                "Between-interval comparisons by ROI (unpaired DeLong p)",
                ["Earlier vs later", *ROIS.values()],
                interval_table(timing.loc[model], model_aucs, windows),
                [0.40, 0.20, 0.20, 0.20],
            )

    notes = (
        [
            "Same day is supplementary. Remaining intervals cover days 1-30, 31-90, and 91-150.",
            "One-sided tests compare the first ROI with the second (first > second) and earlier with later (earlier > later).",
        ]
        if disjoint else
        [
            "Cumulative windows start on day 0 and end on days 30, 90, and 150. Between-window tests are omitted because the samples are nested.",
            "One-sided tests compare the first ROI with the second (first > second).",
        ]
    )
    notes.append("* p < 0.05. P-values are uncorrected for multiple comparisons.")
    for y, note in zip((0.087, 0.063, 0.039), notes):
        figure.text(0.065, y, note, fontsize=9, color="#394754")

    stem = "AUC_by_ROI_and_interval" if disjoint else "AUC_by_ROI_and_cumulative_window"
    for extension in ("png", "pdf"):
        figure.savefig(output / f"{stem}.{extension}", dpi=300)
    plt.close(figure)


def main(workbook=WORKBOOK, output=OUTPUT):
    with pd.ExcelFile(workbook) as source:
        metadata = pd.read_excel(source, "Windows").set_index("window")
        aucs = pd.read_excel(source, "AUC").set_index(
            ["model", "window", "spatial"],
        ).sort_index()
        paired = pd.read_excel(source, "ROI comparisons DeLong").set_index(
            ["model", "window", "comparison"],
        ).sort_index()
        timing = pd.read_excel(source, "Interval vs interval DeLong").set_index(
            ["model", "window_earlier", "window_later", "spatial"],
        ).sort_index()
    if any(not table.index.is_unique for table in (metadata, aucs, paired, timing)):
        raise ValueError("The cohort workbook contains duplicate result keys.")
    for key, row in aucs.iterrows():
        np.testing.assert_equal(
            row[["n", "n_ci12"]].to_numpy(int),
            metadata.loc[key[1], ["n", "n_ci12"]].to_numpy(int),
        )

    plt.rcParams.update({"font.family": "DejaVu Sans", "pdf.fonttype": 42})
    output.mkdir(parents=True, exist_ok=True)
    for window_set in WINDOWS:
        plot(window_set, metadata, aucs, paired, timing, output)


if __name__ == "__main__":
    main()

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from lifelines import KaplanMeierFitter
from matplotlib.gridspec import GridSpec


SUBTYPE_ORDER = ["MSI", "EMT", "MSS/TP53-", "MSS/TP53+"]
SUBTYPE_COLORS = {
    "MSI": "#174c7e",
    "EMT": "#c06c2b",
    "MSS/TP53-": "#1f8a70",
    "MSS/TP53+": "#7b5ea7",
}
MARKER_GROUP_COLORS = {
    "EMT_high": "#c06c2b",
    "epithelial_high": "#174c7e",
    "proliferation_context": "#1f8a70",
}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def draw_km_panel(ax: plt.Axes, km_rows: list[dict[str, str]]) -> None:
    frame = pd.DataFrame(km_rows)
    frame["overall_survival_months"] = pd.to_numeric(frame["overall_survival_months"])
    frame["overall_survival_event"] = pd.to_numeric(frame["overall_survival_event"])

    km_fitter = KaplanMeierFitter()
    for subtype in SUBTYPE_ORDER:
        subframe = frame[frame["molecular_subtype_label"] == subtype]
        if subframe.empty:
            continue
        km_fitter.fit(
            durations=subframe["overall_survival_months"],
            event_observed=subframe["overall_survival_event"],
            label=f"{subtype} (n={len(subframe)})",
        )
        km_fitter.plot(ax=ax, ci_show=False, color=SUBTYPE_COLORS[subtype], linewidth=2.2)

    ax.set_title("A  Derivation-cohort overall survival", loc="left", fontsize=14, fontweight="bold")
    ax.set_xlabel("Months")
    ax.set_ylabel("Survival probability")
    ax.set_facecolor("#fbfaf7")
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False, loc="lower left")



def draw_risk_table(ax: plt.Axes, km_rows: list[dict[str, str]]) -> None:
    frame = pd.DataFrame(km_rows)
    frame["overall_survival_months"] = pd.to_numeric(frame["overall_survival_months"])
    timepoints = [0, 24, 48, 72]
    headers = ["Subtype", *[str(timepoint) for timepoint in timepoints]]
    table_rows = []
    for subtype in SUBTYPE_ORDER:
        subframe = frame[frame["molecular_subtype_label"] == subtype]
        if subframe.empty:
            continue
        counts = [str(int((subframe["overall_survival_months"] >= timepoint).sum())) for timepoint in timepoints]
        table_rows.append([subtype, *counts])

    ax.axis("off")
    ax.text(0.0, 1.08, "Number at risk", fontsize=11.5, fontweight="bold", color="#222222", transform=ax.transAxes)
    table = ax.table(cellText=table_rows, colLabels=headers, loc="center", cellLoc="center", colLoc="center")
    table.auto_set_font_size(False)
    table.set_fontsize(9.5)
    table.scale(1.0, 1.28)
    for (row_index, col_index), cell in table.get_celld().items():
        cell.set_edgecolor("#fbfaf7")
        cell.set_linewidth(0.0)
        if row_index == 0:
            cell.set_facecolor("#ece7dd")
            cell.set_text_props(weight="bold")
        else:
            cell.set_facecolor("#f4efe5" if row_index % 2 == 1 else "#fbfaf7")
            if col_index == 0:
                subtype = table_rows[row_index - 1][0]
                cell.set_text_props(color=SUBTYPE_COLORS[subtype], weight="bold")


def draw_marker_panel(ax: plt.Axes, marker_rows: list[dict[str, str]]) -> None:
    frame = pd.DataFrame(marker_rows)
    subtype_categories = [subtype for subtype in SUBTYPE_ORDER if subtype in frame["subtype"].unique()]
    pivot = (
        frame.assign(marker_label=frame["gene_symbol"])
        .pivot(index="marker_label", columns="subtype", values="subtype_mean_zscore")
        .reindex(columns=subtype_categories)
    )
    pivot = pivot.apply(pd.to_numeric)
    image = ax.imshow(pivot.values, cmap="RdBu_r", aspect="auto", vmin=-1.6, vmax=1.6)
    ax.set_title("B  Expression-level state portrait", loc="left", fontsize=14, fontweight="bold")
    ax.set_xticks(range(len(pivot.columns)), pivot.columns, rotation=20)
    ax.set_yticks(range(len(pivot.index)), pivot.index)
    ax.set_facecolor("#fbfaf7")
    ax.set_xlim(-1.85, len(pivot.columns) - 0.5)
    for spine in ax.spines.values():
        spine.set_visible(False)

    marker_groups = frame[["gene_symbol", "marker_group", "pathway_label", "marker_order"]].drop_duplicates().sort_values("marker_order").reset_index(drop=True)
    for row_index, marker_group in enumerate(marker_groups["marker_group"]):
        ax.add_patch(
            plt.Rectangle(
                (-1.72, row_index - 0.5),
                0.22,
                1.0,
                facecolor=MARKER_GROUP_COLORS.get(marker_group, "#9b9b9b"),
                edgecolor="#fbfaf7",
                linewidth=0.0,
            )
        )

    for row_index in range(pivot.shape[0]):
        for col_index in range(pivot.shape[1]):
            ax.text(col_index, row_index, f"{pivot.iloc[row_index, col_index]:.2f}", ha="center", va="center", fontsize=8.5)

    plt.colorbar(image, ax=ax, fraction=0.04, pad=0.03, label="Subtype-mean z score")


def render_figure2(km_rows: list[dict[str, str]], marker_rows: list[dict[str, str]], output_pdf: Path, output_png: Path) -> None:
    output_pdf.parent.mkdir(parents=True, exist_ok=True)
    output_png.parent.mkdir(parents=True, exist_ok=True)

    figure = plt.figure(figsize=(13.6, 8.6))
    figure.patch.set_facecolor("#fbfaf7")
    grid = GridSpec(2, 2, figure=figure, height_ratios=[1.0, 0.28], width_ratios=[1.1, 1.05], hspace=0.18, wspace=0.3)
    km_axis = figure.add_subplot(grid[0, 0])
    risk_axis = figure.add_subplot(grid[1, 0])
    marker_axis = figure.add_subplot(grid[:, 1])

    draw_km_panel(km_axis, km_rows)
    draw_risk_table(risk_axis, km_rows)
    draw_marker_panel(marker_axis, marker_rows)

    figure.suptitle(
        "Discovery-cohort prognosis and expression portrait",
        x=0.5,
        y=0.98,
        ha="center",
        fontsize=17,
        fontweight="bold",
        color="#111111",
    )
    figure.text(
        0.5,
        0.94,
        "The derivation cohort shows four-state survival structure together with an interpretable marker-level portrait.",
        ha="center",
        fontsize=11,
        color="#4a4a4a",
    )
    figure.savefig(output_pdf, bbox_inches="tight")
    figure.savefig(output_png, dpi=300, bbox_inches="tight")
    plt.close(figure)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Render submission Figure 2 derivation KM and marker portrait")
    parser.add_argument("--km-input", required=True)
    parser.add_argument("--marker-input", required=True)
    parser.add_argument("--pdf-output", required=True)
    parser.add_argument("--png-output", required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    render_figure2(
        read_tsv(Path(args.km_input)),
        read_tsv(Path(args.marker_input)),
        Path(args.pdf_output),
        Path(args.png_output),
    )


if __name__ == "__main__":
    main()
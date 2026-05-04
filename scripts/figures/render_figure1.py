from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.gridspec import GridSpec


STEP_TITLES = {
    "derivation": "Derivation",
    "replication": "Replication",
    "orthogonal_validation": "Orthogonal Validation",
}
STEP_NOTES = {
    "derivation": "Primary dataset",
    "replication": "External validation",
    "orthogonal_validation": "Clinicopathologic context",
}

ROLE_COLORS = {
    "derivation": "#174c7e",
    "replication": "#1f8a70",
    "orthogonal": "#c06c2b",
    "sensitivity": "#7b5ea7",
    "treatment": "#9a3b3b",
}

SUBTYPE_ORDER = ["EMT", "MSI", "MSS/TP53-", "MSS/TP53+"]
SUBTYPE_COLORS = {
    "EMT": "#c06c2b",
    "MSI": "#174c7e",
    "MSS/TP53-": "#1f8a70",
    "MSS/TP53+": "#7b5ea7",
}
CONFIDENCE_ALPHA = {
    "high": 0.92,
    "medium": 0.58,
    "low": 0.28,
}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def draw_workflow(ax: plt.Axes, flow_rows: list[dict[str, str]]) -> None:
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    x_positions = [0.18, 0.5, 0.82]
    box_width = 0.22
    box_height = 0.32
    for x_position, row in zip(x_positions, flow_rows, strict=True):
        role = row["step"]
        color = ROLE_COLORS.get(role, "#4a4a4a")
        rectangle = plt.Rectangle(
            (x_position - box_width / 2, 0.42 - box_height / 2),
            box_width,
            box_height,
            facecolor="#f7f4ef",
            edgecolor=color,
            linewidth=2.2,
            joinstyle="round",
        )
        ax.add_patch(rectangle)
        ax.text(
            x_position,
            0.545,
            STEP_TITLES.get(role, role.replace("_", " ").title()),
            ha="center",
            va="center",
            fontsize=13,
            fontweight="bold",
            color="#1f1f1f",
        )
        ax.text(
            x_position,
            0.485,
            row["dataset_id"],
            ha="center",
            va="center",
            fontsize=12,
            color=color,
        )
        ax.text(
            x_position,
            0.355,
            STEP_NOTES.get(role, row["note"]),
            ha="center",
            va="center",
            fontsize=9.2,
            color="#4a4a4a",
            wrap=True,
        )

    for left_x, right_x in zip(x_positions, x_positions[1:], strict=False):
        if right_x is None:
            continue
        ax.annotate(
            "",
            xy=(right_x - box_width / 2 - 0.02, 0.42),
            xytext=(left_x + box_width / 2 + 0.02, 0.42),
            arrowprops=dict(arrowstyle="->", lw=2, color="#6c6c6c"),
        )

    ax.text(0.0, 0.95, "A", fontsize=18, fontweight="bold", color="#111111")
    ax.text(0.05, 0.95, "Frozen derivation-to-validation workflow", fontsize=15, fontweight="bold", color="#111111")


def draw_cohort_table(ax: plt.Axes, cohort_rows: list[dict[str, str]]) -> None:
    ax.axis("off")
    ax.text(0.0, 1.04, "B", fontsize=18, fontweight="bold", color="#111111", transform=ax.transAxes)
    ax.text(0.1, 1.04, "Cohort roles and endpoint context", fontsize=14, fontweight="bold", color="#111111", transform=ax.transAxes)

    headers = ["Cohort", "Role", "N"]
    role_map = {
        "Primary derivation dataset": "Derivation",
        "External prognosis replication": "Replication",
        "External prognosis replication or optional extension": "Replication",
        "Orthogonal validation and clinicopathologic context": "Orthogonal",
        "Treatment-aware external extension candidate": "Tx extension",
    }
    size_map = {
        "108 series samples; 93 KUGH tumor rows in recovered endpoint workbook": "93",
        "109 tumor samples with KUCM patient-level endpoint linkage": "109",
        "483": "433 used",
        "200 series samples; 192 gastric adenocarcinoma in final publication": "192",
        "443 clinical cases": "443",
        "443 clinical cases; smaller overlap for RNA-seq analyses": "443",
    }
    table_rows = [
        [
            row["dataset_id"],
            role_map.get(row["planned_role"], row["planned_role"]),
            size_map.get(row["sample_size"], row["sample_size"]),
        ]
        for row in cohort_rows
    ]
    table = ax.table(cellText=table_rows, colLabels=headers, loc="center", cellLoc="left", colLoc="left")
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1.0, 1.6)
    table.auto_set_column_width(col=list(range(len(headers))))
    for (row_index, col_index), cell in table.get_celld().items():
        cell.set_edgecolor("#fbfaf7")
        cell.set_linewidth(0.0)
        if row_index == 0:
            cell.set_facecolor("#ece7dd")
            cell.set_text_props(weight="bold", color="#222222")
        else:
            cell.set_facecolor("#f4efe5" if row_index % 2 == 1 else "#fbfaf7")
            if col_index == 1:
                cell.set_text_props(color="#4a4a4a")
            if col_index == 1:
                cell.get_text().set_fontsize(9.1)



def draw_projection_space(ax: plt.Axes, projection_space_rows: list[dict[str, str]]) -> None:
    ax.set_facecolor("#fbfaf7")
    ax.axhline(0.0, color="#d0c7bb", linewidth=1.0, linestyle="--", zorder=0)
    ax.axvline(0.0, color="#d0c7bb", linewidth=1.0, linestyle="--", zorder=0)
    for subtype in SUBTYPE_ORDER:
        subtype_rows = [row for row in projection_space_rows if row["predicted_subtype"] == subtype]
        if not subtype_rows:
            continue
        ax.scatter(
            [float(row["x_score"]) for row in subtype_rows],
            [float(row["y_score"]) for row in subtype_rows],
            s=[28 + 110 * float(row["assignment_margin"]) for row in subtype_rows],
            color=SUBTYPE_COLORS[subtype],
            alpha=[CONFIDENCE_ALPHA.get(row["confidence_tier"], 0.45) for row in subtype_rows],
            edgecolors="white",
            linewidths=0.5,
            label=subtype,
        )

    ax.text(0.02, 1.04, "C", fontsize=18, fontweight="bold", color="#111111", transform=ax.transAxes)
    ax.text(0.12, 1.04, "Centroid space", fontsize=13, fontweight="bold", color="#111111", transform=ax.transAxes)
    ax.set_xlabel("EMT score minus MSI score")
    ax.set_ylabel("MSS/TP53- score minus MSS/TP53+ score")
    ax.spines[["top", "right"]].set_visible(False)
    ax.text(0.05, 0.08, "MSI-favored", transform=ax.transAxes, fontsize=9.5, color="#174c7e")
    ax.text(0.72, 0.08, "EMT-favored", transform=ax.transAxes, fontsize=9.5, color="#c06c2b")
    ax.text(0.05, 0.92, "MSS/TP53-", transform=ax.transAxes, fontsize=9.5, color="#1f8a70")
    ax.text(0.7, 0.92, "MSS/TP53+", transform=ax.transAxes, fontsize=9.5, color="#7b5ea7")
    legend_handles = [
        Line2D([0], [0], marker="o", color="w", label=subtype, markerfacecolor=SUBTYPE_COLORS[subtype], markersize=8)
        for subtype in SUBTYPE_ORDER
    ]
    ax.legend(handles=legend_handles, loc="upper center", bbox_to_anchor=(0.5, -0.18), frameon=False, ncol=2, fontsize=8.8)


def draw_projection_qc(ax: plt.Axes, projection_rows: list[dict[str, str]]) -> None:
    subtype_rows = [row for row in projection_rows if row["predicted_subtype"] in SUBTYPE_ORDER]
    subtype_rows.sort(key=lambda row: SUBTYPE_ORDER.index(row["predicted_subtype"]))

    subtypes = [row["predicted_subtype"] for row in subtype_rows]
    margins = [float(row["median_assignment_margin"]) for row in subtype_rows]
    entropies = [float(row["median_assignment_entropy"]) for row in subtype_rows]

    x_positions = list(range(len(subtypes)))
    margin_bars = ax.bar(x_positions, margins, width=0.58, color="#c06c2b", alpha=0.9, label="Median margin")
    ax.set_ylabel("Median assignment margin", color="#7a3f12", fontsize=11)
    ax.set_xticks(x_positions, subtypes, rotation=20)
    ax.tick_params(axis="y", colors="#7a3f12")
    ax.spines[["top", "right"]].set_visible(False)

    second_axis = ax.twinx()
    second_axis.plot(x_positions, entropies, color="#174c7e", marker="o", linewidth=2.2, label="Median entropy")
    second_axis.set_ylabel("Median normalized entropy", color="#174c7e", fontsize=11)
    second_axis.tick_params(axis="y", colors="#174c7e")
    second_axis.spines["top"].set_visible(False)

    shared_gene_count = subtype_rows[0]["shared_gene_count"] if subtype_rows else "NA"
    low_margin_total = sum(int(row["low_margin_count"]) for row in subtype_rows)
    ax.text(0.02, 1.04, "D", fontsize=18, fontweight="bold", color="#111111", transform=ax.transAxes)
    ax.text(0.12, 1.04, "Projection QC", fontsize=13, fontweight="bold", color="#111111", transform=ax.transAxes)
    ax.text(
        0.02,
        0.98,
        f"{shared_gene_count} shared genes",
        fontsize=9.0,
        color="#4a4a4a",
        transform=ax.transAxes,
    )
    ax.text(
        0.02,
        0.92,
        f"Low-margin = {low_margin_total}",
        fontsize=9.0,
        color="#4a4a4a",
        transform=ax.transAxes,
    )

    lines = [margin_bars, second_axis.lines[0]]
    labels = ["Median margin", "Median entropy"]
    ax.legend(lines, labels, loc="upper left", bbox_to_anchor=(0.0, 0.78), frameon=False, fontsize=8.6)


def render_figure1(
    flow_rows: list[dict[str, str]],
    cohort_rows: list[dict[str, str]],
    projection_rows: list[dict[str, str]],
    projection_space_rows: list[dict[str, str]],
    output_pdf: Path,
    output_png: Path,
) -> None:
    output_pdf.parent.mkdir(parents=True, exist_ok=True)
    output_png.parent.mkdir(parents=True, exist_ok=True)

    figure = plt.figure(figsize=(15.8, 8.6))
    figure.patch.set_facecolor("#fbfaf7")
    grid = GridSpec(2, 3, figure=figure, height_ratios=[0.86, 1.1], width_ratios=[1.2, 0.92, 0.88], hspace=0.42, wspace=0.24)

    workflow_axis = figure.add_subplot(grid[0, :])
    cohort_axis = figure.add_subplot(grid[1, 0])
    projection_space_axis = figure.add_subplot(grid[1, 1])
    qc_axis = figure.add_subplot(grid[1, 2])

    draw_workflow(workflow_axis, flow_rows)
    draw_cohort_table(cohort_axis, cohort_rows)
    draw_projection_space(projection_space_axis, projection_space_rows)
    draw_projection_qc(qc_axis, projection_rows)

    figure.suptitle(
        "Study design, cohort qualification, and technical consistency",
        x=0.5,
        y=0.98,
        ha="center",
        fontsize=17,
        fontweight="bold",
        color="#111111",
    )
    figure.text(
        0.5,
        0.935,
        "Frozen centroids, explicit cohort roles, and projection-confidence summaries define the validation engine.",
        ha="center",
        fontsize=11,
        color="#4a4a4a",
    )
    figure.savefig(output_pdf, bbox_inches="tight")
    figure.savefig(output_png, dpi=300, bbox_inches="tight")
    plt.close(figure)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Render submission Figure 1 workflow and consistency panel")
    parser.add_argument("--flow-input", required=True)
    parser.add_argument("--cohort-map-input", required=True)
    parser.add_argument("--projection-input", required=True)
    parser.add_argument("--projection-space-input", required=True)
    parser.add_argument("--pdf-output", required=True)
    parser.add_argument("--png-output", required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    render_figure1(
        read_tsv(Path(args.flow_input)),
        read_tsv(Path(args.cohort_map_input)),
        read_tsv(Path(args.projection_input)),
        read_tsv(Path(args.projection_space_input)),
        Path(args.pdf_output),
        Path(args.png_output),
    )


if __name__ == "__main__":
    main()
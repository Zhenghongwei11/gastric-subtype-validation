from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from lifelines import KaplanMeierFitter


PALETTE = {
    "GSE62254": "#174c7e",
    "GSE15459": "#1f8a70",
    "GSE84437": "#c06c2b",
    "GSE26899": "#9a3b3b",
    "GSE26901": "#7b5ea7",
    "summary_diamond": "#3a3a3a",
    "NatCommFamily": "#3a3a3a",
}
SUBTYPE_COLORS = {
    "EMT": "#c06c2b",
    "MSI": "#174c7e",
}
TREATMENT_COLORS = {
    "1": "#1f8a70",
    "0": "#7d7d7d",
}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def format_ci(row: dict[str, str]) -> str:
    return f"{row['effect']} ({row['ci_lower']}-{row['ci_upper']})"


def build_figure3_display_rows(
    anchor_rows: list[dict[str, str]],
    summary_rows: list[dict[str, str]],
) -> list[dict[str, str]]:
    rows = list(anchor_rows)
    for row in summary_rows:
        rows.append(
            {
                **row,
                "dataset_id": "summary_diamond",
                "covariates": row.get("covariates", "fixed-effect summary"),
                "display_label": f"Summary (k={row.get('source_dataset_count', '0')})",
                "row_style": "summary",
            }
        )
    return rows


def compact_covariates(value: str) -> str:
    mapping = {
        "pt_stage,pn_stage": "pT + pN",
        "subgroup,adjuvant,subgroup_x_adjuvant": "subgroup + adj + int",
        "fixed-effect summary": "Fixed effect",
    }
    return mapping.get(value, value)


def format_family_summary_lines(family_rows: list[dict[str, str]]) -> list[str]:
    if not family_rows:
        return []

    outcome_labels = {
        "overall_survival": "OS",
        "recurrence_free_survival": "RFS",
    }
    lines = []
    for row in family_rows:
        outcome = outcome_labels.get(row["outcome"], row["outcome"])
        lines.append(
            (
                f"{outcome}: pooled HR {row['pooled_effect']} ({row['pooled_ci_lower']}-{row['pooled_ci_upper']}), "
                f"P={row['pooled_pvalue']}; meta HR {row['meta_pooled_hr']} ({row['meta_ci_lower']}-{row['meta_ci_upper']}), "
                f"P={row['meta_pvalue']}, I2={row['i_squared']}%"
            )
        )
    return lines


def draw_forest_block(
    figure: plt.Figure,
    outer_spec,
    rows: list[dict[str, str]],
) -> None:
    grid = outer_spec.subgridspec(1, 5, width_ratios=[1.7, 2.4, 2.2, 2.2, 1.6], wspace=0.03)
    cohort_axis = figure.add_subplot(grid[0, 0])
    hr_axis = figure.add_subplot(grid[0, 1], sharey=cohort_axis)
    forest_axis = figure.add_subplot(grid[0, 2], sharey=cohort_axis)
    size_axis = figure.add_subplot(grid[0, 3], sharey=cohort_axis)
    covariate_axis = figure.add_subplot(grid[0, 4], sharey=cohort_axis)

    for axis in [cohort_axis, hr_axis, forest_axis, size_axis, covariate_axis]:
        axis.set_facecolor("#fbfaf7")

    y_positions = list(range(len(rows), 0, -1))
    for y_position, row in zip(y_positions, rows, strict=True):
        color = PALETTE.get(row["dataset_id"], "#3a3a3a")
        effect = float(row["effect"])
        ci_lower = float(row["ci_lower"])
        ci_upper = float(row["ci_upper"])
        marker = "D" if row.get("row_style") == "summary" else "s"
        marker_size = 10.5 if row.get("row_style") == "summary" else 9
        forest_axis.errorbar(
            effect,
            y_position,
            xerr=[[effect - ci_lower], [ci_upper - effect]],
            fmt=marker,
            color=color,
            markerfacecolor=color,
            markeredgecolor="white",
            markersize=marker_size,
            elinewidth=2.2,
            capsize=0,
            zorder=3,
        )
        label = row.get("display_label", row["dataset_id"])
        cohort_axis.text(
            0.02,
            y_position,
            label,
            va="center",
            ha="left",
            fontsize=11.3 if row.get("row_style") == "summary" else 12.2,
            fontweight="bold",
            color="#1f1f1f",
        )
        hr_axis.text(0.02, y_position, format_ci(row), va="center", ha="left", fontsize=10.6, color="#2e2e2e")
        size_axis.text(0.02, y_position, f"n={row['n']}, events={row['events']}", va="center", ha="left", fontsize=10.6, color="#4a4a4a")
        covariate_axis.text(0.02, y_position, compact_covariates(row["covariates"]), va="center", ha="left", fontsize=10.2, color="#4a4a4a")

    forest_axis.axvline(1.0, color="#767676", linestyle="--", linewidth=1.1, zorder=1)
    forest_axis.set_xscale("log")
    forest_axis.set_xlim(0.75, max(float(row["ci_upper"]) for row in rows) * 1.18)
    forest_axis.set_ylim(0.4, len(rows) + 1.1)
    forest_axis.set_yticks([])
    forest_axis.set_xlabel("Hazard ratio", fontsize=12, color="#1f1f1f")

    cohort_axis.text(0.02, len(rows) + 0.48, "Cohort", fontsize=10.3, color="#585858", fontweight="bold")
    hr_axis.text(0.02, len(rows) + 0.48, "HR (95% CI)", fontsize=10.3, color="#585858", fontweight="bold")
    size_axis.text(0.02, len(rows) + 0.48, "Cohort size", fontsize=10.3, color="#585858", fontweight="bold")
    covariate_axis.text(0.02, len(rows) + 0.48, "Covariates", fontsize=10.3, color="#585858", fontweight="bold")

    for axis in [cohort_axis, hr_axis, size_axis, covariate_axis]:
        axis.set_xlim(0, 1)
        axis.set_ylim(0.4, len(rows) + 1.1)
        axis.axis("off")

    for spine_name in ["top", "right", "left"]:
        forest_axis.spines[spine_name].set_visible(False)
    forest_axis.spines["bottom"].set_color("#7a7a7a")
    forest_axis.tick_params(axis="x", labelsize=10.5, colors="#333333")


def draw_subtype_km_panel(ax: plt.Axes, rows: list[dict[str, str]], title: str) -> None:
    frame = pd.DataFrame(rows)
    frame["overall_survival_months"] = pd.to_numeric(frame["overall_survival_months"])
    frame["overall_survival_event"] = pd.to_numeric(frame["overall_survival_event"])
    km_fitter = KaplanMeierFitter()
    for subtype in ["EMT", "MSI"]:
        subframe = frame[frame["subtype"] == subtype]
        if subframe.empty:
            continue
        km_fitter.fit(
            durations=subframe["overall_survival_months"],
            event_observed=subframe["overall_survival_event"],
            label=f"{subtype} (n={len(subframe)})",
        )
        km_fitter.plot(ax=ax, ci_show=False, color=SUBTYPE_COLORS[subtype], linewidth=2.2)

    ax.set_title(title, loc="left", fontsize=13, fontweight="bold")
    ax.set_xlabel("Months")
    ax.set_ylabel("Survival probability")
    ax.set_facecolor("#fbfaf7")
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False, loc="lower left", fontsize=9)


def draw_treatment_km_panel(ax: plt.Axes, rows: list[dict[str, str]], subgroup: str, title: str) -> None:
    frame = pd.DataFrame([row for row in rows if row["subgroup"] == subgroup])
    frame["overall_survival_months"] = pd.to_numeric(frame["overall_survival_months"])
    frame["overall_survival_event"] = pd.to_numeric(frame["overall_survival_event"])
    km_fitter = KaplanMeierFitter()
    labels = {"1": "Adjuvant", "0": "No adjuvant"}
    for treatment_flag in ["1", "0"]:
        subframe = frame[frame["adjuvant_chemotherapy_binary"] == treatment_flag]
        if subframe.empty:
            continue
        km_fitter.fit(
            durations=subframe["overall_survival_months"],
            event_observed=subframe["overall_survival_event"],
            label=f"{labels[treatment_flag]} (n={len(subframe)})",
        )
        km_fitter.plot(ax=ax, ci_show=False, color=TREATMENT_COLORS[treatment_flag], linewidth=2.2)

    ax.set_title(title, loc="left", fontsize=13, fontweight="bold")
    ax.set_xlabel("Months")
    ax.set_ylabel("Survival probability")
    ax.set_facecolor("#fbfaf7")
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False, loc="lower left", fontsize=9)


def draw_text_block(ax: plt.Axes, header: str, lines: list[str]) -> None:
    ax.axis("off")
    ax.text(0.5, 0.95, header, ha="center", fontsize=10.6, fontweight="bold", color="#333333", transform=ax.transAxes)
    for index, line in enumerate(lines, start=1):
        ax.text(0.5, 0.95 - 0.28 * index, line, ha="center", fontsize=9.9, color="#4a4a4a", transform=ax.transAxes)


def render_figure3(
    figure3_rows: list[dict[str, str]],
    figure3_summary_rows: list[dict[str, str]],
    figure3_km_rows: list[dict[str, str]],
    output_pdf: Path,
    output_png: Path,
) -> None:
    figure = plt.figure(figsize=(13.8, 10.4))
    figure.patch.set_facecolor("#fbfaf7")
    grid = figure.add_gridspec(2, 2, height_ratios=[1.32, 1.0], hspace=0.42, wspace=0.28)

    draw_forest_block(figure, grid[0, :], build_figure3_display_rows(figure3_rows, figure3_summary_rows))
    draw_subtype_km_panel(
        figure.add_subplot(grid[1, 0]),
        [row for row in figure3_km_rows if row["dataset_id"] == "GSE15459"],
        "B  GSE15459 external validation KM",
    )
    draw_subtype_km_panel(
        figure.add_subplot(grid[1, 1]),
        [row for row in figure3_km_rows if row["dataset_id"] == "GSE84437"],
        "C  GSE84437 projected-state KM",
    )

    figure.text(0.5, 0.98, "External prognosis replication across independent cohorts", ha="center", va="top", fontsize=17, fontweight="bold", color="#111111")
    figure.text(0.5, 0.94, "A  Forest and fixed-effect summary diamond", ha="center", fontsize=14, fontweight="bold", color="#111111")
    figure.text(
        0.5,
        0.915,
        "Frozen shared ACRG state. Cohort-wise HRs remain primary; the fixed-effect summary diamond is subordinate context.",
        ha="center",
        fontsize=10.8,
        color="#4a4a4a",
    )
    figure.savefig(output_pdf, bbox_inches="tight")
    figure.savefig(output_png, dpi=300, bbox_inches="tight")
    plt.close(figure)


def render_figure4(
    figure4_rows: list[dict[str, str]],
    figure4_stage_km_rows: list[dict[str, str]],
    output_pdf: Path,
    output_png: Path,
) -> None:
    figure = plt.figure(figsize=(13.2, 9.8))
    figure.patch.set_facecolor("#fbfaf7")
    grid = figure.add_gridspec(2, 1, height_ratios=[1.3, 0.98], hspace=0.4)

    draw_forest_block(figure, grid[0, 0], figure4_rows)
    draw_subtype_km_panel(
        figure.add_subplot(grid[1, 0]),
        figure4_stage_km_rows,
        "B  Stage III derivation-cohort stratified KM",
    )

    figure.text(0.5, 0.98, "Clinicopathologic sensitivity analyses for the primary EMT-versus-MSI contrast", ha="center", va="top", fontsize=17, fontweight="bold", color="#111111")
    figure.text(0.5, 0.94, "A  Cohort-appropriate multivariable forest", ha="center", fontsize=14, fontweight="bold", color="#111111")
    figure.text(
        0.5,
        0.915,
        "Adjusted models retain the same adverse EMT direction under available stage covariates, while Stage III derivation KM shows within-stage separation.",
        ha="center",
        fontsize=10.8,
        color="#4a4a4a",
    )
    figure.savefig(output_pdf, bbox_inches="tight")
    figure.savefig(output_png, dpi=300, bbox_inches="tight")
    plt.close(figure)


def render_figure5(
    figure5_rows: list[dict[str, str]],
    figure5_family_rows: list[dict[str, str]],
    figure5_cross_km_rows: list[dict[str, str]],
    output_pdf: Path,
    output_png: Path,
) -> None:
    family_summary_lines = format_family_summary_lines(figure5_family_rows)
    figure = plt.figure(figsize=(13.6, 11.0))
    figure.patch.set_facecolor("#fbfaf7")
    grid = figure.add_gridspec(3, 2, height_ratios=[1.26, 0.28, 1.0], hspace=0.34, wspace=0.3)

    draw_forest_block(figure, grid[0, :], figure5_rows)
    draw_text_block(
        figure.add_subplot(grid[1, :]),
        "Family-level pooled/meta summary",
        family_summary_lines,
    )
    draw_treatment_km_panel(
        figure.add_subplot(grid[2, 0]),
        figure5_cross_km_rows,
        subgroup="EP",
        title="B  GSE26899 EP: adjuvant vs no adjuvant",
    )
    draw_treatment_km_panel(
        figure.add_subplot(grid[2, 1]),
        figure5_cross_km_rows,
        subgroup="MP",
        title="C  GSE26899 MP: adjuvant vs no adjuvant",
    )

    figure.text(0.5, 0.98, "Nat Commun cohort-family treatment-extension interaction effects", ha="center", va="top", fontsize=17, fontweight="bold", color="#111111")
    figure.text(0.5, 0.94, "A  Cohort-wise interaction forest", ha="center", fontsize=14, fontweight="bold", color="#111111")
    figure.text(
        0.5,
        0.915,
        "Interaction estimates remain primary. GSE26899 cross-KM panels are exploratory support for differential adjuvant benefit by subgroup.",
        ha="center",
        fontsize=10.8,
        color="#4a4a4a",
    )
    figure.savefig(output_pdf, bbox_inches="tight")
    figure.savefig(output_png, dpi=300, bbox_inches="tight")
    plt.close(figure)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Render balanced manuscript-grade prognosis figures")
    parser.add_argument("--figure3-anchor-input", required=True)
    parser.add_argument("--figure3-summary-input", required=True)
    parser.add_argument("--figure3-km-input", required=True)
    parser.add_argument("--figure4-anchor-input", required=True)
    parser.add_argument("--figure4-stage-km-input", required=True)
    parser.add_argument("--figure3-pdf-output", required=True)
    parser.add_argument("--figure3-png-output", required=True)
    parser.add_argument("--figure4-pdf-output", required=True)
    parser.add_argument("--figure4-png-output", required=True)
    parser.add_argument("--figure5-anchor-input", required=True)
    parser.add_argument("--figure5-family-input", required=True)
    parser.add_argument("--figure5-cross-km-input", required=True)
    parser.add_argument("--figure5-pdf-output", required=True)
    parser.add_argument("--figure5-png-output", required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    render_figure3(
        read_tsv(Path(args.figure3_anchor_input)),
        read_tsv(Path(args.figure3_summary_input)),
        read_tsv(Path(args.figure3_km_input)),
        Path(args.figure3_pdf_output),
        Path(args.figure3_png_output),
    )
    render_figure4(
        read_tsv(Path(args.figure4_anchor_input)),
        read_tsv(Path(args.figure4_stage_km_input)),
        Path(args.figure4_pdf_output),
        Path(args.figure4_png_output),
    )
    render_figure5(
        read_tsv(Path(args.figure5_anchor_input)),
        read_tsv(Path(args.figure5_family_input)),
        read_tsv(Path(args.figure5_cross_km_input)),
        Path(args.figure5_pdf_output),
        Path(args.figure5_png_output),
    )


if __name__ == "__main__":
    main()
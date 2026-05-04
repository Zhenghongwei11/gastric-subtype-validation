from __future__ import annotations

import argparse
import csv
from pathlib import Path

import pandas as pd
from lifelines import CoxPHFitter


DATASET_ORDER = ["GSE26899", "GSE26901"]
SUBGROUP_REFERENCE = "EP"
SUBGROUP_COMPARISON = "MP"
OUTCOME_SPECS = [
    ("overall_survival", "overall_survival_months", "overall_survival_event"),
    ("recurrence_free_survival", "recurrence_free_survival_months", "recurrence_event"),
]


def format_number(value: float) -> str:
    return format(float(value), ".6g")


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def load_frame(path: Path, dataset_id: str) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t", dtype=str)
    frame = frame[frame["dataset_id"] == dataset_id].copy()
    frame = frame[frame["subgroup"].isin([SUBGROUP_REFERENCE, SUBGROUP_COMPARISON])].copy()
    frame = frame[frame["adjuvant_chemotherapy_binary"].isin(["0", "1"])].copy()
    frame["subgroup_binary"] = (frame["subgroup"] == SUBGROUP_COMPARISON).astype(int)
    frame["adjuvant_binary"] = pd.to_numeric(frame["adjuvant_chemotherapy_binary"])
    for _, duration_col, event_col in OUTCOME_SPECS:
        frame[duration_col] = pd.to_numeric(frame[duration_col], errors="coerce")
        frame[event_col] = pd.to_numeric(frame[event_col], errors="coerce")
    return frame


def fit_single_covariate_model(
    frame: pd.DataFrame,
    dataset_id: str,
    outcome: str,
    duration_col: str,
    event_col: str,
    adjuvant_value: int,
) -> dict[str, str]:
    model_name = "cox_subgroup_in_adjuvant_treated" if adjuvant_value == 1 else "cox_subgroup_in_adjuvant_untreated"
    model_frame = frame[frame["adjuvant_binary"] == adjuvant_value].dropna(subset=[duration_col, event_col]).copy()

    fitter = CoxPHFitter(penalizer=0.05)
    fitter.fit(
        model_frame[[duration_col, event_col, "subgroup_binary"]],
        duration_col=duration_col,
        event_col=event_col,
    )
    summary = fitter.summary.loc["subgroup_binary"]
    return {
        "claim_id": "C2",
        "dataset_id": dataset_id,
        "outcome": outcome,
        "model": model_name,
        "reference_subtype": SUBGROUP_REFERENCE,
        "comparison_subtype": SUBGROUP_COMPARISON,
        "effect_type": "hazard_ratio",
        "effect": format_number(summary["exp(coef)"]),
        "ci_lower": format_number(summary["exp(coef) lower 95%"]),
        "ci_upper": format_number(summary["exp(coef) upper 95%"]),
        "pvalue": format_number(summary["p"]),
        "fdr": "",
        "n": str(len(model_frame)),
        "events": str(int(model_frame[event_col].sum())),
        "covariates": "subgroup_only",
        "adjuvant_stratum": str(adjuvant_value),
    }


def fit_interaction_model(
    frame: pd.DataFrame,
    dataset_id: str,
    outcome: str,
    duration_col: str,
    event_col: str,
) -> dict[str, str]:
    model_frame = frame.dropna(subset=[duration_col, event_col]).copy()
    model_frame["subgroup_adjuvant_interaction"] = model_frame["subgroup_binary"] * model_frame["adjuvant_binary"]

    fitter = CoxPHFitter(penalizer=0.05)
    fitter.fit(
        model_frame[
            [duration_col, event_col, "subgroup_binary", "adjuvant_binary", "subgroup_adjuvant_interaction"]
        ],
        duration_col=duration_col,
        event_col=event_col,
    )
    summary = fitter.summary.loc["subgroup_adjuvant_interaction"]
    return {
        "claim_id": "C2",
        "dataset_id": dataset_id,
        "outcome": outcome,
        "model": "cox_subgroup_adjuvant_interaction",
        "reference_subtype": SUBGROUP_REFERENCE,
        "comparison_subtype": SUBGROUP_COMPARISON,
        "effect_type": "interaction_hazard_ratio",
        "effect": format_number(summary["exp(coef)"]),
        "ci_lower": format_number(summary["exp(coef) lower 95%"]),
        "ci_upper": format_number(summary["exp(coef) upper 95%"]),
        "pvalue": format_number(summary["p"]),
        "fdr": "",
        "n": str(len(model_frame)),
        "events": str(int(model_frame[event_col].sum())),
        "covariates": "subgroup,adjuvant,subgroup_x_adjuvant",
        "adjuvant_stratum": "interaction",
    }


def build_effect_rows(gse26899_input: Path, gse26901_input: Path) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    frames = {
        "GSE26899": load_frame(gse26899_input, "GSE26899"),
        "GSE26901": load_frame(gse26901_input, "GSE26901"),
    }
    for dataset_id in DATASET_ORDER:
        frame = frames[dataset_id]
        for outcome, duration_col, event_col in OUTCOME_SPECS:
            rows.append(fit_single_covariate_model(frame, dataset_id, outcome, duration_col, event_col, 1))
            rows.append(fit_single_covariate_model(frame, dataset_id, outcome, duration_col, event_col, 0))
            rows.append(fit_interaction_model(frame, dataset_id, outcome, duration_col, event_col))
    return rows


def build_figure_rows(effect_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    return [row for row in effect_rows if row["model"] == "cox_subgroup_adjuvant_interaction"]


def build_summary_rows(gse26899_input: Path, gse26901_input: Path) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for dataset_id, input_path in [("GSE26899", gse26899_input), ("GSE26901", gse26901_input)]:
        frame = load_frame(input_path, dataset_id)
        grouped = frame.groupby(["subgroup", "adjuvant_chemotherapy_binary"], sort=True)
        for (subgroup, adjuvant_binary), subframe in grouped:
            rows.append(
                {
                    "dataset_id": dataset_id,
                    "subgroup": str(subgroup),
                    "adjuvant_chemotherapy_binary": str(adjuvant_binary),
                    "n": str(len(subframe)),
                    "overall_survival_events": str(int(subframe["overall_survival_event"].sum())),
                    "recurrence_free_survival_events": str(int(subframe["recurrence_event"].sum())),
                    "median_overall_survival_months": format_number(subframe["overall_survival_months"].median()),
                    "median_recurrence_free_survival_months": format_number(subframe["recurrence_free_survival_months"].median()),
                }
            )
    return rows


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run Nat Commun cohort-family C2 survival and treatment-heterogeneity analyses")
    parser.add_argument("--gse26899-input", required=True)
    parser.add_argument("--gse26901-input", required=True)
    parser.add_argument("--effect-output", required=True)
    parser.add_argument("--figure-output", required=True)
    parser.add_argument("--summary-output", required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    gse26899_input = Path(args.gse26899_input)
    gse26901_input = Path(args.gse26901_input)
    effect_rows = build_effect_rows(gse26899_input, gse26901_input)
    write_tsv(Path(args.effect_output), effect_rows)
    write_tsv(Path(args.figure_output), build_figure_rows(effect_rows))
    write_tsv(Path(args.summary_output), build_summary_rows(gse26899_input, gse26901_input))


if __name__ == "__main__":
    main()
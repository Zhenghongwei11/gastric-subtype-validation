from __future__ import annotations

import argparse
import csv
from pathlib import Path

import pandas as pd
from lifelines import CoxPHFitter


REFERENCE_SUBTYPE = "Metabolic"
COMPARISON_SUBTYPES = ["Invasive", "Proliferative", "Unstable"]


def format_number(value: float) -> str:
    return format(float(value), ".6g")


def load_outcomes(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t")
    frame["overall_survival_months"] = pd.to_numeric(frame["overall_survival_months"])
    frame["overall_survival_event"] = pd.to_numeric(frame["overall_survival_event"])
    frame["stage"] = pd.to_numeric(frame["stage"])
    frame = frame[frame["subtype"].isin([REFERENCE_SUBTYPE, *COMPARISON_SUBTYPES])].copy()
    frame["subtype"] = pd.Categorical(
        frame["subtype"], categories=[REFERENCE_SUBTYPE, *COMPARISON_SUBTYPES], ordered=False
    )
    return frame


def fit_model(frame: pd.DataFrame, model_name: str) -> list[dict[str, str]]:
    covariate_columns = []
    model_frame = frame[["overall_survival_months", "overall_survival_event", "subtype", "stage"]].copy()
    subtype_dummies = pd.get_dummies(model_frame["subtype"], prefix="subtype")
    reference_column = f"subtype_{REFERENCE_SUBTYPE}"
    if reference_column in subtype_dummies.columns:
        subtype_dummies = subtype_dummies.drop(columns=[reference_column])
    covariate_columns.extend(sorted(subtype_dummies.columns))
    model_frame = pd.concat([model_frame, subtype_dummies], axis=1)

    if model_name == "cox_stage_adjusted":
        covariate_columns.append("stage")

    cox_frame = model_frame[["overall_survival_months", "overall_survival_event", *covariate_columns]].copy()
    fitter = CoxPHFitter()
    fitter.fit(
        cox_frame,
        duration_col="overall_survival_months",
        event_col="overall_survival_event",
    )

    rows: list[dict[str, str]] = []
    for comparison in COMPARISON_SUBTYPES:
        coefficient_name = f"subtype_{comparison}"
        summary = fitter.summary.loc[coefficient_name]
        rows.append(
            {
                "dataset_id": "GSE15459",
                "outcome": "overall_survival",
                "model": model_name,
                "reference_subtype": REFERENCE_SUBTYPE,
                "comparison_subtype": comparison,
                "effect_type": "hazard_ratio",
                "effect": format_number(summary["exp(coef)"]),
                "ci_lower": format_number(summary["exp(coef) lower 95%"]),
                "ci_upper": format_number(summary["exp(coef) upper 95%"]),
                "pvalue": format_number(summary["p"]),
                "n": str(len(frame)),
                "events": str(int(frame["overall_survival_event"].sum())),
                "covariates": "stage" if model_name == "cox_stage_adjusted" else "none",
            }
        )
    return rows


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run GSE15459 subtype survival models")
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    frame = load_outcomes(Path(args.input))
    rows = [
        *fit_model(frame, "cox_unadjusted"),
        *fit_model(frame, "cox_stage_adjusted"),
    ]
    write_tsv(Path(args.output), rows)


if __name__ == "__main__":
    main()
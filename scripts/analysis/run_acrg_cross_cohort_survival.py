from __future__ import annotations

import argparse
import csv
from pathlib import Path

import pandas as pd
import xlrd
from lifelines import CoxPHFitter


REFERENCE_SUBTYPE = "MSI"
COMPARISON_SUBTYPES = ["EMT", "MSS/TP53-", "MSS/TP53+"]
SUBTYPE_LABELS = {
    "0": "MSS/TP53-",
    "1": "MSS/TP53+",
    "2": "MSI",
    "3": "EMT",
}


def format_number(value: float) -> str:
    return format(float(value), ".6g")


def normalize_stage_value(raw_value: str) -> int | None:
    value = str(raw_value).strip().upper()
    if not value:
        return None

    roman_map = {
        "I": 1,
        "II": 2,
        "III": 3,
        "IV": 4,
    }
    for token in ("IV", "III", "II", "I"):
        if token in value:
            return roman_map[token]

    for token in ("4", "3", "2", "1"):
        if token in value:
            return int(token)

    return None


def parse_source_sheet(path: Path, sheet_name: str) -> dict[str, dict[str, str]]:
    workbook = xlrd.open_workbook(str(path))
    sheet = workbook.sheet_by_name(sheet_name)
    headers = [str(value).strip() for value in sheet.row_values(0)]
    sample_idx = headers.index("Sample ID")
    os_idx = headers.index("OS (mos)")
    censor_idx = headers.index("OS censor")
    subtype_idx = headers.index("Mol subtype")

    rows: dict[str, dict[str, str]] = {}
    for row_index in range(1, sheet.nrows):
        sample_id = str(sheet.cell_value(row_index, sample_idx)).strip()
        if sample_id.endswith(".0"):
            sample_id = sample_id[:-2]
        subtype_code = str(int(sheet.cell_value(row_index, subtype_idx)))
        rows[sample_id] = {
            "overall_survival_months": format_number(sheet.cell_value(row_index, os_idx)),
            "overall_survival_event": str(int(1 - float(sheet.cell_value(row_index, censor_idx)))),
            "molecular_subtype_code": subtype_code,
            "molecular_subtype_label": SUBTYPE_LABELS[subtype_code],
        }
    return rows


def load_acrg_frame(acrg_input: Path, source_data_input: Path) -> pd.DataFrame:
    frame = pd.read_csv(acrg_input, sep="\t", dtype=str)
    source_rows = parse_source_sheet(source_data_input, "ACRG")

    frame["patient_id"] = frame["patient_id"].astype(str)
    frame = frame[frame["patient_id"].isin(source_rows)].copy()
    frame["overall_survival_months"] = frame["patient_id"].map(
        lambda patient_id: source_rows[patient_id]["overall_survival_months"]
    )
    frame["overall_survival_event"] = frame["patient_id"].map(
        lambda patient_id: source_rows[patient_id]["overall_survival_event"]
    )
    frame["subtype"] = frame["patient_id"].map(
        lambda patient_id: source_rows[patient_id]["molecular_subtype_label"]
    )
    frame["overall_survival_months"] = pd.to_numeric(frame["overall_survival_months"])
    frame["overall_survival_event"] = pd.to_numeric(frame["overall_survival_event"])
    frame["stage"] = frame["pathologic_stage"].map(normalize_stage_value)
    frame["subtype"] = pd.Categorical(
        frame["subtype"], categories=[REFERENCE_SUBTYPE, *COMPARISON_SUBTYPES], ordered=False
    )
    return frame


def load_gse15459_frame(gse15459_input: Path, source_data_input: Path) -> pd.DataFrame:
    frame = pd.read_csv(gse15459_input, sep="\t", dtype=str)
    source_rows = parse_source_sheet(source_data_input, "Singapore")

    frame["sample_id"] = frame["cel_file"].str.replace(".CEL", "", regex=False)
    frame = frame[frame["sample_id"].isin(source_rows)].copy()
    frame["subtype"] = frame["sample_id"].map(
        lambda sample_id: source_rows[sample_id]["molecular_subtype_label"]
    )
    frame["overall_survival_months"] = pd.to_numeric(frame["overall_survival_months"])
    frame["overall_survival_event"] = pd.to_numeric(frame["overall_survival_event"])
    frame["stage"] = frame["stage"].map(normalize_stage_value)
    frame["subtype"] = pd.Categorical(
        frame["subtype"], categories=[REFERENCE_SUBTYPE, *COMPARISON_SUBTYPES], ordered=False
    )
    return frame


def fit_model(frame: pd.DataFrame, dataset_id: str, model_name: str) -> list[dict[str, str]]:
    covariate_columns: list[str] = []
    model_frame = frame[["overall_survival_months", "overall_survival_event", "subtype", "stage"]].copy()
    subtype_dummies = pd.get_dummies(model_frame["subtype"], prefix="subtype")
    reference_column = f"subtype_{REFERENCE_SUBTYPE}"
    if reference_column in subtype_dummies.columns:
        subtype_dummies = subtype_dummies.drop(columns=[reference_column])
    covariate_columns.extend(sorted(subtype_dummies.columns))
    model_frame = pd.concat([model_frame, subtype_dummies], axis=1)

    if model_name == "cox_stage_adjusted":
        model_frame = model_frame.dropna(subset=["stage"]).copy()
        model_frame["stage"] = pd.to_numeric(model_frame["stage"])
        covariate_columns.append("stage")

    cox_frame = model_frame[["overall_survival_months", "overall_survival_event", *covariate_columns]].copy()
    fitter = CoxPHFitter()
    fitter.fit(cox_frame, duration_col="overall_survival_months", event_col="overall_survival_event")

    rows: list[dict[str, str]] = []
    for comparison in COMPARISON_SUBTYPES:
        coefficient_name = f"subtype_{comparison}"
        summary = fitter.summary.loc[coefficient_name]
        rows.append(
            {
                "dataset_id": dataset_id,
                "outcome": "overall_survival",
                "model": model_name,
                "reference_subtype": REFERENCE_SUBTYPE,
                "comparison_subtype": comparison,
                "effect_type": "hazard_ratio",
                "effect": format_number(summary["exp(coef)"]),
                "ci_lower": format_number(summary["exp(coef) lower 95%"]),
                "ci_upper": format_number(summary["exp(coef) upper 95%"]),
                "pvalue": format_number(summary["p"]),
                "n": str(len(model_frame)),
                "events": str(int(model_frame["overall_survival_event"].sum())),
                "covariates": "stage" if model_name == "cox_stage_adjusted" else "none",
            }
        )
    return rows


def build_model_frame(frame: pd.DataFrame, model_name: str) -> pd.DataFrame:
    model_frame = frame[["overall_survival_months", "overall_survival_event", "subtype", "stage"]].copy()
    if model_name == "cox_stage_adjusted":
        model_frame = model_frame.dropna(subset=["stage"]).copy()
        model_frame["stage"] = pd.to_numeric(model_frame["stage"])
    return model_frame


def subtype_counts(model_frame: pd.DataFrame, dataset_id: str, model_name: str) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for subtype, subframe in model_frame.dropna(subset=["subtype"]).groupby("subtype", observed=True):
        rows.append(
            {
                "dataset_id": dataset_id,
                "outcome": "overall_survival",
                "model": model_name,
                "subtype": str(subtype),
                "n": str(len(subframe)),
                "events": str(int(pd.to_numeric(subframe["overall_survival_event"]).sum())),
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
    parser = argparse.ArgumentParser(description="Run aligned ACRG subtype survival models across cohorts")
    parser.add_argument("--acrg-input", required=True)
    parser.add_argument("--gse15459-input", required=True)
    parser.add_argument("--source-data-input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--subtype-counts-output")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    acrg_frame = load_acrg_frame(Path(args.acrg_input), Path(args.source_data_input))
    singapore_frame = load_gse15459_frame(Path(args.gse15459_input), Path(args.source_data_input))
    rows = [
        *fit_model(acrg_frame, "GSE62254", "cox_unadjusted"),
        *fit_model(acrg_frame, "GSE62254", "cox_stage_adjusted"),
        *fit_model(singapore_frame, "GSE15459", "cox_unadjusted"),
        *fit_model(singapore_frame, "GSE15459", "cox_stage_adjusted"),
    ]
    write_tsv(Path(args.output), rows)
    if args.subtype_counts_output:
        counts_rows = [
            *subtype_counts(build_model_frame(acrg_frame, "cox_unadjusted"), "GSE62254", "cox_unadjusted"),
            *subtype_counts(
                build_model_frame(acrg_frame, "cox_stage_adjusted"),
                "GSE62254",
                "cox_stage_adjusted",
            ),
            *subtype_counts(
                build_model_frame(singapore_frame, "cox_unadjusted"),
                "GSE15459",
                "cox_unadjusted",
            ),
            *subtype_counts(
                build_model_frame(singapore_frame, "cox_stage_adjusted"),
                "GSE15459",
                "cox_stage_adjusted",
            ),
        ]
        if counts_rows:
            write_tsv(Path(args.subtype_counts_output), counts_rows)


if __name__ == "__main__":
    main()

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import xlrd


DATASET_ID = "GSE15459"
SHEET_NAME = "GSE15460 Outcome Data"


def normalize_value(value: object) -> str:
    if value in (None, ""):
        return ""
    if isinstance(value, float):
        if value.is_integer():
            return str(int(value))
        return format(value, ".10g")
    return str(value).strip()


def is_sample_row(gsm_id: str) -> bool:
    return gsm_id.startswith("GSM")


def parse_rows(input_path: Path) -> tuple[list[dict[str, str]], dict[str, str]]:
    workbook = xlrd.open_workbook(str(input_path))
    sheet = workbook.sheet_by_name(SHEET_NAME)
    headers = [normalize_value(value) for value in sheet.row_values(0)]
    column_index = {name: index for index, name in enumerate(headers)}

    sample_rows: list[dict[str, str]] = []
    raw_sheet_rows = sheet.nrows - 1

    for row_index in range(1, sheet.nrows):
        values = [normalize_value(value) for value in sheet.row_values(row_index)]
        gsm_id = values[column_index["GSM ID"]]
        if not is_sample_row(gsm_id):
            continue

        sample_rows.append(
            {
                "dataset_id": DATASET_ID,
                "gsm_id": gsm_id,
                "patient_id": values[column_index["ID"]],
                "cel_file": values[column_index["Expression CEL file"]],
                "subtype": values[column_index["Subtype"]],
                "age_at_surgery": values[column_index["Age_at_surgery"]],
                "gender": values[column_index["Gender"]],
                "lauren_classification": values[column_index["Laurenclassification"]],
                "stage": values[column_index["Stage"]],
                "overall_survival_months": values[column_index["Overall.Survival (Months)**"]],
                "overall_survival_event": values[column_index["Outcome (1=dead)"]],
            }
        )

    summary_row = {
        "dataset_id": DATASET_ID,
        "endpoint_name": "overall_survival",
        "raw_sheet_rows": str(raw_sheet_rows),
        "valid_sample_rows": str(len(sample_rows)),
        "time_nonmissing_rows": str(sum(1 for row in sample_rows if row["overall_survival_months"])),
        "event_nonmissing_rows": str(sum(1 for row in sample_rows if row["overall_survival_event"])),
        "complete_case_rows": str(
            sum(
                1
                for row in sample_rows
                if row["overall_survival_months"] and row["overall_survival_event"]
            )
        ),
        "stage_nonmissing_rows": str(sum(1 for row in sample_rows if row["stage"])),
        "source_file": input_path.name,
    }
    return sample_rows, summary_row


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Parse GSE15459 clinical outcome workbook")
    parser.add_argument("--input", required=True)
    parser.add_argument("--outcomes-output", required=True)
    parser.add_argument("--summary-output", required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    outcome_rows, summary_row = parse_rows(Path(args.input))
    if not outcome_rows:
        raise SystemExit("No valid GSE15459 outcome rows were parsed")
    write_tsv(Path(args.outcomes_output), outcome_rows)
    write_tsv(Path(args.summary_output), [summary_row])


if __name__ == "__main__":
    main()
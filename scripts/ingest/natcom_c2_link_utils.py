from __future__ import annotations

import csv
import gzip
from pathlib import Path
from typing import Callable


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return [{key: (value or "") for key, value in row.items()} for row in reader]


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def parse_matrix_fields(line: str) -> list[str]:
    return next(csv.reader([line], delimiter="\t", quotechar='"'))


def load_matrix_metadata(matrix_input: Path, dataset_id: str) -> tuple[list[str], list[str], list[str]]:
    sample_titles: list[str] = []
    sample_accessions: list[str] = []
    patient_ids: list[str] = []

    with gzip.open(matrix_input, "rt", encoding="utf-8", errors="replace", newline="") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if line.startswith("!Sample_title"):
                sample_titles = parse_matrix_fields(line)[1:]
            elif line.startswith("!Sample_geo_accession"):
                sample_accessions = parse_matrix_fields(line)[1:]
            elif line.startswith("!Sample_characteristics_ch"):
                fields = parse_matrix_fields(line)[1:]
                if fields and fields[0].startswith("patient: "):
                    patient_ids = [value.replace("patient:", "", 1).strip() for value in fields]
            elif line == "!series_matrix_table_begin":
                break

    if not sample_titles or not sample_accessions or not patient_ids:
        raise SystemExit(f"Failed to parse required sample metadata from {dataset_id} matrix")
    if not (len(sample_titles) == len(sample_accessions) == len(patient_ids)):
        raise SystemExit(f"{dataset_id} matrix metadata lengths do not match")
    return sample_titles, sample_accessions, patient_ids


def load_clinical_rows(clinical_input: Path, dataset_id: str) -> dict[str, dict[str, str]]:
    rows = [row for row in read_tsv(clinical_input) if row.get("dataset_id") == dataset_id]
    if not rows:
        raise SystemExit(f"No standardized clinical rows available for {dataset_id}")

    rows_by_patient: dict[str, dict[str, str]] = {}
    for row in rows:
        patient_id = row["patient_id"]
        if patient_id in rows_by_patient:
            raise SystemExit(f"Duplicate patient_id {patient_id} in standardized clinical rows for {dataset_id}")
        rows_by_patient[patient_id] = row
    return rows_by_patient


def build_link_rows(
    dataset_id: str,
    sample_titles: list[str],
    sample_accessions: list[str],
    matrix_patient_ids: list[str],
    clinical_rows: dict[str, dict[str, str]],
    patient_normalizer: Callable[[str], str],
    link_method: str,
) -> tuple[list[dict[str, str]], dict[str, str], list[int]]:
    link_rows: list[dict[str, str]] = []
    selected_positions: list[int] = []

    for position, (gsm_id, sample_title, matrix_patient_id) in enumerate(
        zip(sample_accessions, sample_titles, matrix_patient_ids),
        start=1,
    ):
        patient_id = patient_normalizer(matrix_patient_id)
        clinical_row = clinical_rows.get(patient_id)
        if clinical_row is None:
            continue

        selected_positions.append(position)
        link_rows.append(
            {
                "dataset_id": dataset_id,
                "gsm_id": gsm_id,
                "patient_id": clinical_row["patient_id"],
                "matrix_patient_id": matrix_patient_id,
                "matrix_sample_title": sample_title,
                "source_sheet": clinical_row["source_sheet"],
                "subgroup": clinical_row["subgroup"],
                "sex": clinical_row["sex"],
                "age": clinical_row["age"],
                "location": clinical_row["location"],
                "lauren": clinical_row["lauren"],
                "ajcc_stage": clinical_row["ajcc_stage"],
                "m_stage": clinical_row["m_stage"],
                "overall_survival_months": clinical_row["overall_survival_months"],
                "overall_survival_event": clinical_row["overall_survival_event"],
                "recurrence_free_survival_months": clinical_row["recurrence_free_survival_months"],
                "recurrence_event": clinical_row["recurrence_event"],
                "adjuvant_chemotherapy": clinical_row["adjuvant_chemotherapy"],
                "adjuvant_chemotherapy_binary": clinical_row["adjuvant_chemotherapy_binary"],
                "matrix_column_index": str(position),
                "link_method": link_method,
            }
        )

    linked_patient_ids = {row["patient_id"] for row in link_rows}
    unmatched_matrix = [
        matrix_patient_id
        for matrix_patient_id in matrix_patient_ids
        if patient_normalizer(matrix_patient_id) not in linked_patient_ids
    ]
    unmatched_clinical = [patient_id for patient_id in clinical_rows if patient_id not in linked_patient_ids]
    summary_row = {
        "dataset_id": dataset_id,
        "matrix_sample_count": str(len(sample_accessions)),
        "clinical_sample_count": str(len(clinical_rows)),
        "linked_sample_count": str(len(link_rows)),
        "unmatched_matrix_samples": str(len(unmatched_matrix)),
        "unmatched_clinical_samples": str(len(unmatched_clinical)),
        "os_nonmissing_rows": str(sum(1 for row in link_rows if row["overall_survival_months"])),
        "rfs_nonmissing_rows": str(sum(1 for row in link_rows if row["recurrence_free_survival_months"])),
        "adjuvant_nonmissing_rows": str(sum(1 for row in link_rows if row["adjuvant_chemotherapy_binary"])),
    }
    return link_rows, summary_row, selected_positions


def write_filtered_matrix(matrix_input: Path, matrix_output: Path, selected_positions: list[int]) -> None:
    matrix_output.parent.mkdir(parents=True, exist_ok=True)
    selected_indices = [0] + selected_positions
    expected_row_length = max(selected_positions, default=0) + 1

    with gzip.open(matrix_input, "rt", encoding="utf-8", newline="") as source, gzip.open(
        matrix_output, "wt", encoding="utf-8", newline=""
    ) as target:
        writer = csv.writer(target, delimiter="\t", lineterminator="\n")
        in_table = False
        for raw_line in source:
            line = raw_line.rstrip("\n")
            if line == "!series_matrix_table_begin":
                in_table = True
                continue
            if line == "!series_matrix_table_end":
                break
            if not in_table:
                continue

            fields = parse_matrix_fields(line)
            if len(fields) < expected_row_length:
                fields.extend([""] * (expected_row_length - len(fields)))
            writer.writerow([fields[index] for index in selected_indices])
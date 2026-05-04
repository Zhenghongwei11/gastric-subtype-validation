from __future__ import annotations

import argparse
import csv
import gzip
from pathlib import Path

import xlrd


DATASET_ID = "GSE62254"
SHEET_NAME = "FINAL"
SUBTYPE_LABELS = {
    "0": "MSS/TP53-",
    "1": "MSS/TP53+",
    "2": "MSI",
    "3": "EMT",
}


def normalize_value(value: object) -> str:
    if value in (None, ""):
        return ""
    if isinstance(value, float):
        if value.is_integer():
            return str(int(value))
        return format(value, ".10g")
    return str(value).strip()


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def parse_matrix_fields(line: str) -> list[str]:
    return next(csv.reader([line], delimiter="\t", quotechar='"'))


def load_matrix_metadata(matrix_input: Path) -> tuple[list[str], list[str], list[str]]:
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
        raise SystemExit("Failed to parse required sample metadata from GSE62254 matrix")
    if not (len(sample_titles) == len(sample_accessions) == len(patient_ids)):
        raise SystemExit("GSE62254 matrix metadata lengths do not match")
    return sample_titles, sample_accessions, patient_ids


def load_clinical_rows(clinical_input: Path) -> dict[str, dict[str, str]]:
    workbook = xlrd.open_workbook(str(clinical_input))
    sheet = workbook.sheet_by_name(SHEET_NAME)
    headers = [normalize_value(value) for value in sheet.row_values(0)]
    column_index = {name: index for index, name in enumerate(headers)}

    rows: dict[str, dict[str, str]] = {}
    for row_index in range(1, sheet.nrows):
        values = [normalize_value(value) for value in sheet.row_values(row_index)]
        patient_id = values[column_index["Tumor ID"]]
        if not patient_id:
            continue
        subtype_code = values[column_index["Mol. Subtype: 0=MSS/TP53-, 1=MSS/TP53+, 2 = MSI, 3= EMT"]]
        rows[patient_id] = {
            "patient_id": patient_id,
            "clinical_sample_name": values[column_index["Sample\nName"]],
            "sex": values[column_index["sex"]],
            "age": values[column_index["age"]],
            "pathology": values[column_index["Pathology"]],
            "stage_tnm": values[column_index["stage(TNM)"]],
            "t_stage": values[column_index["T"]],
            "n_stage": values[column_index["N"]],
            "m_stage": values[column_index["M"]],
            "pathologic_stage": values[column_index["pStage"]],
            "lauren_code": values[column_index["LAUREN 1=intestinal, 2=diffuse, 3=mixed"]],
            "lauren": values[column_index["Lauren"]],
            "documented_recurrence": values[column_index["documented recurrence no=0  yes=1 unknown=2"]],
            "follow_up_status": values[column_index["FU status0=alive without ds, 1=alive with recurren ds, 2=dead without ds, 3=dead d/t recurrent ds, 4=dead, unknown, 5= FU loss"]],
            "date_of_death_or_last_follow_up": values[column_index["date of death or last follow up"]],
            "disease_free_survival_months": values[column_index["DFS\n(months)"]],
            "overall_survival_months": values[column_index["OS\n(months)"]],
            "molecular_subtype_code": subtype_code,
            "molecular_subtype_label": SUBTYPE_LABELS.get(subtype_code, "unknown"),
        }
    return rows


def build_link_rows(
    sample_titles: list[str],
    sample_accessions: list[str],
    patient_ids: list[str],
    clinical_rows: dict[str, dict[str, str]],
) -> tuple[list[dict[str, str]], dict[str, str], list[int]]:
    link_rows: list[dict[str, str]] = []
    selected_positions: list[int] = []

    for position, (gsm_id, sample_title, patient_id) in enumerate(
        zip(sample_accessions, sample_titles, patient_ids),
        start=1,
    ):
        clinical_row = clinical_rows.get(patient_id)
        if clinical_row is None:
            continue

        selected_positions.append(position)
        link_rows.append(
            {
                "dataset_id": DATASET_ID,
                "gsm_id": gsm_id,
                "patient_id": patient_id,
                "matrix_sample_title": sample_title,
                "clinical_sample_name": clinical_row["clinical_sample_name"],
                "sex": clinical_row["sex"],
                "age": clinical_row["age"],
                "pathology": clinical_row["pathology"],
                "stage_tnm": clinical_row["stage_tnm"],
                "t_stage": clinical_row["t_stage"],
                "n_stage": clinical_row["n_stage"],
                "m_stage": clinical_row["m_stage"],
                "pathologic_stage": clinical_row["pathologic_stage"],
                "lauren_code": clinical_row["lauren_code"],
                "lauren": clinical_row["lauren"],
                "documented_recurrence": clinical_row["documented_recurrence"],
                "follow_up_status": clinical_row["follow_up_status"],
                "date_of_death_or_last_follow_up": clinical_row["date_of_death_or_last_follow_up"],
                "disease_free_survival_months": clinical_row["disease_free_survival_months"],
                "overall_survival_months": clinical_row["overall_survival_months"],
                "molecular_subtype_code": clinical_row["molecular_subtype_code"],
                "molecular_subtype_label": clinical_row["molecular_subtype_label"],
                "matrix_column_index": str(position),
                "link_method": "patient_id_exact",
            }
        )

    linked_patient_ids = {row["patient_id"] for row in link_rows}
    unmatched_matrix = [patient_id for patient_id in patient_ids if patient_id not in linked_patient_ids]
    unmatched_clinical = [patient_id for patient_id in clinical_rows if patient_id not in linked_patient_ids]
    summary_row = {
        "dataset_id": DATASET_ID,
        "matrix_sample_count": str(len(sample_accessions)),
        "clinical_sample_count": str(len(clinical_rows)),
        "linked_sample_count": str(len(link_rows)),
        "unmatched_matrix_samples": str(len(unmatched_matrix)),
        "unmatched_clinical_samples": str(len(unmatched_clinical)),
        "os_nonmissing_rows": str(sum(1 for row in link_rows if row["overall_survival_months"])),
        "dfs_nonmissing_rows": str(sum(1 for row in link_rows if row["disease_free_survival_months"])),
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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Link GSE62254 series matrix to ACRG clinical supplement")
    parser.add_argument("--matrix-input", required=True)
    parser.add_argument("--clinical-input", required=True)
    parser.add_argument("--link-output", required=True)
    parser.add_argument("--summary-output", required=True)
    parser.add_argument("--matrix-output", required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    sample_titles, sample_accessions, patient_ids = load_matrix_metadata(Path(args.matrix_input))
    clinical_rows = load_clinical_rows(Path(args.clinical_input))
    link_rows, summary_row, selected_positions = build_link_rows(
        sample_titles,
        sample_accessions,
        patient_ids,
        clinical_rows,
    )
    if len(link_rows) != len(sample_accessions):
        raise SystemExit(
            f"Only linked {len(link_rows)} of {len(sample_accessions)} GSE62254 matrix samples"
        )

    write_tsv(Path(args.link_output), link_rows)
    write_tsv(Path(args.summary_output), [summary_row])
    write_filtered_matrix(Path(args.matrix_input), Path(args.matrix_output), selected_positions)


if __name__ == "__main__":
    main()
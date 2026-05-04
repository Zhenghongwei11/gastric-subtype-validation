from __future__ import annotations

import argparse
import csv
import gzip
from pathlib import Path


DATASET_ID = "GSE84437"
REQUIRED_FIELDS = [
    "age",
    "sex",
    "ptstage",
    "pnstage",
    "death",
    "duration overall survival",
]


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def parse_matrix_fields(line: str) -> list[str]:
    return next(csv.reader([line], delimiter="\t", quotechar='"'))


def parse_characteristic(value: str) -> tuple[str, str]:
    cleaned = value.strip()
    if ": " in cleaned:
        field_name, field_value = cleaned.split(": ", 1)
    elif ":" in cleaned:
        field_name, field_value = cleaned.split(":", 1)
    else:
        field_name, field_value = cleaned, ""
    return field_name.strip().lower(), field_value.strip()


def load_matrix_annotations(matrix_input: Path) -> tuple[list[str], list[str], list[dict[str, str]]]:
    sample_titles: list[str] = []
    sample_accessions: list[str] = []
    sample_annotations: list[dict[str, str]] = []

    with gzip.open(matrix_input, "rt", encoding="utf-8", errors="replace", newline="") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if line.startswith("!Sample_title"):
                sample_titles = parse_matrix_fields(line)[1:]
            elif line.startswith("!Sample_geo_accession"):
                sample_accessions = parse_matrix_fields(line)[1:]
                sample_annotations = [{} for _ in sample_accessions]
            elif line.startswith("!Sample_characteristics_ch"):
                values = parse_matrix_fields(line)[1:]
                for index, value in enumerate(values):
                    if index >= len(sample_annotations):
                        continue
                    field_name, field_value = parse_characteristic(value)
                    if field_name:
                        sample_annotations[index][field_name] = field_value
            elif line == "!series_matrix_table_begin":
                break

    if not sample_titles or not sample_accessions or not sample_annotations:
        raise SystemExit("Failed to parse sample metadata from GSE84437 series matrix")
    if not (len(sample_titles) == len(sample_accessions) == len(sample_annotations)):
        raise SystemExit("Sample metadata lengths do not match in GSE84437 matrix")
    return sample_titles, sample_accessions, sample_annotations


def build_link_rows(
    sample_titles: list[str],
    sample_accessions: list[str],
    sample_annotations: list[dict[str, str]],
) -> tuple[list[dict[str, str]], dict[str, str], list[int]]:
    link_rows: list[dict[str, str]] = []
    selected_positions: list[int] = []

    for position, (gsm_id, sample_title, annotations) in enumerate(
        zip(sample_accessions, sample_titles, sample_annotations),
        start=1,
    ):
        if any(not annotations.get(field, "") for field in REQUIRED_FIELDS):
            continue

        selected_positions.append(position)
        link_rows.append(
            {
                "dataset_id": DATASET_ID,
                "gsm_id": gsm_id,
                "sample_title": sample_title,
                "age": annotations["age"],
                "sex": annotations["sex"],
                "pt_stage": annotations["ptstage"],
                "pn_stage": annotations["pnstage"],
                "overall_survival_months": annotations["duration overall survival"],
                "overall_survival_event": annotations["death"],
                "matrix_column_index": str(position),
                "link_method": "matrix_annotation_direct",
            }
        )

    summary_row = {
        "dataset_id": DATASET_ID,
        "matrix_sample_count": str(len(sample_accessions)),
        "linked_sample_count": str(len(link_rows)),
        "unmatched_matrix_samples": str(len(sample_accessions) - len(link_rows)),
        "age_nonmissing_rows": str(sum(1 for row in link_rows if row["age"])),
        "ptstage_nonmissing_rows": str(sum(1 for row in link_rows if row["pt_stage"])),
        "pnstage_nonmissing_rows": str(sum(1 for row in link_rows if row["pn_stage"])),
        "os_nonmissing_rows": str(sum(1 for row in link_rows if row["overall_survival_months"])),
        "event_nonmissing_rows": str(sum(1 for row in link_rows if row["overall_survival_event"])),
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
    parser = argparse.ArgumentParser(description="Extract survival-ready GSE84437 rows from matrix annotations")
    parser.add_argument("--matrix-input", required=True)
    parser.add_argument("--link-output", required=True)
    parser.add_argument("--summary-output", required=True)
    parser.add_argument("--matrix-output", required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    sample_titles, sample_accessions, sample_annotations = load_matrix_annotations(Path(args.matrix_input))
    link_rows, summary_row, selected_positions = build_link_rows(
        sample_titles,
        sample_accessions,
        sample_annotations,
    )
    if not link_rows:
        raise SystemExit("No survival-ready rows found in GSE84437 series matrix")

    write_tsv(Path(args.link_output), link_rows)
    write_tsv(Path(args.summary_output), [summary_row])
    write_filtered_matrix(Path(args.matrix_input), Path(args.matrix_output), selected_positions)


if __name__ == "__main__":
    main()
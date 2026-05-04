from __future__ import annotations

import argparse
import csv
import gzip
from pathlib import Path


DATASET_ID = "GSE15459"


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


def normalize_sample_title(value: str) -> str:
    return value.replace(" [EXCLUDED]", "").strip()


def load_matrix_sample_metadata(matrix_input: Path) -> tuple[list[str], list[str]]:
    sample_titles: list[str] = []
    sample_accessions: list[str] = []

    with gzip.open(matrix_input, "rt", encoding="utf-8", newline="") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if line.startswith("!Sample_title"):
                sample_titles = [normalize_sample_title(value) for value in parse_matrix_fields(line)[1:]]
            elif line.startswith("!Sample_geo_accession"):
                sample_accessions = parse_matrix_fields(line)[1:]
            elif line == "!series_matrix_table_begin":
                break

    if not sample_titles or not sample_accessions:
        raise SystemExit("Failed to parse sample metadata from GSE15459 series matrix")
    if len(sample_titles) != len(sample_accessions):
        raise SystemExit("Sample title count does not match sample accession count in GSE15459 matrix")
    return sample_titles, sample_accessions


def build_link_rows(
    clinical_rows: list[dict[str, str]],
    sample_titles: list[str],
    sample_accessions: list[str],
) -> tuple[list[dict[str, str]], dict[str, str], list[int]]:
    clinical_by_gsm = {row["gsm_id"]: row for row in clinical_rows}
    link_rows: list[dict[str, str]] = []
    selected_positions: list[int] = []

    for position, (gsm_id, sample_title) in enumerate(zip(sample_accessions, sample_titles), start=1):
        clinical_row = clinical_by_gsm.get(gsm_id)
        if clinical_row is None:
            continue

        selected_positions.append(position)
        link_rows.append(
            {
                "dataset_id": DATASET_ID,
                "gsm_id": gsm_id,
                "sample_title": sample_title,
                "cel_file": clinical_row["cel_file"],
                "patient_id": clinical_row["patient_id"],
                "subtype": clinical_row["subtype"],
                "stage": clinical_row["stage"],
                "matrix_column_index": str(position),
                "link_method": "gsm_exact",
            }
        )

    linked_gsm_ids = {row["gsm_id"] for row in link_rows}
    unmatched_clinical = [row["gsm_id"] for row in clinical_rows if row["gsm_id"] not in linked_gsm_ids]

    summary_row = {
        "dataset_id": DATASET_ID,
        "matrix_sample_count": str(len(sample_accessions)),
        "clinical_sample_count": str(len(clinical_rows)),
        "linked_sample_count": str(len(link_rows)),
        "unmatched_matrix_samples": str(len(sample_accessions) - len(link_rows)),
        "unmatched_clinical_samples": str(len(unmatched_clinical)),
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
    parser = argparse.ArgumentParser(description="Link GSE15459 expression matrix columns to clinical rows")
    parser.add_argument("--clinical-input", required=True)
    parser.add_argument("--matrix-input", required=True)
    parser.add_argument("--link-output", required=True)
    parser.add_argument("--summary-output", required=True)
    parser.add_argument("--matrix-output", required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    clinical_rows = read_tsv(Path(args.clinical_input))
    if not clinical_rows:
        raise SystemExit("No clinical rows available for GSE15459 expression linking")

    sample_titles, sample_accessions = load_matrix_sample_metadata(Path(args.matrix_input))
    link_rows, summary_row, selected_positions = build_link_rows(clinical_rows, sample_titles, sample_accessions)
    if len(link_rows) != len(clinical_rows):
        raise SystemExit(
            f"Only linked {len(link_rows)} of {len(clinical_rows)} clinical rows from GSE15459 matrix"
        )

    write_tsv(Path(args.link_output), link_rows)
    write_tsv(Path(args.summary_output), [summary_row])
    write_filtered_matrix(Path(args.matrix_input), Path(args.matrix_output), selected_positions)


if __name__ == "__main__":
    main()
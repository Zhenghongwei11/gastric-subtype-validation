from __future__ import annotations

import argparse
from pathlib import Path

from natcom_c2_link_utils import (
    build_link_rows,
    load_clinical_rows,
    load_matrix_metadata,
    write_filtered_matrix,
    write_tsv,
)


DATASET_ID = "GSE26899"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Link GSE26899 series matrix to standardized Nat Commun clinical rows")
    parser.add_argument("--clinical-input", required=True)
    parser.add_argument("--matrix-input", required=True)
    parser.add_argument("--link-output", required=True)
    parser.add_argument("--summary-output", required=True)
    parser.add_argument("--matrix-output", required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    clinical_rows = load_clinical_rows(Path(args.clinical_input), DATASET_ID)
    sample_titles, sample_accessions, patient_ids = load_matrix_metadata(Path(args.matrix_input), DATASET_ID)
    link_rows, summary_row, selected_positions = build_link_rows(
        DATASET_ID,
        sample_titles,
        sample_accessions,
        patient_ids,
        clinical_rows,
        patient_normalizer=lambda value: value.strip(),
        link_method="patient_id_exact",
    )
    if len(link_rows) != len(clinical_rows):
        raise SystemExit(f"Only linked {len(link_rows)} of {len(clinical_rows)} clinical rows from {DATASET_ID}")

    write_tsv(Path(args.link_output), link_rows)
    write_tsv(Path(args.summary_output), [summary_row])
    write_filtered_matrix(Path(args.matrix_input), Path(args.matrix_output), selected_positions)


if __name__ == "__main__":
    main()
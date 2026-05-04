from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path


def parse_sample_count(text: str) -> str:
    patterns = [r"Samples\s*\(([^)]+)\)", r"Cases:\s*([0-9]+)"]
    for pattern in patterns:
        match = re.search(pattern, text, flags=re.IGNORECASE)
        if match:
            return match.group(1).strip()
    return "tbd"


def parse_platform(text: str) -> str:
    match = re.search(r"Platforms?\s*\([^)]*\):\s*(.+)", text, flags=re.IGNORECASE)
    if match:
        return match.group(1).strip()
    return "tbd"


def parse_title(text: str) -> str:
    for prefix in ["Title:", "Project Name:", "Project ID:"]:
        match = re.search(rf"{re.escape(prefix)}\s*(.+)", text)
        if match:
            return match.group(1).strip()
    return "tbd"


def parse_endpoint_hint(text: str) -> str:
    lower_text = text.lower()
    hints = []
    if "survival" in lower_text:
        hints.append("survival")
    if "recurrence" in lower_text:
        hints.append("recurrence")
    if "outcome" in lower_text:
        hints.append("outcome_file")
    return ";".join(hints) if hints else "tbd"


def parse_treatment_hint(text: str) -> str:
    lower_text = text.lower()
    hints = []
    for token in ["chemotherapy", "chemoradiotherapy", "treatment", "adjuvant"]:
        if token in lower_text:
            hints.append(token)
    return ";".join(dict.fromkeys(hints)) if hints else "none_detected"


def parse_outcome_file(text: str) -> str:
    matches = re.findall(r"Supplementary file:\s*(.+)", text, flags=re.IGNORECASE)
    for candidate in matches:
        file_name = candidate.strip()
        lower_name = file_name.lower()
        if any(token in lower_name for token in ["outcome", "survival", "clinical", "endpoint", "response"]):
            return file_name
    return "none"


def parse_summary_row(path: Path) -> dict[str, str]:
    text = path.read_text(encoding="utf-8")
    dataset_id = path.stem
    return {
        "dataset_id": dataset_id,
        "title": parse_title(text),
        "sample_count": parse_sample_count(text),
        "platform": parse_platform(text),
        "endpoint_hint": parse_endpoint_hint(text),
        "treatment_hint": parse_treatment_hint(text),
        "outcome_file_hint": parse_outcome_file(text),
    }


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Parse cached cohort accession metadata summaries")
    parser.add_argument("--input-dir", required=True)
    parser.add_argument("--output", required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    input_dir = Path(args.input_dir)
    rows = [parse_summary_row(path) for path in sorted(input_dir.glob("*.txt"))]
    if not rows:
        raise SystemExit(f"No accession summary files found in {input_dir}")
    write_tsv(Path(args.output), rows)


if __name__ == "__main__":
    main()
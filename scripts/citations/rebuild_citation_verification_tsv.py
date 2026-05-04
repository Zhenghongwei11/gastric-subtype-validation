#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path


REF_LINE_RE = re.compile(r"^(\d+)\.\s+(.*)$")
DOI_RE = re.compile(r"\bdoi:\s*([^\s]+)", re.IGNORECASE)


def extract_references(manuscript_text: str) -> list[tuple[int, str]]:
    if "## References" not in manuscript_text:
        raise SystemExit("No '## References' section found.")
    _, refs = manuscript_text.split("## References", 1)
    items: list[tuple[int, str]] = []
    for raw in refs.splitlines():
        m = REF_LINE_RE.match(raw.strip())
        if not m:
            continue
        items.append((int(m.group(1)), m.group(2).strip()))
    if not items:
        raise SystemExit("No numbered reference lines found under '## References'.")
    return items


def parse_doi(reference: str) -> str:
    m = DOI_RE.search(reference)
    if not m:
        return ""
    doi = m.group(1).rstrip(".").strip()
    # Strip trailing punctuation that sometimes sneaks in.
    doi = doi.rstrip(").,;")
    return doi


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description="Rebuild a CITATION_VERIFICATION.tsv from the manuscript reference list.")
    ap.add_argument("--manuscript", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--notes", default="TODO: verify (2 sources + full text)")
    return ap.parse_args()


def main() -> int:
    args = parse_args()
    manuscript_path = Path(args.manuscript)
    out_path = Path(args.out)
    text = manuscript_path.read_text(encoding="utf-8", errors="replace")
    refs = extract_references(text)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["number", "reference", "doi", "pmid", "pmcid", "notes"],
            delimiter="\t",
        )
        writer.writeheader()
        for number, reference in refs:
            writer.writerow(
                {
                    "number": str(number),
                    "reference": reference,
                    "doi": parse_doi(reference),
                    "pmid": "",
                    "pmcid": "",
                    "notes": args.notes,
                }
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


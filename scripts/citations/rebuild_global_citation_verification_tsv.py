#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import re
from dataclasses import dataclass
from pathlib import Path


BRACKET_CITATION_RE = re.compile(r"\[(\d[\d,\s-]*)\]")
REF_LINE_RE = re.compile(r"^(\d+)\.\s+(.*)$")


def expand_citation_blob(blob: str) -> list[int]:
    nums: list[int] = []
    parts = [p.strip() for p in blob.strip().split(",") if p.strip()]
    for part in parts:
        if "-" in part:
            a, b = part.split("-", 1)
            try:
                start = int(a.strip())
                end = int(b.strip())
            except ValueError:
                continue
            if start > end:
                start, end = end, start
            nums.extend(range(start, end + 1))
        else:
            try:
                nums.append(int(part))
            except ValueError:
                continue
    return [n for n in nums if n > 0]


def extract_references(text: str) -> list[tuple[int, str]]:
    if "## References" not in text:
        raise SystemExit("No '## References' section found.")
    _, refs = text.split("## References", 1)
    items: list[tuple[int, str]] = []
    for raw in refs.splitlines():
        m = REF_LINE_RE.match(raw.strip())
        if not m:
            continue
        items.append((int(m.group(1)), m.group(2).strip()))
    if not items:
        raise SystemExit("No numbered reference lines found under '## References'.")
    return items


def build_first_appearance_sections(text_before_refs: str) -> dict[int, str]:
    current_section = "Unknown"
    first_section: dict[int, str] = {}
    for raw in text_before_refs.splitlines():
        line = raw.strip()
        if line.startswith("#"):
            # Keep a compact section label; prefer the deepest heading.
            current_section = line.lstrip("#").strip() or current_section
            continue
        for m in BRACKET_CITATION_RE.finditer(line):
            for n in expand_citation_blob(m.group(1)):
                if n not in first_section:
                    first_section[n] = current_section
    return first_section


def short_citation_from_reference(reference: str) -> str:
    # Heuristic: "FirstAuthor et al. Journal YEAR"
    first_author = reference.split(".")[0].split(",")[0].strip()
    year_match = re.search(r"\b(19|20)\d{2}\b", reference)
    year = year_match.group(0) if year_match else ""
    journal = ""
    # Try to capture "Journal. YEAR;..." or "Journal. YEAR;" patterns
    m = re.search(r"\.\s*([^.\n]+)\.\s*(19|20)\d{2}\s*;", reference)
    if m:
        journal = m.group(1).strip()
    parts = [p for p in [first_author + " et al." if "et al" in reference else first_author, journal, year] if p]
    return " ".join(parts).strip()


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description="Rebuild docs/CITATION_VERIFICATION.tsv from a manuscript reference list.")
    ap.add_argument("--manuscript", required=True)
    ap.add_argument("--out", required=True)
    return ap.parse_args()


def main() -> int:
    args = parse_args()
    text = Path(args.manuscript).read_text(encoding="utf-8", errors="replace")
    before, _ = text.split("## References", 1)
    ref_items = extract_references(text)
    first_section = build_first_appearance_sections(before)

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "reference_number",
                "short_citation",
                "first_appearance_section",
                "verification_source_1",
                "verification_source_2",
                "full_text_check_status",
                "bibliographic_fields_checked",
                "verification_status",
                "notes",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for number, reference in ref_items:
            writer.writerow(
                {
                    "reference_number": str(number),
                    "short_citation": short_citation_from_reference(reference),
                    "first_appearance_section": first_section.get(number, "Unknown"),
                    "verification_source_1": "pending",
                    "verification_source_2": "pending",
                    "full_text_check_status": "pending",
                    "bibliographic_fields_checked": "authors,title,journal,year,volume,issue,pages,doi",
                    "verification_status": "draft_initialized",
                    "notes": "Current manuscript entry captured; external two-source verification and full-text check still required before final submission package lock.",
                }
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


#!/usr/bin/env python3
from __future__ import annotations

import argparse
import re
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class ReferenceItem:
    number: int
    body: str


BRACKET_CITATION_RE = re.compile(r"\[(\d[\d,\s-]*)\]")
REF_LINE_RE = re.compile(r"^(\d+)\.\s+(.*)$")


def expand_citation_blob(blob: str) -> list[int]:
    nums: list[int] = []
    blob = blob.strip()
    if not blob:
        return nums
    parts = [p.strip() for p in blob.split(",")]
    for part in parts:
        if not part:
            continue
        if "-" in part:
            a, b = part.split("-", 1)
            try:
                start = int(a.strip())
                end = int(b.strip())
            except ValueError:
                continue
            if start > end:
                start, end = end, start
            nums.extend(list(range(start, end + 1)))
        else:
            try:
                nums.append(int(part))
            except ValueError:
                continue
    return nums


def compress_numbers(nums: list[int]) -> str:
    if not nums:
        return ""
    nums = sorted(dict.fromkeys(nums))
    ranges: list[str] = []
    start = prev = nums[0]
    for n in nums[1:]:
        if n == prev + 1:
            prev = n
            continue
        ranges.append(f"{start}-{prev}" if start != prev else f"{start}")
        start = prev = n
    ranges.append(f"{start}-{prev}" if start != prev else f"{start}")
    return ", ".join(ranges)


def parse_reference_items(ref_text: str) -> list[ReferenceItem]:
    items: list[ReferenceItem] = []
    for raw in ref_text.splitlines():
        m = REF_LINE_RE.match(raw.strip())
        if not m:
            continue
        items.append(ReferenceItem(number=int(m.group(1)), body=m.group(2).strip()))
    return items


def build_first_appearance_order(text: str) -> list[int]:
    order: list[int] = []
    seen: set[int] = set()
    for m in BRACKET_CITATION_RE.finditer(text):
        nums = expand_citation_blob(m.group(1))
        for n in nums:
            if n <= 0:
                continue
            if n not in seen:
                seen.add(n)
                order.append(n)
    return order


def renumber_citations(text: str, mapping: dict[int, int]) -> str:
    def repl(match: re.Match[str]) -> str:
        nums = expand_citation_blob(match.group(1))
        mapped = [mapping.get(n, n) for n in nums if n > 0]
        return "[" + compress_numbers(mapped) + "]"

    return BRACKET_CITATION_RE.sub(repl, text)


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description="Renumber Vancouver bracket citations and reference list by first appearance.")
    ap.add_argument("--manuscript", required=True, help="Markdown manuscript containing a '## References' section.")
    ap.add_argument("--in-place", action="store_true", help="Rewrite the manuscript in place.")
    ap.add_argument("--out", help="Output path (default: stdout) when not using --in-place.")
    ap.add_argument("--require-all-cited", action="store_true", help="Fail if any reference list entry is uncited.")
    return ap.parse_args()


def main() -> int:
    args = parse_args()
    path = Path(args.manuscript)
    text = path.read_text(encoding="utf-8", errors="replace")
    if "## References" not in text:
        raise SystemExit("No '## References' section found.")

    pre, refs = text.split("## References", 1)
    pre = pre.rstrip() + "\n\n## References\n\n"

    ref_items = parse_reference_items(refs)
    if not ref_items:
        raise SystemExit("No numbered references found under '## References'.")

    first_order = build_first_appearance_order(pre)
    if not first_order:
        raise SystemExit("No bracketed citations like [1] found in manuscript text.")

    mapping: dict[int, int] = {}
    for new_num, old_num in enumerate(first_order, start=1):
        mapping[old_num] = new_num

    ref_by_old = {item.number: item for item in ref_items}
    missing_refs = [n for n in mapping if n not in ref_by_old]
    if missing_refs:
        raise SystemExit(f"Citations refer to missing reference numbers: {missing_refs[:20]}")

    uncited = sorted([item.number for item in ref_items if item.number not in mapping])
    if args.require_all_cited and uncited:
        raise SystemExit(f"Uncited references present: {uncited[:30]}")

    # Rebuild references in new order; append any uncited refs at the end (stable).
    new_refs: list[ReferenceItem] = []
    for old_num in first_order:
        item = ref_by_old[old_num]
        new_refs.append(ReferenceItem(number=mapping[old_num], body=item.body))
    tail_start = len(new_refs) + 1
    for i, old_num in enumerate(uncited, start=tail_start):
        item = ref_by_old[old_num]
        new_refs.append(ReferenceItem(number=i, body=item.body))

    rewritten_pre = renumber_citations(pre, mapping)
    rewritten_refs = "\n".join(f"{item.number}. {item.body}" for item in new_refs).rstrip() + "\n"
    rewritten = rewritten_pre + rewritten_refs

    if args.in_place:
        path.write_text(rewritten, encoding="utf-8")
        return 0

    if args.out:
        Path(args.out).write_text(rewritten, encoding="utf-8")
        return 0

    print(rewritten, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


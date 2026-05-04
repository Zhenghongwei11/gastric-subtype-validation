from __future__ import annotations

import argparse
import csv
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def format_effect(row: dict[str, str]) -> str:
    return f"HR {row['effect']} ({row['ci_lower']}-{row['ci_upper']}), P={row['pvalue']}"


def build_count_string(summary_by_key: dict[tuple[str, str, str], dict[str, str]], dataset_id: str, adjuvant: str) -> str:
    ep = summary_by_key[(dataset_id, "EP", adjuvant)]
    mp = summary_by_key[(dataset_id, "MP", adjuvant)]
    return (
        f"EP n={ep['n']}, MP n={mp['n']}; "
        f"OS events {ep['overall_survival_events']} vs {mp['overall_survival_events']}; "
        f"RFS events {ep['recurrence_free_survival_events']} vs {mp['recurrence_free_survival_events']}"
    )


def classify_takeaway(treated_row: dict[str, str], untreated_row: dict[str, str], interaction_row: dict[str, str]) -> str:
    treated_effect = float(treated_row["effect"])
    untreated_effect = float(untreated_row["effect"])
    interaction_effect = float(interaction_row["effect"])
    interaction_pvalue = float(interaction_row["pvalue"])

    if treated_effect > 1.5 and untreated_effect > 1.5 and 0.8 <= interaction_effect <= 1.5 and interaction_pvalue >= 0.2:
        lead = "similar adverse MP-versus-EP hazards in both adjuvant strata"
    elif abs(treated_effect - untreated_effect) < 0.5:
        lead = "similar adverse MP-versus-EP hazards in both adjuvant strata"
    elif treated_effect > untreated_effect:
        lead = "stronger MP-versus-EP hazard under adjuvant exposure"
    else:
        lead = "similar or stronger MP-versus-EP hazard without adjuvant exposure"

    if interaction_pvalue < 0.05:
        tail = "interaction supported at nominal significance"
    elif 0.05 <= interaction_pvalue < 0.2:
        tail = "interaction remained imprecise"
    else:
        tail = "near-null interaction estimate"
        if interaction_effect > 1.25:
            tail = "interaction remained imprecise"

    return f"{lead}; {tail}."


def build_rows(effect_rows: list[dict[str, str]], summary_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    summary_by_key = {
        (row["dataset_id"], row["subgroup"], row["adjuvant_chemotherapy_binary"]): row for row in summary_rows
    }
    by_key = {(row["dataset_id"], row["outcome"], row["model"]): row for row in effect_rows}

    rows: list[dict[str, str]] = []
    dataset_outcomes = sorted({(row["dataset_id"], row["outcome"]) for row in effect_rows})
    for dataset_id, outcome in dataset_outcomes:
        treated_row = by_key[(dataset_id, outcome, "cox_subgroup_in_adjuvant_treated")]
        untreated_row = by_key[(dataset_id, outcome, "cox_subgroup_in_adjuvant_untreated")]
        interaction_row = by_key[(dataset_id, outcome, "cox_subgroup_adjuvant_interaction")]
        rows.append(
            {
                "dataset_id": dataset_id,
                "outcome": outcome,
                "treated_counts": build_count_string(summary_by_key, dataset_id, "1"),
                "untreated_counts": build_count_string(summary_by_key, dataset_id, "0"),
                "treated_effect": format_effect(treated_row),
                "untreated_effect": format_effect(untreated_row),
                "interaction_effect": format_effect(interaction_row),
                "cohort_takeaway": classify_takeaway(treated_row, untreated_row, interaction_row),
            }
        )
    return rows


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build cohort-wise Nat Commun C2 supplement interpretation table")
    parser.add_argument("--effect-input", required=True)
    parser.add_argument("--summary-input", required=True)
    parser.add_argument("--output", required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    effect_rows = read_tsv(Path(args.effect_input))
    summary_rows = read_tsv(Path(args.summary_input))
    write_tsv(Path(args.output), build_rows(effect_rows, summary_rows))


if __name__ == "__main__":
    main()
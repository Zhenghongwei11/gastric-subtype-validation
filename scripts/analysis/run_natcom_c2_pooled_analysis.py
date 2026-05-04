from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import pandas as pd
from lifelines import CoxPHFitter

from run_natcom_c2_survival import (
    DATASET_ORDER,
    OUTCOME_SPECS,
    SUBGROUP_COMPARISON,
    SUBGROUP_REFERENCE,
    fit_interaction_model,
    format_number,
    load_frame,
)


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def load_pooled_frame(gse26899_input: Path, gse26901_input: Path) -> pd.DataFrame:
    frame = pd.concat(
        [
            load_frame(gse26899_input, "GSE26899"),
            load_frame(gse26901_input, "GSE26901"),
        ],
        ignore_index=True,
    )
    frame["dataset_binary"] = (frame["dataset_id"] == "GSE26901").astype(int)
    return frame


def fit_pooled_interaction_model(frame: pd.DataFrame, outcome: str, duration_col: str, event_col: str) -> dict[str, str]:
    model_frame = frame.dropna(subset=[duration_col, event_col]).copy()
    model_frame["subgroup_adjuvant_interaction"] = model_frame["subgroup_binary"] * model_frame["adjuvant_binary"]

    fitter = CoxPHFitter(penalizer=0.05)
    fitter.fit(
        model_frame[
            [duration_col, event_col, "subgroup_binary", "adjuvant_binary", "subgroup_adjuvant_interaction", "dataset_binary"]
        ],
        duration_col=duration_col,
        event_col=event_col,
    )
    summary = fitter.summary.loc["subgroup_adjuvant_interaction"]
    return {
        "claim_id": "C2",
        "dataset_id": "NatCommFamily",
        "outcome": outcome,
        "model": "cox_pooled_dataset_adjusted_interaction",
        "reference_subtype": SUBGROUP_REFERENCE,
        "comparison_subtype": SUBGROUP_COMPARISON,
        "effect_type": "interaction_hazard_ratio",
        "effect": format_number(summary["exp(coef)"]),
        "ci_lower": format_number(summary["exp(coef) lower 95%"]),
        "ci_upper": format_number(summary["exp(coef) upper 95%"]),
        "pvalue": format_number(summary["p"]),
        "fdr": "",
        "n": str(len(model_frame)),
        "events": str(int(model_frame[event_col].sum())),
        "covariates": "subgroup,adjuvant,subgroup_x_adjuvant,dataset",
        "adjuvant_stratum": "interaction",
    }


def log_hr_and_se_from_row(row: dict[str, str]) -> tuple[float, float]:
    effect = float(row["effect"])
    ci_lower = float(row["ci_lower"])
    ci_upper = float(row["ci_upper"])
    log_hr = math.log(effect)
    se = (math.log(ci_upper) - math.log(ci_lower)) / (2 * 1.96)
    return log_hr, se


def two_sided_normal_pvalue(z_value: float) -> float:
    return math.erfc(abs(z_value) / math.sqrt(2.0))


def chi_square_df1_survival_function(q_value: float) -> float:
    return math.erfc(math.sqrt(max(q_value, 0.0) / 2.0))


def build_meta_rows(gse26899_input: Path, gse26901_input: Path) -> list[dict[str, str]]:
    interaction_rows: list[dict[str, str]] = []
    for dataset_id, input_path in [("GSE26899", gse26899_input), ("GSE26901", gse26901_input)]:
        frame = load_frame(input_path, dataset_id)
        for outcome, duration_col, event_col in OUTCOME_SPECS:
            interaction_rows.append(fit_interaction_model(frame, dataset_id, outcome, duration_col, event_col))

    rows: list[dict[str, str]] = []
    for outcome in sorted({row["outcome"] for row in interaction_rows}):
        outcome_rows = [row for row in interaction_rows if row["outcome"] == outcome]
        log_hr_and_ses = [log_hr_and_se_from_row(row) for row in outcome_rows]
        weights = [1.0 / (se**2) for _, se in log_hr_and_ses]
        pooled_log_hr = sum(weight * log_hr for weight, (log_hr, _) in zip(weights, log_hr_and_ses)) / sum(weights)
        pooled_se = math.sqrt(1.0 / sum(weights))
        z_value = pooled_log_hr / pooled_se
        pooled_pvalue = two_sided_normal_pvalue(z_value)

        q_stat = sum(weight * ((log_hr - pooled_log_hr) ** 2) for weight, (log_hr, _) in zip(weights, log_hr_and_ses))
        degrees_of_freedom = max(len(outcome_rows) - 1, 1)
        heterogeneity_pvalue = chi_square_df1_survival_function(q_stat) if degrees_of_freedom == 1 else "nan"
        i_squared = max(0.0, (q_stat - degrees_of_freedom) / q_stat) * 100 if q_stat > 0 else 0.0

        rows.append(
            {
                "outcome": outcome,
                "meta_model": "fixed_effect_inverse_variance",
                "dataset_count": str(len(outcome_rows)),
                "pooled_hr": format_number(math.exp(pooled_log_hr)),
                "ci_lower": format_number(math.exp(pooled_log_hr - 1.96 * pooled_se)),
                "ci_upper": format_number(math.exp(pooled_log_hr + 1.96 * pooled_se)),
                "pooled_pvalue": format_number(pooled_pvalue),
                "cochran_q": format_number(q_stat),
                "heterogeneity_pvalue": format_number(float(heterogeneity_pvalue)) if heterogeneity_pvalue != "nan" else "nan",
                "i_squared": format_number(i_squared),
            }
        )
    return rows


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run pooled and fixed-effect meta Nat Commun C2 interaction analyses")
    parser.add_argument("--gse26899-input", required=True)
    parser.add_argument("--gse26901-input", required=True)
    parser.add_argument("--effect-output", required=True)
    parser.add_argument("--meta-output", required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    gse26899_input = Path(args.gse26899_input)
    gse26901_input = Path(args.gse26901_input)
    pooled_frame = load_pooled_frame(gse26899_input, gse26901_input)
    effect_rows = [
        fit_pooled_interaction_model(pooled_frame, outcome, duration_col, event_col)
        for outcome, duration_col, event_col in OUTCOME_SPECS
    ]
    write_tsv(Path(args.effect_output), effect_rows)
    write_tsv(Path(args.meta_output), build_meta_rows(gse26899_input, gse26901_input))


if __name__ == "__main__":
    main()
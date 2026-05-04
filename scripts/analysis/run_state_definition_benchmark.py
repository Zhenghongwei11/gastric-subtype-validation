from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import pandas as pd
from lifelines import CoxPHFitter

from run_acrg_cross_cohort_survival import parse_source_sheet
from run_gse84437_projection_survival import (
    REFERENCE_SUBTYPE,
    collapse_to_gene_matrix,
    load_expression_matrix,
    load_or_build_probe_gene_map,
    parse_stage_token,
    zscore_by_gene,
)


TOP_GENES_PER_DIRECTION = 25
COMPACT_METHOD = "compact_program_score"
BASELINE_METHOD = "published_subtype_mapping"


def format_number(value: float) -> str:
    return format(float(value), ".6g")


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def normalize_stage(raw_value: object) -> int | None:
    value = str(raw_value).strip()
    if not value:
        return None
    return parse_stage_token(value)


def normalize_stage_roman(raw_value: object) -> int | None:
    value = str(raw_value).strip().upper()
    if not value:
        return None
    roman_map = {"I": 1, "II": 2, "III": 3, "IV": 4}
    for token in ("IV", "III", "II", "I"):
        if token in value:
            return roman_map[token]
    return normalize_stage(value)


def build_gse15459_acrg_labels(outcomes_path: Path, source_data_path: Path) -> pd.DataFrame:
    outcomes = pd.read_csv(outcomes_path, sep="\t", dtype=str)
    source_rows = parse_source_sheet(source_data_path, "Singapore")
    outcomes["sample_id"] = outcomes["cel_file"].str.replace(".CEL", "", regex=False)
    outcomes = outcomes[outcomes["sample_id"].isin(source_rows)].copy()
    outcomes["subtype"] = outcomes["sample_id"].map(lambda sample_id: source_rows[sample_id]["molecular_subtype_label"])
    outcomes["overall_survival_months"] = pd.to_numeric(outcomes["overall_survival_months"])
    outcomes["overall_survival_event"] = pd.to_numeric(outcomes["overall_survival_event"])
    outcomes["stage_numeric"] = outcomes["stage"].map(normalize_stage_roman)
    return outcomes


def build_gse84437_state_frame(clinical_path: Path, assignment_path: Path) -> pd.DataFrame:
    clinical = pd.read_csv(clinical_path, sep="\t", dtype=str)
    assignments = pd.read_csv(assignment_path, sep="\t", dtype=str)
    merged = clinical.merge(assignments[["gsm_id", "predicted_subtype"]], on="gsm_id", how="inner")
    merged["subtype"] = merged["predicted_subtype"]
    merged["overall_survival_months"] = pd.to_numeric(merged["overall_survival_months"])
    merged["overall_survival_event"] = pd.to_numeric(merged["overall_survival_event"])
    merged["pt_stage_numeric"] = merged["pt_stage"].map(normalize_stage)
    merged["pn_stage_numeric"] = merged["pn_stage"].map(normalize_stage)
    return merged


def build_gse62254_state_frame(clinical_path: Path) -> pd.DataFrame:
    clinical = pd.read_csv(clinical_path, sep="\t", dtype=str)
    clinical["overall_survival_months"] = pd.to_numeric(clinical["overall_survival_months"])
    clinical["overall_survival_event"] = pd.to_numeric(clinical["overall_survival_months"].notna() & clinical["overall_survival_months"].astype(float).ge(0)).astype(int)
    clinical["overall_survival_event"] = pd.to_numeric(clinical["follow_up_status"].fillna("0").replace({"0": 0, "3": 1}), errors="coerce").fillna(0).astype(int)
    clinical["stage_numeric"] = clinical["pathologic_stage"].map(normalize_stage_roman)
    clinical["subtype"] = clinical["molecular_subtype_label"]
    return clinical


def derive_gene_panel(training_gene_matrix: pd.DataFrame, training_clinical: pd.DataFrame, shared_genes: pd.Index) -> tuple[list[dict[str, str]], list[str], list[str]]:
    subtype_by_sample = training_clinical.set_index("gsm_id")["molecular_subtype_label"]
    subtype_by_sample = subtype_by_sample.loc[
        subtype_by_sample.index.isin(training_gene_matrix.columns)
        & subtype_by_sample.isin([REFERENCE_SUBTYPE, "EMT"])
    ]
    standardized = zscore_by_gene(training_gene_matrix.loc[shared_genes, subtype_by_sample.index])
    emt_mean = standardized.loc[:, subtype_by_sample[subtype_by_sample == "EMT"].index].mean(axis=1)
    msi_mean = standardized.loc[:, subtype_by_sample[subtype_by_sample == REFERENCE_SUBTYPE].index].mean(axis=1)
    delta = (emt_mean - msi_mean).sort_values(ascending=False)

    positive = delta.head(TOP_GENES_PER_DIRECTION)
    negative = delta.tail(TOP_GENES_PER_DIRECTION).sort_values(ascending=True)
    rows: list[dict[str, str]] = []
    for rank, (gene_symbol, value) in enumerate(positive.items(), start=1):
        rows.append(
            {
                "gene_symbol": gene_symbol,
                "direction": "EMT_high",
                "weight": "1",
                "derivation_delta": format_number(value),
                "rank_within_direction": str(rank),
            }
        )
    for rank, (gene_symbol, value) in enumerate(negative.items(), start=1):
        rows.append(
            {
                "gene_symbol": gene_symbol,
                "direction": "MSI_high",
                "weight": "-1",
                "derivation_delta": format_number(value),
                "rank_within_direction": str(rank),
            }
        )
    return rows, list(positive.index), list(negative.index)


def score_gene_matrix(gene_matrix: pd.DataFrame, positive_genes: list[str], negative_genes: list[str]) -> pd.Series:
    panel_genes = positive_genes + negative_genes
    standardized = zscore_by_gene(gene_matrix.loc[panel_genes])
    score = standardized.loc[positive_genes].mean(axis=0) - standardized.loc[negative_genes].mean(axis=0)
    score = (score - score.mean()) / score.std(ddof=0)
    return score.rename("compact_score_z")


def fit_state_model(frame: pd.DataFrame, dataset_id: str, model_name: str) -> dict[str, str]:
    covariate_columns: list[str] = []
    base_columns = ["overall_survival_months", "overall_survival_event", "subtype"]
    if model_name == "cox_stage_adjusted":
        base_columns.append("stage_numeric")
    elif model_name == "cox_pt_pn_adjusted":
        base_columns.extend(["pt_stage_numeric", "pn_stage_numeric"])
    model_frame = frame[base_columns].copy()
    subtype_dummies = pd.get_dummies(model_frame["subtype"], prefix="subtype")
    reference_column = f"subtype_{REFERENCE_SUBTYPE}"
    if reference_column in subtype_dummies.columns:
        subtype_dummies = subtype_dummies.drop(columns=[reference_column])
    if "subtype_EMT" not in subtype_dummies.columns:
        raise ValueError(f"EMT comparison missing for {dataset_id} {model_name}")
    covariate_columns.extend(sorted(subtype_dummies.columns))
    model_frame = pd.concat([model_frame, subtype_dummies], axis=1)

    covariates_label = "none"
    if model_name == "cox_stage_adjusted":
        model_frame = model_frame.dropna(subset=["stage_numeric"]).copy()
        model_frame["stage_numeric"] = pd.to_numeric(model_frame["stage_numeric"])
        covariate_columns.append("stage_numeric")
        covariates_label = "stage"
    elif model_name == "cox_pt_pn_adjusted":
        model_frame = model_frame.dropna(subset=["pt_stage_numeric", "pn_stage_numeric"]).copy()
        model_frame["pt_stage_numeric"] = pd.to_numeric(model_frame["pt_stage_numeric"])
        model_frame["pn_stage_numeric"] = pd.to_numeric(model_frame["pn_stage_numeric"])
        covariate_columns.extend(["pt_stage_numeric", "pn_stage_numeric"])
        covariates_label = "pt_stage,pn_stage"

    cox_frame = model_frame[["overall_survival_months", "overall_survival_event", *covariate_columns]].copy()
    fitter = CoxPHFitter()
    fitter.fit(cox_frame, duration_col="overall_survival_months", event_col="overall_survival_event")
    summary = fitter.summary.loc["subtype_EMT"]
    return {
        "method_name": BASELINE_METHOD,
        "method_class": "baseline",
        "dataset_id": dataset_id,
        "model": model_name,
        "feature_definition": "4-state frozen shared ACRG label",
        "contrast": "EMT_vs_MSI",
        "primary_metric": "hazard_ratio",
        "primary_metric_value": format_number(summary["exp(coef)"]),
        "ci_lower": format_number(summary["exp(coef) lower 95%"]),
        "ci_upper": format_number(summary["exp(coef) upper 95%"]),
        "pvalue": format_number(summary["p"]),
        "secondary_metric": "c_index",
        "secondary_metric_value": format_number(fitter.concordance_index_),
        "n": str(len(model_frame)),
        "events": str(int(model_frame["overall_survival_event"].sum())),
        "covariates": covariates_label,
    }


def fit_compact_score_model(frame: pd.DataFrame, dataset_id: str, model_name: str) -> dict[str, str]:
    covariate_columns = ["compact_score_z"]
    model_frame = frame[["overall_survival_months", "overall_survival_event", "compact_score_z"]].copy()

    covariates_label = "compact_score"
    if model_name == "cox_score_stage_adjusted":
        model_frame = frame[["overall_survival_months", "overall_survival_event", "compact_score_z", "stage_numeric"]].dropna().copy()
        covariate_columns.append("stage_numeric")
        covariates_label = "compact_score,stage"
    elif model_name == "cox_score_pt_pn_adjusted":
        model_frame = frame[
            ["overall_survival_months", "overall_survival_event", "compact_score_z", "pt_stage_numeric", "pn_stage_numeric"]
        ].dropna().copy()
        covariate_columns.extend(["pt_stage_numeric", "pn_stage_numeric"])
        covariates_label = "compact_score,pt_stage,pn_stage"

    fitter = CoxPHFitter()
    fitter.fit(model_frame[["overall_survival_months", "overall_survival_event", *covariate_columns]], duration_col="overall_survival_months", event_col="overall_survival_event")
    summary = fitter.summary.loc["compact_score_z"]
    return {
        "method_name": COMPACT_METHOD,
        "method_class": "candidate",
        "dataset_id": dataset_id,
        "model": model_name,
        "feature_definition": f"{TOP_GENES_PER_DIRECTION}+{TOP_GENES_PER_DIRECTION} derivation-only EMT-vs-MSI gene program score",
        "contrast": "per_sd_increase",
        "primary_metric": "hazard_ratio",
        "primary_metric_value": format_number(summary["exp(coef)"]),
        "ci_lower": format_number(summary["exp(coef) lower 95%"]),
        "ci_upper": format_number(summary["exp(coef) upper 95%"]),
        "pvalue": format_number(summary["p"]),
        "secondary_metric": "c_index",
        "secondary_metric_value": format_number(fitter.concordance_index_),
        "n": str(len(model_frame)),
        "events": str(int(model_frame["overall_survival_event"].sum())),
        "covariates": covariates_label,
    }


def attach_compact_delta(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    baseline_map = {
        (row["dataset_id"], row["model"]): float(row["secondary_metric_value"])
        for row in rows
        if row["method_name"] == BASELINE_METHOD
    }
    enriched: list[dict[str, str]] = []
    for row in rows:
        updated = dict(row)
        baseline_model = row["model"].replace("cox_score_", "cox_") if row["method_name"] == COMPACT_METHOD else row["model"]
        baseline_cindex = baseline_map.get((row["dataset_id"], baseline_model))
        if baseline_cindex is None:
            updated["c_index_delta_vs_baseline"] = "0"
        else:
            updated["c_index_delta_vs_baseline"] = format_number(float(row["secondary_metric_value"]) - baseline_cindex)
        enriched.append(updated)
    return enriched


def build_decision_rows(rows: list[dict[str, str]], panel_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    external_compact = [
        row
        for row in rows
        if row["method_name"] == COMPACT_METHOD and row["dataset_id"] in {"GSE15459", "GSE84437"}
    ]
    direction_pass = all(float(row["primary_metric_value"]) > 1.0 for row in external_compact)
    deltas = [float(row["c_index_delta_vs_baseline"]) for row in external_compact]
    noninferior_pass = all(delta >= -0.02 for delta in deltas)
    gain_pass = any(delta >= 0.01 for delta in deltas) and (sum(deltas) / len(deltas) >= 0)
    replace_pass = direction_pass and noninferior_pass and gain_pass
    return [
        {
            "candidate_method": COMPACT_METHOD,
            "gene_panel_size": str(len(panel_rows)),
            "replace_threshold": "all external HR>1 in matched models; no external c_index drop worse than -0.02 vs shared state; at least one external c_index gain >=0.01 and mean external delta >=0",
            "external_direction_pass": str(direction_pass).lower(),
            "external_cindex_noninferior_pass": str(noninferior_pass).lower(),
            "external_gain_pass": str(gain_pass).lower(),
            "replace_decision": "replace_with_compact_score" if replace_pass else "keep_shared_acrg_state",
            "decision_rationale": "compact score stays secondary until it beats the shared-state baseline under predeclared external criteria",
        }
    ]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Benchmark shared ACRG state against a derivation-only compact score")
    parser.add_argument("--training-clinical-input", required=True)
    parser.add_argument("--training-matrix-input", required=True)
    parser.add_argument("--gse15459-outcomes-input", required=True)
    parser.add_argument("--gse15459-matrix-input", required=True)
    parser.add_argument("--gse84437-clinical-input", required=True)
    parser.add_argument("--gse84437-matrix-input", required=True)
    parser.add_argument("--gse84437-assignment-input", required=True)
    parser.add_argument("--source-data-input", required=True)
    parser.add_argument("--gpl570-map-input", required=True)
    parser.add_argument("--gpl6947-map-input", required=True)
    parser.add_argument("--benchmark-output", required=True)
    parser.add_argument("--decision-output", required=True)
    parser.add_argument("--gene-panel-output", required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    training_clinical_path = Path(args.training_clinical_input)
    training_matrix_path = Path(args.training_matrix_input)
    gse15459_outcomes_path = Path(args.gse15459_outcomes_input)
    gse15459_matrix_path = Path(args.gse15459_matrix_input)
    gse84437_clinical_path = Path(args.gse84437_clinical_input)
    gse84437_matrix_path = Path(args.gse84437_matrix_input)
    gse84437_assignment_path = Path(args.gse84437_assignment_input)
    source_data_path = Path(args.source_data_input)
    gpl570_map_path = Path(args.gpl570_map_input)
    gpl6947_map_path = Path(args.gpl6947_map_input)

    training_clinical = pd.read_csv(training_clinical_path, sep="\t", dtype=str)
    gse15459_state = build_gse15459_acrg_labels(gse15459_outcomes_path, source_data_path)
    gse84437_state = build_gse84437_state_frame(gse84437_clinical_path, gse84437_assignment_path)
    gse62254_state = build_gse62254_state_frame(training_clinical_path)

    training_matrix = load_expression_matrix(training_matrix_path)
    gse15459_matrix = load_expression_matrix(gse15459_matrix_path)
    gse84437_matrix = load_expression_matrix(gse84437_matrix_path)

    gpl570_probe_ids = set(training_matrix.index.astype(str)).union(gse15459_matrix.index.astype(str))
    gpl6947_probe_ids = set(gse84437_matrix.index.astype(str))
    gpl570_map = load_or_build_probe_gene_map(gpl570_map_path, "GPL570", "Gene Symbol", gpl570_probe_ids)
    gpl6947_map = load_or_build_probe_gene_map(gpl6947_map_path, "GPL6947", "Symbol", gpl6947_probe_ids)

    training_gene_matrix = collapse_to_gene_matrix(training_matrix, gpl570_map)
    gse15459_gene_matrix = collapse_to_gene_matrix(gse15459_matrix, gpl570_map)
    gse84437_gene_matrix = collapse_to_gene_matrix(gse84437_matrix, gpl6947_map)

    shared_genes = training_gene_matrix.index.intersection(gse15459_gene_matrix.index).intersection(gse84437_gene_matrix.index)
    panel_rows, positive_genes, negative_genes = derive_gene_panel(training_gene_matrix, training_clinical, shared_genes)

    gse62254_score = score_gene_matrix(training_gene_matrix, positive_genes, negative_genes)
    gse15459_score = score_gene_matrix(gse15459_gene_matrix, positive_genes, negative_genes)
    gse84437_score = score_gene_matrix(gse84437_gene_matrix, positive_genes, negative_genes)

    gse62254_state = gse62254_state.merge(gse62254_score.rename_axis("gsm_id").reset_index(), on="gsm_id", how="inner")
    gse15459_state = gse15459_state.merge(gse15459_score.rename_axis("gsm_id").reset_index(), on="gsm_id", how="inner")
    gse84437_state = gse84437_state.merge(gse84437_score.rename_axis("gsm_id").reset_index(), on="gsm_id", how="inner")

    rows = [
        fit_state_model(gse62254_state, "GSE62254", "cox_unadjusted"),
        fit_state_model(gse62254_state, "GSE62254", "cox_stage_adjusted"),
        fit_state_model(gse15459_state, "GSE15459", "cox_unadjusted"),
        fit_state_model(gse15459_state, "GSE15459", "cox_stage_adjusted"),
        fit_state_model(gse84437_state, "GSE84437", "cox_unadjusted"),
        fit_state_model(gse84437_state, "GSE84437", "cox_pt_pn_adjusted"),
        fit_compact_score_model(gse62254_state, "GSE62254", "cox_score_unadjusted"),
        fit_compact_score_model(gse62254_state, "GSE62254", "cox_score_stage_adjusted"),
        fit_compact_score_model(gse15459_state, "GSE15459", "cox_score_unadjusted"),
        fit_compact_score_model(gse15459_state, "GSE15459", "cox_score_stage_adjusted"),
        fit_compact_score_model(gse84437_state, "GSE84437", "cox_score_unadjusted"),
        fit_compact_score_model(gse84437_state, "GSE84437", "cox_score_pt_pn_adjusted"),
    ]
    benchmark_rows = attach_compact_delta(rows)
    decision_rows = build_decision_rows(benchmark_rows, panel_rows)

    write_tsv(Path(args.benchmark_output), benchmark_rows)
    write_tsv(Path(args.decision_output), decision_rows)
    write_tsv(Path(args.gene_panel_output), panel_rows)


if __name__ == "__main__":
    main()
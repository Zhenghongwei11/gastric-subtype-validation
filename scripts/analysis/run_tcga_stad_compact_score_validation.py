from __future__ import annotations

import argparse
import csv
from pathlib import Path

import pandas as pd
import xenaPython as xena
from lifelines import CoxPHFitter


XENA_HUB = "https://gdc.xenahubs.net"
EXPRESSION_DATASET = "TCGA-STAD.star_fpkm-uq.tsv"
SURVIVAL_DATASET = "TCGA-STAD.survival.tsv"
CLINICAL_DATASET = "TCGA-STAD.clinical.tsv"
AGE_FIELD = "age_at_earliest_diagnosis_in_years.diagnoses.xena_derived"
SEX_FIELD = "gender.demographic"
STAGE_FIELD = "ajcc_pathologic_stage.diagnoses"


def format_number(value: float) -> str:
    return format(float(value), ".6g")


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def parse_stage_label(raw_value: str) -> int | None:
    value = str(raw_value).strip().upper()
    if not value:
        return None
    roman_map = {"I": 1, "II": 2, "III": 3, "IV": 4}
    for token in ("IV", "III", "II", "I"):
        if token in value:
            return roman_map[token]
    return None


def decode_codes(code_string: str) -> list[str]:
    return code_string.split("\t") if code_string else []


def decode_value(raw_value: object, labels: list[str]) -> str:
    try:
        index = int(raw_value)
    except (TypeError, ValueError):
        return str(raw_value)
    return labels[index] if 0 <= index < len(labels) else ""


def primary_tumor_samples(samples: list[str]) -> list[str]:
    selected = []
    for sample in samples:
        tokens = sample.split("-")
        if len(tokens) < 4:
            continue
        sample_code = tokens[3][:2]
        if sample_code == "01":
            selected.append(sample)
    return selected


def load_gene_panel(path: Path) -> tuple[list[str], list[str]]:
    frame = pd.read_csv(path, sep="\t", dtype=str)
    positive = frame.loc[frame["direction"] == "EMT_high", "gene_symbol"].tolist()
    negative = frame.loc[frame["direction"] == "MSI_high", "gene_symbol"].tolist()
    return positive, negative


def query_tcga_frame(positive_genes: list[str], negative_genes: list[str]) -> pd.DataFrame:
    tumor_samples = primary_tumor_samples(xena.dataset_samples(XENA_HUB, EXPRESSION_DATASET, None))
    genes = positive_genes + negative_genes
    expression = xena.dataset_gene_probe_avg(XENA_HUB, EXPRESSION_DATASET, tumor_samples, genes)
    expression_frame = pd.DataFrame(
        {entry["gene"]: entry["scores"][0] for entry in expression},
        index=tumor_samples,
    )
    expression_frame = expression_frame.apply(pd.to_numeric, errors="coerce")

    _, survival_values = xena.dataset_probe_values(XENA_HUB, SURVIVAL_DATASET, tumor_samples, ["OS", "OS.time"])
    survival_frame = pd.DataFrame(
        {
            "sample_id": tumor_samples,
            "overall_survival_event": pd.to_numeric(survival_values[0], errors="coerce"),
            "overall_survival_months": pd.to_numeric(survival_values[1], errors="coerce"),
        }
    ).set_index("sample_id")

    clinical_fields = [SEX_FIELD, AGE_FIELD, STAGE_FIELD]
    _, clinical_values = xena.dataset_probe_values(XENA_HUB, CLINICAL_DATASET, tumor_samples, clinical_fields)
    codebook = xena.field_codes(XENA_HUB, CLINICAL_DATASET, [SEX_FIELD, STAGE_FIELD])
    sex_labels = decode_codes(codebook[0]["code"])
    stage_labels = decode_codes(codebook[1]["code"])
    clinical_frame = pd.DataFrame(
        {
            "sample_id": tumor_samples,
            "sex": [decode_value(value, sex_labels) for value in clinical_values[0]],
            "age_years": pd.to_numeric(clinical_values[1], errors="coerce"),
            "stage_label": [decode_value(value, stage_labels) for value in clinical_values[2]],
        }
    ).set_index("sample_id")
    clinical_frame["stage_numeric"] = clinical_frame["stage_label"].map(parse_stage_label)
    clinical_frame["female_indicator"] = clinical_frame["sex"].map({"female": 1.0, "male": 0.0})

    standardized = expression_frame.sub(expression_frame.mean(axis=0), axis=1).div(expression_frame.std(axis=0, ddof=0), axis=1)
    standardized = standardized.dropna(axis=1, how="all")
    score = standardized[positive_genes].mean(axis=1) - standardized[negative_genes].mean(axis=1)
    score = (score - score.mean()) / score.std(ddof=0)

    merged = pd.concat([survival_frame, clinical_frame], axis=1)
    merged["compact_score_z"] = score
    merged = merged.dropna(subset=["overall_survival_event", "overall_survival_months", "compact_score_z"]).copy()
    return merged.reset_index().rename(columns={"index": "sample_id"})


def fit_model(frame: pd.DataFrame, model_name: str) -> dict[str, str]:
    covariate_columns = ["compact_score_z"]
    model_frame = frame[["overall_survival_months", "overall_survival_event", "compact_score_z"]].copy()
    covariates_label = "compact_score"

    if model_name == "cox_score_stage_age_sex_adjusted":
        model_frame = frame[
            ["overall_survival_months", "overall_survival_event", "compact_score_z", "stage_numeric", "age_years", "female_indicator"]
        ].dropna().copy()
        covariate_columns.extend(["stage_numeric", "age_years", "female_indicator"])
        covariates_label = "compact_score,stage,age,sex"

    fitter = CoxPHFitter()
    fitter.fit(
        model_frame[["overall_survival_months", "overall_survival_event", *covariate_columns]],
        duration_col="overall_survival_months",
        event_col="overall_survival_event",
    )
    summary = fitter.summary.loc["compact_score_z"]
    return {
        "dataset_id": "TCGA-STAD",
        "outcome": "overall_survival",
        "model": model_name,
        "reference_subtype": "",
        "comparison_subtype": "compact_score",
        "effect_type": "hazard_ratio",
        "effect": format_number(summary["exp(coef)"]),
        "ci_lower": format_number(summary["exp(coef) lower 95%"]),
        "ci_upper": format_number(summary["exp(coef) upper 95%"]),
        "pvalue": format_number(summary["p"]),
        "n": str(len(model_frame)),
        "events": str(int(model_frame["overall_survival_event"].sum())),
        "covariates": covariates_label,
        "secondary_metric": "c_index",
        "secondary_metric_value": format_number(fitter.concordance_index_),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run TCGA-STAD compact-score orthogonal validation via UCSC Xena")
    parser.add_argument("--gene-panel-input", required=True)
    parser.add_argument("--output", required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    positive_genes, negative_genes = load_gene_panel(Path(args.gene_panel_input))
    frame = query_tcga_frame(positive_genes, negative_genes)
    rows = [
        fit_model(frame, "cox_score_unadjusted"),
        fit_model(frame, "cox_score_stage_age_sex_adjusted"),
    ]
    write_tsv(Path(args.output), rows)


if __name__ == "__main__":
    main()
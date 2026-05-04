from __future__ import annotations

import argparse
import csv
import gzip
import io
import urllib.request
from pathlib import Path

import pandas as pd
from lifelines import CoxPHFitter

from run_gse84437_projection_survival import (
    collapse_to_gene_matrix,
    load_or_build_probe_gene_map,
    load_probe_gene_map,
    parse_stage_token,
    zscore_by_gene,
)


GPL8432 = "GPL8432"
GPL8432_SYMBOL_CANDIDATES = ["Symbol", "ILMN_Gene", "Gene Symbol", "SYMBOL"]


def format_number(value: float) -> str:
    return format(float(value), ".6g")


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def load_gene_panel(path: Path) -> tuple[list[str], list[str]]:
    frame = pd.read_csv(path, sep="\t", dtype=str)
    positive = frame.loc[frame["direction"] == "EMT_high", "gene_symbol"].tolist()
    negative = frame.loc[frame["direction"] == "MSI_high", "gene_symbol"].tolist()
    return positive, negative


def parse_characteristic(raw_value: str) -> tuple[str, str]:
    value = str(raw_value).strip().strip('"')
    if ":" not in value:
        return value.lower(), ""
    label, parsed_value = value.split(":", 1)
    return label.strip().lower(), parsed_value.strip().strip('"')


def load_gse26253_expression_and_clinical(path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    sample_ids: list[str] = []
    sample_records: dict[str, dict[str, object]] = {}
    table_lines: list[str] = []
    in_table = False

    with gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t", quotechar='"')
        for row in reader:
            if not row:
                continue
            tag = str(row[0]).strip()
            if tag == "!series_matrix_table_begin":
                in_table = True
                continue
            if tag == "!series_matrix_table_end":
                break

            if in_table:
                table_lines.append("\t".join(row))
                continue

            if tag == "!Sample_geo_accession":
                sample_ids = [str(value).strip() for value in row[1:] if str(value).strip()]
                sample_records = {sample_id: {"gsm_id": sample_id} for sample_id in sample_ids}
                continue

            if not tag.startswith("!Sample_characteristics_ch") or not sample_ids:
                continue

            values = [str(value) for value in row[1 : 1 + len(sample_ids)]]
            for sample_id, raw_value in zip(sample_ids, values, strict=True):
                label, parsed_value = parse_characteristic(raw_value)
                if label.startswith("status "):
                    sample_records[sample_id]["recurrence_event"] = pd.to_numeric(parsed_value, errors="coerce")
                elif label.startswith("recurrence free survival time"):
                    sample_records[sample_id]["recurrence_free_survival_months"] = pd.to_numeric(
                        parsed_value,
                        errors="coerce",
                    )
                elif label.startswith("pathological stage"):
                    sample_records[sample_id]["stage_label"] = parsed_value

    matrix = pd.read_csv(io.StringIO("\n".join(table_lines)), sep="\t", index_col=0)
    matrix.index = matrix.index.astype(str)
    matrix.columns = matrix.columns.astype(str)
    matrix = matrix.apply(pd.to_numeric, errors="coerce")

    clinical = pd.DataFrame(sample_records.values())
    clinical["stage_label"] = clinical.get("stage_label", "")
    clinical["stage_numeric"] = clinical["stage_label"].map(parse_stage_token)
    clinical["recurrence_event"] = pd.to_numeric(clinical["recurrence_event"], errors="coerce")
    clinical["recurrence_free_survival_months"] = pd.to_numeric(
        clinical["recurrence_free_survival_months"],
        errors="coerce",
    )
    return matrix, clinical


def resolve_platform_symbol_column(platform_accession: str) -> str:
    url = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={platform_accession}&targ=self&form=text&view=data"
    with urllib.request.urlopen(url, timeout=180) as response:
        for raw in response:
            line = raw.decode("utf-8", "replace").strip()
            if line != "!platform_table_begin":
                continue
            header = next(response).decode("utf-8", "replace").strip().split("\t")
            for candidate in GPL8432_SYMBOL_CANDIDATES:
                if candidate in header:
                    return candidate
            for field in header:
                lowered = field.lower()
                if "symbol" in lowered or lowered == "ilmn_gene":
                    return field
            raise ValueError(f"Could not identify a gene-symbol column for {platform_accession}")
    raise ValueError(f"Could not locate platform table for {platform_accession}")


def load_gpl8432_probe_gene_map(path: Path, probe_ids: set[str]) -> dict[str, str]:
    if path.exists():
        mapping = load_probe_gene_map(path)
        if mapping:
            return mapping

    symbol_column = resolve_platform_symbol_column(GPL8432)
    return load_or_build_probe_gene_map(path, GPL8432, symbol_column, probe_ids)


def build_validation_frame(
    matrix_path: Path,
    gene_panel_path: Path,
    gpl8432_map_path: Path,
) -> tuple[pd.DataFrame, int, int]:
    positive_genes, negative_genes = load_gene_panel(gene_panel_path)
    probe_matrix, clinical = load_gse26253_expression_and_clinical(matrix_path)
    probe_gene_map = load_gpl8432_probe_gene_map(gpl8432_map_path, set(probe_matrix.index.astype(str)))
    gene_matrix = collapse_to_gene_matrix(probe_matrix, probe_gene_map)

    available_positive = [gene for gene in positive_genes if gene in gene_matrix.index]
    available_negative = [gene for gene in negative_genes if gene in gene_matrix.index]
    panel_genes = available_positive + available_negative
    if not available_positive or not available_negative:
        raise ValueError("Insufficient GSE26253 overlap with frozen compact-score gene panel")

    standardized = zscore_by_gene(gene_matrix.loc[panel_genes])
    score = standardized.loc[available_positive].mean(axis=0) - standardized.loc[available_negative].mean(axis=0)
    score = (score - score.mean()) / score.std(ddof=0)

    frame = clinical.merge(score.rename("compact_score_z").rename_axis("gsm_id").reset_index(), on="gsm_id", how="inner")
    frame = frame.dropna(subset=["recurrence_free_survival_months", "recurrence_event", "compact_score_z"]).copy()
    return frame, len(available_positive), len(available_negative)


def fit_model(frame: pd.DataFrame, model_name: str, positive_count: int, negative_count: int) -> dict[str, str]:
    covariate_columns = ["compact_score_z"]
    model_frame = frame[["recurrence_free_survival_months", "recurrence_event", "compact_score_z"]].copy()
    covariates_label = "compact_score"

    if model_name == "cox_score_stage_adjusted":
        model_frame = frame[
            ["recurrence_free_survival_months", "recurrence_event", "compact_score_z", "stage_numeric"]
        ].dropna().copy()
        covariate_columns.append("stage_numeric")
        covariates_label = "compact_score,stage"

    fitter = CoxPHFitter()
    fitter.fit(
        model_frame[["recurrence_free_survival_months", "recurrence_event", *covariate_columns]],
        duration_col="recurrence_free_survival_months",
        event_col="recurrence_event",
    )
    summary = fitter.summary.loc["compact_score_z"]
    return {
        "dataset_id": "GSE26253",
        "outcome": "recurrence_free_survival",
        "model": model_name,
        "reference_subtype": "",
        "comparison_subtype": "compact_score",
        "effect_type": "hazard_ratio",
        "effect": format_number(summary["exp(coef)"]),
        "ci_lower": format_number(summary["exp(coef) lower 95%"]),
        "ci_upper": format_number(summary["exp(coef) upper 95%"]),
        "pvalue": format_number(summary["p"]),
        "n": str(len(model_frame)),
        "events": str(int(model_frame["recurrence_event"].sum())),
        "covariates": covariates_label,
        "positive_genes_available": str(positive_count),
        "negative_genes_available": str(negative_count),
        "secondary_metric": "c_index",
        "secondary_metric_value": format_number(fitter.concordance_index_),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run compact-score validation in the recurrence-oriented GSE26253 cohort")
    parser.add_argument("--gene-panel-input", required=True)
    parser.add_argument("--matrix-input", required=True)
    parser.add_argument("--gpl8432-map-input", required=True)
    parser.add_argument("--output", required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    frame, positive_count, negative_count = build_validation_frame(
        Path(args.matrix_input),
        Path(args.gene_panel_input),
        Path(args.gpl8432_map_input),
    )
    rows = [
        fit_model(frame, "cox_score_unadjusted", positive_count, negative_count),
        fit_model(frame, "cox_score_stage_adjusted", positive_count, negative_count),
    ]
    write_tsv(Path(args.output), rows)


if __name__ == "__main__":
    main()
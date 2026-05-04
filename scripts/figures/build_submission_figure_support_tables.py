from __future__ import annotations

import argparse
import csv
import gzip
import math
from pathlib import Path

import pandas as pd
import xlrd


SUBTYPE_ORDER = ["EMT", "MSI", "MSS/TP53-", "MSS/TP53+"]
MARKER_DEFINITIONS = [
    ("VIM", "EMT_high"),
    ("ZEB1", "EMT_high"),
    ("CDH1", "epithelial_high"),
    ("EPCAM", "epithelial_high"),
    ("CLDN18", "epithelial_high"),
    ("MKI67", "proliferation_context"),
]
COHORT_ORDER = {
    "GSE62254": 1,
    "GSE15459": 2,
    "GSE84437": 3,
    "TCGA-STAD": 4,
    "GSE26899": 5,
    "GSE26901": 6,
}
SCORE_COLUMN_TO_SUBTYPE = {
    "score_emt": "EMT",
    "score_msi": "MSI",
    "score_mss_tp53_minus": "MSS/TP53-",
    "score_mss_tp53_plus": "MSS/TP53+",
}
PATHWAY_LABELS = {
    "EMT_high": "EMT pathway",
    "epithelial_high": "Epithelial markers",
    "proliferation_context": "Proliferation context",
}


def format_number(value: float) -> str:
    return format(float(value), ".6g")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return [{key: (value or "") for key, value in row.items()} for row in reader]


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    if not rows:
        raise ValueError(f"Cannot write empty TSV for {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def parse_source_sheet(path: Path, sheet_name: str) -> dict[str, dict[str, str]]:
    workbook = xlrd.open_workbook(str(path))
    sheet = workbook.sheet_by_name(sheet_name)
    headers = [str(value).strip() for value in sheet.row_values(0)]
    sample_idx = headers.index("Sample ID")
    subtype_idx = headers.index("Mol subtype")

    subtype_labels = {
        "0": "MSS/TP53-",
        "1": "MSS/TP53+",
        "2": "MSI",
        "3": "EMT",
    }

    rows: dict[str, dict[str, str]] = {}
    for row_index in range(1, sheet.nrows):
        sample_id = str(sheet.cell_value(row_index, sample_idx)).strip()
        if sample_id.endswith(".0"):
            sample_id = sample_id[:-2]
        subtype_code = str(int(sheet.cell_value(row_index, subtype_idx)))
        rows[sample_id] = {
            "molecular_subtype_code": subtype_code,
            "molecular_subtype_label": subtype_labels[subtype_code],
        }
    return rows


def normalize_stage_group(raw_value: str) -> str:
    value = str(raw_value).strip().upper()
    if "III" in value or value == "3":
        return "Stage III"
    return ""


def load_probe_gene_map(path: Path) -> dict[str, str]:
    frame = pd.read_csv(path, sep="\t", dtype=str)
    frame = frame.dropna(subset=["probe_id", "gene_symbol"]).copy()
    frame["probe_id"] = frame["probe_id"].str.strip()
    frame["gene_symbol"] = frame["gene_symbol"].str.strip()
    frame = frame[(frame["probe_id"] != "") & (frame["gene_symbol"] != "")].copy()
    return dict(zip(frame["probe_id"], frame["gene_symbol"]))


def load_expression_matrix(path: Path) -> pd.DataFrame:
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        frame = pd.read_csv(handle, sep="\t", index_col=0)
    frame.index = frame.index.astype(str)
    frame.columns = frame.columns.astype(str)
    return frame.apply(pd.to_numeric)


def collapse_to_gene_matrix(matrix: pd.DataFrame, probe_gene_map: dict[str, str]) -> pd.DataFrame:
    probe_index = [probe for probe in matrix.index if probe in probe_gene_map]
    if not probe_index:
        raise ValueError("No overlapping probes found between matrix and GPL570 gene-symbol map")
    collapsed = matrix.loc[probe_index].copy()
    collapsed.insert(0, "gene_symbol", [probe_gene_map[probe] for probe in probe_index])
    collapsed = collapsed.groupby("gene_symbol", sort=True).mean(numeric_only=True)
    collapsed.index = collapsed.index.astype(str)
    return collapsed


def build_figure1_cohort_map_rows(dataset_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for dataset_id in [
        "GSE62254",
        "GSE15459",
        "GSE84437",
        "TCGA-STAD",
        "GSE26899",
        "GSE26901",
    ]:
        dataset_row = next((row for row in dataset_rows if row["dataset_id"] == dataset_id), None)
        if dataset_row is None:
            continue
        rows.append(
            {
                "cohort_order": str(COHORT_ORDER[dataset_id]),
                "dataset_id": dataset_row["dataset_id"],
                "sample_size": dataset_row["sample_size"],
                "planned_role": dataset_row["planned_role"],
                "platform": dataset_row["platform"],
                "endpoint_availability": dataset_row["endpoint_availability"],
                "treatment_metadata": dataset_row["treatment_metadata"],
                "resource_type": dataset_row["resource_type"],
            }
        )
    return rows


def build_figure1_projection_rows(projection_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    return [dict(row) for row in projection_rows]


def build_figure1_projection_space_rows(projected_assignment_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for row in projected_assignment_rows:
        score_pairs = sorted(
            [
                (subtype, float(row[column_name]))
                for column_name, subtype in SCORE_COLUMN_TO_SUBTYPE.items()
            ],
            key=lambda item: item[1],
            reverse=True,
        )
        rows.append(
            {
                "dataset_id": row["dataset_id"],
                "gsm_id": row["gsm_id"],
                "predicted_subtype": row["predicted_subtype"],
                "confidence_tier": row["assignment_margin_rank"],
                "assignment_margin": row["assignment_margin"],
                "x_score": format_number(float(row["score_emt"]) - float(row["score_msi"])),
                "y_score": format_number(float(row["score_mss_tp53_minus"]) - float(row["score_mss_tp53_plus"])),
                "primary_axis_subtype": score_pairs[0][0],
                "secondary_axis_subtype": score_pairs[1][0],
                "primary_axis_score": format_number(score_pairs[0][1]),
                "secondary_axis_score": format_number(score_pairs[1][1]),
            }
        )
    return rows


def overall_survival_event_from_status(follow_up_status: str) -> str:
    return "1" if follow_up_status in {"2", "3", "4"} else "0"


def build_figure2_derivation_km_rows(derivation_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for row in derivation_rows:
        subtype = row.get("molecular_subtype_label", "")
        survival_months = row.get("overall_survival_months", "")
        if not subtype or not survival_months:
            continue
        rows.append(
            {
                "dataset_id": row["dataset_id"],
                "sample_id": row.get("patient_id") or row.get("gsm_id", ""),
                "gsm_id": row.get("gsm_id", ""),
                "molecular_subtype_label": subtype,
                "overall_survival_months": survival_months,
                "overall_survival_event": overall_survival_event_from_status(row.get("follow_up_status", "")),
                "pathologic_stage": row.get("pathologic_stage", ""),
                "follow_up_status": row.get("follow_up_status", ""),
            }
        )
    return rows


def build_figure2_marker_expression_rows(
    derivation_clinical_input: Path,
    derivation_matrix_input: Path,
    project_root: Path,
) -> list[dict[str, str]]:
    clinical_frame = pd.read_csv(derivation_clinical_input, sep="\t", dtype=str)
    subtype_by_sample = clinical_frame.set_index("gsm_id")["molecular_subtype_label"]
    subtype_by_sample = subtype_by_sample[subtype_by_sample.isin(SUBTYPE_ORDER)]

    probe_gene_map = load_probe_gene_map(
        project_root / "data" / "cache" / "platform_annotations" / "GPL570_gene_symbols.tsv"
    )
    gene_matrix = collapse_to_gene_matrix(load_expression_matrix(derivation_matrix_input), probe_gene_map)

    aligned_sample_ids = [sample_id for sample_id in subtype_by_sample.index if sample_id in gene_matrix.columns]
    if not aligned_sample_ids:
        raise ValueError("No overlapping GSE62254 samples found between linked matrix and clinical linkage table")
    subtype_by_sample = subtype_by_sample.loc[aligned_sample_ids]
    gene_matrix = gene_matrix.loc[:, aligned_sample_ids]

    rows: list[dict[str, str]] = []
    for marker_order, (gene_symbol, marker_group) in enumerate(MARKER_DEFINITIONS, start=1):
        if gene_symbol not in gene_matrix.index:
            continue

        subtype_summaries: list[tuple[str, list[str], float, float]] = []
        for subtype in SUBTYPE_ORDER:
            sample_ids = subtype_by_sample[subtype_by_sample == subtype].index.tolist()
            if not sample_ids:
                continue
            values = pd.to_numeric(gene_matrix.loc[gene_symbol, sample_ids])
            subtype_summaries.append((subtype, sample_ids, float(values.mean()), float(values.median())))

        if not subtype_summaries:
            continue

        mean_values = [mean_value for _, _, mean_value, _ in subtype_summaries]
        center = sum(mean_values) / len(mean_values)
        variance = sum((mean_value - center) ** 2 for mean_value in mean_values) / len(mean_values)
        scale = math.sqrt(variance) if variance > 0 else 0.0
        zscores = {
            subtype: ((mean_value - center) / scale if scale else 0.0)
            for subtype, _, mean_value, _ in subtype_summaries
        }
        emt_mean = next((mean_value for subtype, _, mean_value, _ in subtype_summaries if subtype == "EMT"), None)
        msi_mean = next((mean_value for subtype, _, mean_value, _ in subtype_summaries if subtype == "MSI"), None)
        emt_minus_msi = emt_mean - msi_mean if emt_mean is not None and msi_mean is not None else None

        for subtype, sample_ids, mean_value, median_value in subtype_summaries:
            rows.append(
                {
                    "dataset_id": "GSE62254",
                    "gene_symbol": gene_symbol,
                    "marker_group": marker_group,
                    "pathway_label": PATHWAY_LABELS.get(marker_group, marker_group.replace("_", " ").title()),
                    "marker_order": str(marker_order),
                    "subtype": subtype,
                    "sample_count": str(len(sample_ids)),
                    "mean_expression": format_number(mean_value),
                    "median_expression": format_number(median_value),
                    "subtype_mean_zscore": format_number(zscores[subtype]),
                    "emt_minus_msi_mean_delta": format_number(emt_minus_msi) if emt_minus_msi is not None else "",
                }
            )

    if not rows:
        raise ValueError("No marker-expression rows could be generated from the GSE62254 linked matrix")
    return rows


def load_gse15459_acrg_rows(gse15459_outcomes_rows: list[dict[str, str]], source_data_input: Path) -> list[dict[str, str]]:
    source_rows = parse_source_sheet(source_data_input, "Singapore")
    rows: list[dict[str, str]] = []
    for row in gse15459_outcomes_rows:
        sample_id = row["cel_file"].replace(".CEL", "")
        if sample_id not in source_rows:
            continue
        rows.append(
            {
                "dataset_id": "GSE15459",
                "sample_id": sample_id,
                "subtype": source_rows[sample_id]["molecular_subtype_label"],
                "overall_survival_months": row["overall_survival_months"],
                "overall_survival_event": row["overall_survival_event"],
            }
        )
    return rows


def build_projected_survival_rows(
    projected_assignment_rows: list[dict[str, str]],
    gse84437_clinical_rows: list[dict[str, str]],
) -> list[dict[str, str]]:
    clinical_by_gsm = {row["gsm_id"]: row for row in gse84437_clinical_rows}
    rows: list[dict[str, str]] = []
    for row in projected_assignment_rows:
        clinical_row = clinical_by_gsm.get(row["gsm_id"])
        if clinical_row is None:
            continue
        rows.append(
            {
                "dataset_id": "GSE84437",
                "sample_id": row["gsm_id"],
                "subtype": row["predicted_subtype"],
                "overall_survival_months": clinical_row["overall_survival_months"],
                "overall_survival_event": clinical_row["overall_survival_event"],
            }
        )
    return rows


def build_figure3_external_km_rows(
    gse15459_rows: list[dict[str, str]],
    projected_survival_rows: list[dict[str, str]],
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for default_dataset_id, source_rows in (("GSE15459", gse15459_rows), ("GSE84437", projected_survival_rows)):
        for row in source_rows:
            subtype = row.get("subtype", row.get("predicted_subtype", ""))
            if subtype not in {"EMT", "MSI"}:
                continue
            rows.append(
                {
                    "dataset_id": row.get("dataset_id", default_dataset_id),
                    "sample_id": row.get("sample_id", row.get("gsm_id", "")),
                    "subtype": subtype,
                    "overall_survival_months": row["overall_survival_months"],
                    "overall_survival_event": row["overall_survival_event"],
                }
            )
    return rows


def build_figure4_stage_stratified_km_rows(derivation_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for row in derivation_rows:
        subtype = row.get("molecular_subtype_label", "")
        if subtype not in {"EMT", "MSI"}:
            continue
        stage_group = normalize_stage_group(row.get("pathologic_stage", ""))
        if stage_group != "Stage III":
            continue
        rows.append(
            {
                "dataset_id": row.get("dataset_id", "GSE62254"),
                "sample_id": row.get("sample_id") or row.get("gsm_id", ""),
                "stage_group": stage_group,
                "subtype": subtype,
                "overall_survival_months": row["overall_survival_months"],
                "overall_survival_event": row["overall_survival_event"],
            }
        )
    return rows


def build_figure5_cross_km_rows(
    natcom_rows: list[dict[str, str]],
    outcome: str,
    dataset_id: str = "GSE26899",
) -> list[dict[str, str]]:
    time_key = "overall_survival_months" if outcome == "overall_survival" else "recurrence_free_survival_months"
    event_key = "overall_survival_event" if outcome == "overall_survival" else "recurrence_event"
    rows: list[dict[str, str]] = []
    for row in natcom_rows:
        subgroup = row.get("subgroup", "")
        adjuvant_binary = row.get("adjuvant_chemotherapy_binary", "")
        if row.get("dataset_id") != dataset_id or subgroup not in {"EP", "MP"} or adjuvant_binary not in {"0", "1"}:
            continue
        rows.append(
            {
                "dataset_id": row["dataset_id"],
                "patient_id": row.get("patient_id", ""),
                "outcome": outcome,
                "subgroup": subgroup,
                "adjuvant_chemotherapy_binary": adjuvant_binary,
                "overall_survival_months": row[time_key],
                "overall_survival_event": row[event_key],
            }
        )
    return rows


def build_figure3_summary_rows(replication_anchor_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    real_rows = [row for row in replication_anchor_rows if row.get("effect", "") and row["effect"] != "tbd"]
    if not real_rows:
        return [
            {
                "dataset_id": "summary_diamond",
                "outcome": "overall_survival",
                "comparison_subtype": "EMT",
                "reference_subtype": "MSI",
                "effect": "tbd",
                "ci_lower": "tbd",
                "ci_upper": "tbd",
                "model": "fixed_effect_summary",
                "covariates": "none",
                "n": "0",
                "events": "0",
                "source_dataset_count": "0",
            }
        ]

    log_effects = [math.log(float(row["effect"])) for row in real_rows]
    standard_errors = [
        (math.log(float(row["ci_upper"])) - math.log(float(row["ci_lower"]))) / (2 * 1.96)
        for row in real_rows
    ]
    weights = [1.0 / (standard_error**2) for standard_error in standard_errors]
    pooled_log_effect = sum(weight * log_effect for weight, log_effect in zip(weights, log_effects, strict=True)) / sum(weights)
    pooled_standard_error = math.sqrt(1.0 / sum(weights))

    return [
        {
            "dataset_id": "summary_diamond",
            "outcome": real_rows[0].get("outcome", "overall_survival"),
            "comparison_subtype": real_rows[0].get("comparison_subtype", "EMT"),
            "reference_subtype": real_rows[0].get("reference_subtype", "MSI"),
            "effect": format_number(math.exp(pooled_log_effect)),
            "ci_lower": format_number(math.exp(pooled_log_effect - 1.96 * pooled_standard_error)),
            "ci_upper": format_number(math.exp(pooled_log_effect + 1.96 * pooled_standard_error)),
            "model": "fixed_effect_summary",
            "covariates": "none",
            "n": str(sum(int(row["n"]) for row in real_rows)),
            "events": str(sum(int(row["events"]) for row in real_rows)),
            "source_dataset_count": str(len(real_rows)),
        }
    ]


def build_figure5_family_summary_rows(
    pooled_rows: list[dict[str, str]],
    meta_rows: list[dict[str, str]],
) -> list[dict[str, str]]:
    meta_by_outcome = {row["outcome"]: row for row in meta_rows}
    rows: list[dict[str, str]] = []
    for pooled_row in pooled_rows:
        meta_row = meta_by_outcome.get(pooled_row["outcome"])
        if meta_row is None:
            continue
        rows.append(
            {
                "claim_id": pooled_row["claim_id"],
                "dataset_id": pooled_row["dataset_id"],
                "outcome": pooled_row["outcome"],
                "pooled_model": pooled_row["model"],
                "pooled_effect": pooled_row["effect"],
                "pooled_ci_lower": pooled_row["ci_lower"],
                "pooled_ci_upper": pooled_row["ci_upper"],
                "pooled_pvalue": pooled_row["pvalue"],
                "pooled_n": pooled_row["n"],
                "pooled_events": pooled_row["events"],
                "meta_model": meta_row["meta_model"],
                "meta_dataset_count": meta_row["dataset_count"],
                "meta_pooled_hr": meta_row["pooled_hr"],
                "meta_ci_lower": meta_row["ci_lower"],
                "meta_ci_upper": meta_row["ci_upper"],
                "meta_pvalue": meta_row["pooled_pvalue"],
                "cochran_q": meta_row["cochran_q"],
                "heterogeneity_pvalue": meta_row["heterogeneity_pvalue"],
                "i_squared": meta_row["i_squared"],
            }
        )
    return rows


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build panel-ready submission support tables from workflow outputs")
    parser.add_argument("--dataset-summary-input", required=True)
    parser.add_argument("--projection-qc-input", required=True)
    parser.add_argument("--projected-subtypes-input", required=True)
    parser.add_argument("--derivation-clinical-input", required=True)
    parser.add_argument("--derivation-matrix-input", required=True)
    parser.add_argument("--gse15459-outcomes-input", required=True)
    parser.add_argument("--source-data-input", required=True)
    parser.add_argument("--gse84437-clinical-input", required=True)
    parser.add_argument("--gse26899-input", required=True)
    parser.add_argument("--figure3-anchor-input", required=True)
    parser.add_argument("--figure5-pooled-input", required=True)
    parser.add_argument("--figure5-meta-input", required=True)
    parser.add_argument("--figure1-cohort-map-output", required=True)
    parser.add_argument("--figure1-projection-output", required=True)
    parser.add_argument("--figure1-projection-space-output", required=True)
    parser.add_argument("--figure2-km-output", required=True)
    parser.add_argument("--figure2-marker-output", required=True)
    parser.add_argument("--figure3-summary-output", required=True)
    parser.add_argument("--figure3-km-output", required=True)
    parser.add_argument("--figure4-stage-km-output", required=True)
    parser.add_argument("--figure5-family-output", required=True)
    parser.add_argument("--figure5-cross-km-output", required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    dataset_summary_input = Path(args.dataset_summary_input)
    project_root = dataset_summary_input.parents[1]

    dataset_rows = read_tsv(dataset_summary_input)
    projection_rows = read_tsv(Path(args.projection_qc_input))
    projected_subtype_rows = read_tsv(Path(args.projected_subtypes_input))
    derivation_rows = read_tsv(Path(args.derivation_clinical_input))
    gse15459_rows = read_tsv(Path(args.gse15459_outcomes_input))
    gse84437_rows = read_tsv(Path(args.gse84437_clinical_input))
    gse26899_rows = read_tsv(Path(args.gse26899_input))
    figure3_anchor_rows = read_tsv(Path(args.figure3_anchor_input))
    pooled_rows = read_tsv(Path(args.figure5_pooled_input))
    meta_rows = read_tsv(Path(args.figure5_meta_input))
    gse15459_acrg_rows = load_gse15459_acrg_rows(gse15459_rows, Path(args.source_data_input))
    projected_survival_rows = build_projected_survival_rows(projected_subtype_rows, gse84437_rows)

    write_tsv(Path(args.figure1_cohort_map_output), build_figure1_cohort_map_rows(dataset_rows))
    write_tsv(Path(args.figure1_projection_output), build_figure1_projection_rows(projection_rows))
    write_tsv(Path(args.figure1_projection_space_output), build_figure1_projection_space_rows(projected_subtype_rows))
    write_tsv(Path(args.figure2_km_output), build_figure2_derivation_km_rows(derivation_rows))
    write_tsv(
        Path(args.figure2_marker_output),
        build_figure2_marker_expression_rows(
            Path(args.derivation_clinical_input),
            Path(args.derivation_matrix_input),
            project_root,
        ),
    )
    write_tsv(Path(args.figure3_summary_output), build_figure3_summary_rows(figure3_anchor_rows))
    write_tsv(Path(args.figure3_km_output), build_figure3_external_km_rows(gse15459_acrg_rows, projected_survival_rows))
    write_tsv(Path(args.figure4_stage_km_output), build_figure4_stage_stratified_km_rows(build_figure2_derivation_km_rows(derivation_rows)))
    write_tsv(Path(args.figure5_family_output), build_figure5_family_summary_rows(pooled_rows, meta_rows))
    write_tsv(Path(args.figure5_cross_km_output), build_figure5_cross_km_rows(gse26899_rows, outcome="overall_survival"))


if __name__ == "__main__":
    main()
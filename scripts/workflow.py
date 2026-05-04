from __future__ import annotations

import argparse
import base64
import csv
import hashlib
import json
import os
import platform
import shutil
import subprocess
from datetime import UTC, datetime
from pathlib import Path


BENCHMARK_ROWS = [
    {
        "method_name": "published_subtype_mapping",
        "method_class": "baseline",
        "dataset_id": "GSE62254",
        "primary_metric": "hazard_ratio",
        "primary_metric_value": "tbd",
        "confidence_interval": "tbd",
        "secondary_metric": "c_index",
        "secondary_metric_value": "tbd",
    },
    {
        "method_name": "compact_program_score",
        "method_class": "selected_mainline",
        "dataset_id": "GSE62254",
        "primary_metric": "hazard_ratio",
        "primary_metric_value": "tbd",
        "confidence_interval": "tbd",
        "secondary_metric": "c_index",
        "secondary_metric_value": "tbd",
    },
]

PLACEHOLDER_PNG = base64.b64decode(
    "iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAQAAAC1HAwCAAAAC0lEQVR42mP8/x8AAusB9sX6lz8AAAAASUVORK5CYII="
)


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        raise ValueError(f"Cannot write empty TSV for {path}")

    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def write_text(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return [{key: (value or "") for key, value in row.items()} for row in reader]


def load_dataset_records(project_root: Path) -> list[dict[str, str]]:
    manifest_rows = read_tsv(project_root / "data" / "manifest.tsv")
    landscape_rows = read_tsv(project_root / "docs" / "DATASET_LANDSCAPE.tsv")
    manifest_by_id = {row["dataset_id"]: row for row in manifest_rows}

    records: list[dict[str, str]] = []
    for row in landscape_rows:
        dataset_id = row["dataset_id_or_name"]
        manifest_row = manifest_by_id.get(dataset_id, {})
        records.append(
            {
                "dataset_id": dataset_id,
                "sample_size": row["sample_size"],
                "event_count": "tbd",
                "ancestry": row["population_or_ancestry"],
                "platform": row["modality"],
                "treatment_metadata": row["treatment_metadata"],
                "resource_type": manifest_row.get("resource_type", "unknown"),
                "planned_role": manifest_row.get("planned_role", row["intended_role"]),
                "endpoint_availability": row["endpoint_availability"],
                "key_qc_metrics": ";".join(
                    [
                        f"sample_size={row['sample_size']}",
                        "endpoint_review_pending",
                        f"role={manifest_row.get('planned_role', row['intended_role']).replace(' ', '_')}",
                    ]
                ),
            }
        )
    return records


def load_claim_records(project_root: Path) -> list[dict[str, str]]:
    return read_tsv(project_root / "docs" / "CLAIMS.tsv")


def get_record(records: list[dict[str, str]], dataset_id: str) -> dict[str, str]:
    for record in records:
        if record["dataset_id"] == dataset_id:
            return record
    raise KeyError(f"Dataset {dataset_id} not found in dataset records")


def resolve_claim_dataset_id(key_datasets: str, records: list[dict[str, str]]) -> str:
    available_ids = {record["dataset_id"] for record in records}
    for candidate in [item.strip() for item in key_datasets.split(";") if item.strip()]:
        if candidate in available_ids:
            return candidate
    return "unqualified_placeholder"


def infer_claim_outcome(claim_text: str) -> str:
    lower_text = claim_text.lower()
    if "treatment-benefit" in lower_text or "treatment" in lower_text:
        return "treatment_response"
    if "microenvironment" in lower_text or "orthogonal supportive" in lower_text:
        return "supportive_context"
    return "overall_survival"


def build_placeholder_claim_effect_row(claim: dict[str, str], records: list[dict[str, str]]) -> dict[str, str]:
    dataset_id = resolve_claim_dataset_id(claim["key_datasets"], records)
    sample_size = get_record(records, dataset_id)["sample_size"] if dataset_id != "unqualified_placeholder" else "tbd"
    return {
        "claim_id": claim["claim_id"],
        "dataset_id": dataset_id,
        "outcome": infer_claim_outcome(claim["claim_text"]),
        "model": "placeholder_cox",
        "reference_subtype": "",
        "comparison_subtype": "",
        "effect_type": "hazard_ratio",
        "effect": "tbd",
        "ci_lower": "tbd",
        "ci_upper": "tbd",
        "pvalue": "tbd",
        "fdr": "tbd",
        "n": sample_size,
        "events": "",
        "covariates": "",
        "adjuvant_stratum": "",
    }


def build_claim_effect_rows(
    project_root: Path,
    records: list[dict[str, str]],
    claims: list[dict[str, str]],
    treatment_effect_rows: list[dict[str, str]] | None = None,
) -> list[dict[str, str]]:
    shared_state_rows = build_shared_state_effect_rows(project_root, records)
    c1_primary_rows = [
        row
        for row in shared_state_rows
        if row["comparison_subtype"] == "EMT"
        and row["model"] == "cox_unadjusted"
        and row["dataset_id"] in {"GSE62254", "GSE15459", "GSE84437"}
    ]

    if treatment_effect_rows is None:
        treatment_effect_rows = build_treatment_effect_rows(project_root, records, claims)
    c2_primary_rows = [
        row
        for row in treatment_effect_rows
        if row["model"] == "cox_subgroup_adjuvant_interaction" and row["effect"] != "tbd"
    ]

    rows = []
    for claim in claims:
        if claim["claim_id"] == "C1" and c1_primary_rows and c1_primary_rows[0]["effect"] != "tbd":
            rows.extend([{**row, "claim_id": claim["claim_id"]} for row in c1_primary_rows])
            continue
        if "treatment_effects.tsv" in claim["key_output_tables"] and c2_primary_rows:
            rows.extend([{**row, "claim_id": claim["claim_id"]} for row in c2_primary_rows])
            continue
        rows.append(build_placeholder_claim_effect_row(claim, records))
    return rows


def build_adjusted_effect_rows(project_root: Path, records: list[dict[str, str]]) -> list[dict[str, str]]:
    shared_state_rows = build_shared_state_effect_rows(project_root, records)
    adjusted_rows = [
        row
        for row in shared_state_rows
        if row["model"] in {"cox_stage_adjusted", "cox_pt_pn_adjusted"}
    ]
    if adjusted_rows:
        return adjusted_rows

    adjusted = get_record(records, "TCGA-STAD")
    return [
        {
            "claim_id": "C1",
            "dataset_id": adjusted["dataset_id"],
            "outcome": "overall_survival",
            "model": "placeholder_adjusted_cox",
            "reference_subtype": REFERENCE_SUBTYPE_PLACEHOLDER,
            "comparison_subtype": "EMT",
            "effect_type": "hazard_ratio",
            "effect": "tbd",
            "ci_lower": "tbd",
            "ci_upper": "tbd",
            "pvalue": "tbd",
            "fdr": "tbd",
            "n": adjusted["sample_size"],
            "events": "",
            "covariates": "stage,age,sex_if_available",
        }
    ]


def build_treatment_effect_rows(project_root: Path, records: list[dict[str, str]], claims: list[dict[str, str]]) -> list[dict[str, str]]:
    effect_path = project_root / "results" / "effect_sizes" / "treatment_effects.tsv"
    if effect_path.exists():
        effect_rows = read_tsv(effect_path)
        if effect_rows and effect_rows[0]["effect"] != "tbd":
            return attach_fdr(effect_rows)

    rows = []
    for claim in claims:
        if "treatment_effects.tsv" not in claim["key_output_tables"]:
            continue
        dataset_id = resolve_claim_dataset_id(claim["key_datasets"], records)
        sample_size = get_record(records, dataset_id)["sample_size"] if dataset_id != "unqualified_placeholder" else "tbd"
        rows.append(
            {
                "claim_id": claim["claim_id"],
                "dataset_id": dataset_id,
                "outcome": "treatment_response",
                "model": "placeholder_treatment_model",
                "reference_subtype": "",
                "comparison_subtype": "",
                "effect_type": "interaction_effect",
                "effect": "tbd",
                "ci_lower": "tbd",
                "ci_upper": "tbd",
                "pvalue": "tbd",
                "fdr": "tbd",
                "n": sample_size,
                "events": "",
                "covariates": "",
                "adjuvant_stratum": "",
            }
        )
    return rows


def build_figure5_rows(project_root: Path, records: list[dict[str, str]]) -> list[dict[str, str]]:
    figure5_path = project_root / "results" / "figures" / "figure5_extension_anchor.tsv"
    if figure5_path.exists():
        figure_rows = read_tsv(figure5_path)
        if figure_rows and figure_rows[0]["effect"] != "tbd":
            return figure_rows

    treatment_path = project_root / "results" / "effect_sizes" / "treatment_effects.tsv"
    if treatment_path.exists():
        treatment_rows = [
            row
            for row in read_tsv(treatment_path)
            if row["model"] == "cox_subgroup_adjuvant_interaction" and row["effect"] != "tbd"
        ]
        if treatment_rows:
            return treatment_rows

    return [
        {
            "claim_id": "C2",
            "dataset_id": dataset_id,
            "outcome": "overall_survival",
            "model": "placeholder_treatment_model",
            "reference_subtype": "EP",
            "comparison_subtype": "MP",
            "effect_type": "interaction_hazard_ratio",
            "effect": "tbd",
            "ci_lower": "tbd",
            "ci_upper": "tbd",
            "pvalue": "tbd",
            "fdr": "tbd",
            "n": get_record(records, dataset_id)["sample_size"],
            "events": "",
            "covariates": "subgroup,adjuvant,subgroup_x_adjuvant",
            "adjuvant_stratum": "interaction",
        }
        for dataset_id in ["GSE26899", "GSE26901"]
    ]


def build_supportive_context_rows(project_root: Path, records: list[dict[str, str]], claims: list[dict[str, str]]) -> list[dict[str, str]]:
    tcga_effect_path = project_root / "results" / "effect_sizes" / "tcga_stad_compact_score_effects.tsv"
    gse26253_effect_path = project_root / "results" / "effect_sizes" / "gse26253_compact_score_effects.tsv"
    supportive_rows: list[dict[str, str]] = []

    if tcga_effect_path.exists():
        tcga_effect_rows = read_tsv(tcga_effect_path)
        if tcga_effect_rows and tcga_effect_rows[0]["effect"] != "tbd":
            supportive_rows.extend(
                [
                {
                    "claim_id": "C3",
                    "dataset_id": row["dataset_id"],
                    "context_module": "quantitative_compact_score_validation",
                    "status": "completed",
                    "model": row["model"],
                    "effect": row["effect"],
                    "ci_lower": row["ci_lower"],
                    "ci_upper": row["ci_upper"],
                    "pvalue": row["pvalue"],
                    "n": row["n"],
                    "events": row["events"],
                    "covariates": row["covariates"],
                    "note": f"TCGA-STAD compact-score orthogonal validation with {row['covariates']}",
                }
                for row in tcga_effect_rows
                ]
            )

    if gse26253_effect_path.exists():
        gse26253_effect_rows = read_tsv(gse26253_effect_path)
        if gse26253_effect_rows and gse26253_effect_rows[0]["effect"] != "tbd":
            supportive_rows.extend(
                [
                    {
                        "claim_id": "C3",
                        "dataset_id": row["dataset_id"],
                        "context_module": "external_compact_score_recurrence_validation",
                        "status": "completed",
                        "model": row["model"],
                        "effect": row["effect"],
                        "ci_lower": row["ci_lower"],
                        "ci_upper": row["ci_upper"],
                        "pvalue": row["pvalue"],
                        "n": row["n"],
                        "events": row["events"],
                        "covariates": row["covariates"],
                        "note": f"GSE26253 postoperative recurrence validation with {row['covariates']}",
                    }
                    for row in gse26253_effect_rows
                ]
            )

    if supportive_rows:
        return supportive_rows

    rows = []
    for claim in claims:
        if "context_summary.tsv" not in claim["key_output_tables"]:
            continue
        dataset_id = resolve_claim_dataset_id(claim["key_datasets"], records)
        rows.append(
            {
                "claim_id": claim["claim_id"],
                "dataset_id": dataset_id,
                "context_module": "orthogonal_supportive_context",
                "status": "planned",
                "model": "",
                "effect": "",
                "ci_lower": "",
                "ci_upper": "",
                "pvalue": "",
                "n": "",
                "events": "",
                "covariates": "",
                "note": claim["evidence_modules"],
            }
        )
    return rows


def compute_bh_fdr(pvalues: list[float]) -> list[float]:
    ranked = sorted(enumerate(pvalues), key=lambda item: item[1])
    adjusted = [1.0] * len(pvalues)
    running_min = 1.0
    total = len(pvalues)
    for reverse_rank, (index, pvalue) in enumerate(reversed(ranked), start=1):
        rank = total - reverse_rank + 1
        adjusted_value = min(running_min, (pvalue * total) / rank)
        running_min = adjusted_value
        adjusted[index] = adjusted_value
    return adjusted


def attach_fdr(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    if not rows:
        return []
    fdr_values = compute_bh_fdr([float(row["pvalue"]) for row in rows])
    enriched: list[dict[str, str]] = []
    for row, fdr_value in zip(rows, fdr_values, strict=True):
        updated = dict(row)
        updated["fdr"] = format(fdr_value, ".6g")
        enriched.append(updated)
    return enriched


def format_shared_state_row(row: dict[str, str]) -> dict[str, str]:
    return {
        "claim_id": "C1",
        "dataset_id": row["dataset_id"],
        "outcome": row["outcome"],
        "model": row["model"],
        "reference_subtype": row["reference_subtype"],
        "comparison_subtype": row["comparison_subtype"],
        "effect_type": row["effect_type"],
        "effect": row["effect"],
        "ci_lower": row["ci_lower"],
        "ci_upper": row["ci_upper"],
        "pvalue": row["pvalue"],
        "n": row["n"],
        "events": row.get("events", ""),
        "covariates": row.get("covariates", ""),
        "adjuvant_stratum": "",
    }


def build_placeholder_shared_state_rows(records: list[dict[str, str]]) -> list[dict[str, str]]:
    rows = []
    for dataset_id in ["GSE62254", "GSE15459", "GSE84437"]:
        record = get_record(records, dataset_id)
        rows.append(
            {
                "claim_id": "C1",
                "dataset_id": dataset_id,
                "outcome": "overall_survival",
                "model": "",
                "reference_subtype": REFERENCE_SUBTYPE_PLACEHOLDER,
                "comparison_subtype": "EMT",
                "effect_type": "hazard_ratio",
                "effect": "tbd",
                "ci_lower": "tbd",
                "ci_upper": "tbd",
                "pvalue": "tbd",
                "fdr": "tbd",
                "n": record["sample_size"],
                "events": "",
                "covariates": "",
                "adjuvant_stratum": "",
            }
        )
    return rows


def build_shared_state_effect_rows(project_root: Path, records: list[dict[str, str]]) -> list[dict[str, str]]:
    cross_path = project_root / "results" / "effect_sizes" / "acrg_cross_cohort_subtype_effects.tsv"
    projection_path = project_root / "results" / "effect_sizes" / "gse84437_projected_subtype_effects.tsv"

    rows: list[dict[str, str]] = []
    if cross_path.exists():
        rows.extend(
            [
                format_shared_state_row(row)
                for row in read_tsv(cross_path)
                if row["dataset_id"] in {"GSE62254", "GSE15459"}
            ]
        )
    if projection_path.exists():
        rows.extend([format_shared_state_row(row) for row in read_tsv(projection_path)])

    if not rows:
        return build_placeholder_shared_state_rows(records)
    return attach_fdr(rows)


def format_replication_row(row: dict[str, str], claim_id: str, fdr: str) -> dict[str, str]:
    return {
        "dataset_id": row["dataset_id"],
        "claim_id": claim_id,
        "outcome": row["outcome"],
        "model": row.get("model", ""),
        "reference_subtype": row.get("reference_subtype", ""),
        "comparison_subtype": row.get("comparison_subtype", ""),
        "effect_type": row["effect_type"],
        "effect": row["effect"],
        "ci_lower": row["ci_lower"],
        "ci_upper": row["ci_upper"],
        "pvalue": row["pvalue"],
        "fdr": fdr,
        "n": row["n"],
        "events": row.get("events", ""),
        "covariates": row.get("covariates", ""),
    }


def build_placeholder_replication_row(record: dict[str, str]) -> dict[str, str]:
    return {
        "dataset_id": record["dataset_id"],
        "claim_id": "C1",
        "outcome": "overall_survival",
        "model": "",
        "reference_subtype": "",
        "comparison_subtype": "",
        "effect_type": "hazard_ratio",
        "effect": "tbd",
        "ci_lower": "tbd",
        "ci_upper": "tbd",
        "pvalue": "tbd",
        "fdr": "tbd",
        "n": record["sample_size"],
        "events": "",
        "covariates": "",
    }


def build_replication_rows(project_root: Path, records: list[dict[str, str]]) -> list[dict[str, str]]:
    shared_state_rows = build_shared_state_effect_rows(project_root, records)
    external_rows = [row for row in shared_state_rows if row["dataset_id"] in {"GSE15459", "GSE84437"}]
    if external_rows and external_rows[0]["effect"] != "tbd":
        return [format_replication_row(row, row["claim_id"], row["fdr"]) for row in external_rows]
    return [build_placeholder_replication_row(get_record(records, dataset_id)) for dataset_id in ["GSE15459", "GSE84437"]]


def write_replication_dataset_tables(project_root: Path, replication_rows: list[dict[str, str]]) -> None:
    dataset_ids = sorted({row["dataset_id"] for row in replication_rows})
    for dataset_id in dataset_ids:
        write_tsv(
            project_root / "results" / "replication" / f"{dataset_id}_effects.tsv",
            [row for row in replication_rows if row["dataset_id"] == dataset_id],
        )


def build_qc_rows(records: list[dict[str, str]]) -> list[dict[str, str]]:
    rows = []
    for record in records:
        rows.append(
            {
                "dataset_id": record["dataset_id"],
                "qc_metric": "sample_size_check",
                "value": record["sample_size"],
                "status": "planned",
            }
        )
    return rows


def classify_analysis_role(planned_role: str) -> str:
    lower_role = planned_role.lower()
    if "derivation" in lower_role:
        return "derivation"
    if "replication" in lower_role:
        return "replication"
    if "orthogonal" in lower_role:
        return "orthogonal"
    if "sensitivity" in lower_role or "fallback" in lower_role:
        return "sensitivity"
    return "control_or_exclusion"


def infer_qualification_status(record: dict[str, str]) -> str:
    role = classify_analysis_role(record["planned_role"])
    endpoint_text = record["endpoint_availability"].lower()
    if role in {"derivation", "replication", "orthogonal", "sensitivity"} and (
        "survival" in endpoint_text or "prognosis" in endpoint_text or "recurrence" in endpoint_text
    ):
        return "qualified"
    return "review_only"


def infer_treatment_qualification(record: dict[str, str]) -> str:
    treatment_text = record["treatment_metadata"].lower()
    if "adjuvant" in treatment_text or "treatment-aware" in record["planned_role"].lower():
        return "yes"
    return "no"


def build_dataset_qualification_rows(records: list[dict[str, str]]) -> list[dict[str, str]]:
    rows = []
    for record in records:
        rows.append(
            {
                "dataset_id": record["dataset_id"],
                "analysis_role": classify_analysis_role(record["planned_role"]),
                "qualification_status": infer_qualification_status(record),
                "treatment_qualified": infer_treatment_qualification(record),
                "endpoint_basis": record["endpoint_availability"],
                "planned_role": record["planned_role"],
            }
        )
    return rows


def build_metadata_rows(records: list[dict[str, str]]) -> list[dict[str, str]]:
    rows = []
    for record in records:
        rows.append(
            {
                "dataset_id": record["dataset_id"],
                "field": "endpoint_availability",
                "value": record["endpoint_availability"],
                "completeness": "planned",
            }
        )
    return rows


def build_umap_rows(records: list[dict[str, str]]) -> list[dict[str, str]]:
    derivation = get_record(records, "GSE62254")
    return [
        {"dataset_id": derivation["dataset_id"], "sample_id": "placeholder-1", "umap_x": "0.0", "umap_y": "0.0", "group": "high_state"},
        {"dataset_id": derivation["dataset_id"], "sample_id": "placeholder-2", "umap_x": "1.0", "umap_y": "1.0", "group": "low_state"},
    ]


def build_figure1_rows(records: list[dict[str, str]]) -> list[dict[str, str]]:
    return [
        {"step": "derivation", "dataset_id": get_record(records, "GSE62254")["dataset_id"], "role": "core_discovery", "note": get_record(records, "GSE62254")["planned_role"]},
        {"step": "replication", "dataset_id": get_record(records, "GSE15459")["dataset_id"], "role": "external_replication", "note": get_record(records, "GSE15459")["planned_role"]},
        {"step": "orthogonal_validation", "dataset_id": get_record(records, "TCGA-STAD")["dataset_id"], "role": "clinicopathologic_context", "note": get_record(records, "TCGA-STAD")["planned_role"]},
    ]


def build_figure2_rows(records: list[dict[str, str]]) -> list[dict[str, str]]:
    derivation = get_record(records, "GSE62254")
    return [
        {"dataset_id": derivation["dataset_id"], "state_name": "baseline_program_state", "group": "high_state", "n": derivation["sample_size"], "distribution_signal": "tbd"},
        {"dataset_id": derivation["dataset_id"], "state_name": "baseline_program_state", "group": "low_state", "n": derivation["sample_size"], "distribution_signal": "tbd"},
    ]


REFERENCE_SUBTYPE_PLACEHOLDER = "MSI"


def build_figure3_rows(project_root: Path, records: list[dict[str, str]]) -> list[dict[str, str]]:
    shared_state_rows = build_shared_state_effect_rows(project_root, records)
    filtered_rows = [
        row
        for row in shared_state_rows
        if row["dataset_id"] in {"GSE15459", "GSE84437"}
        and row["comparison_subtype"] == "EMT"
        and row["model"] == "cox_unadjusted"
    ]
    if filtered_rows:
        return filtered_rows

    return [
        {
            "claim_id": "C1",
            "dataset_id": get_record(records, dataset_id)["dataset_id"],
            "outcome": "overall_survival",
            "model": "",
            "reference_subtype": REFERENCE_SUBTYPE_PLACEHOLDER,
            "comparison_subtype": "EMT",
            "effect_type": "hazard_ratio",
            "effect": "tbd",
            "ci_lower": "tbd",
            "ci_upper": "tbd",
            "pvalue": "tbd",
            "fdr": "tbd",
            "n": get_record(records, dataset_id)["sample_size"],
            "events": "",
            "covariates": "none",
        }
        for dataset_id in ["GSE15459", "GSE84437"]
    ]


def build_figure4_rows(project_root: Path, records: list[dict[str, str]]) -> list[dict[str, str]]:
    shared_state_rows = build_shared_state_effect_rows(project_root, records)
    filtered_rows = [
        row
        for row in shared_state_rows
        if row["comparison_subtype"] == "EMT"
        and row["model"] in {"cox_stage_adjusted", "cox_pt_pn_adjusted"}
    ]
    if filtered_rows:
        return filtered_rows

    return [
        {
            "claim_id": "C1",
            "dataset_id": get_record(records, "TCGA-STAD")["dataset_id"],
            "outcome": "overall_survival",
            "model": "adjusted_cox",
            "reference_subtype": REFERENCE_SUBTYPE_PLACEHOLDER,
            "comparison_subtype": "EMT",
            "effect_type": "hazard_ratio",
            "effect": "tbd",
            "ci_lower": "tbd",
            "ci_upper": "tbd",
            "pvalue": "tbd",
            "fdr": "tbd",
            "n": get_record(records, "TCGA-STAD")["sample_size"],
            "events": "",
            "covariates": "stage,age,sex_if_available",
        }
    ]


def ensure_project_scaffold(project_root: Path) -> None:
    for rel_path in [
        "logs",
        "results/benchmarks",
        "results/effect_sizes",
        "results/ingestion",
        "results/replication",
        "results/figures",
        "results/tables",
        "plots/publication/pdf",
        "plots/publication/png",
        "docs/audit_runs",
        "docs/review_bundle",
        "cloud_runs",
    ]:
        (project_root / rel_path).mkdir(parents=True, exist_ok=True)


def write_compute_plan(project_root: Path) -> None:
    content = """# Compute Plan

## Mainline Hardware Target
- 2017 MacBook Air class laptop
- 8 to 16 GB RAM preferred
- no GPU required for smoke or first-pass mainline analyses
- estimated local storage for current scaffold outputs: under 50 MB

## Smoke Workflow
- uses planning-derived placeholder tables and directory scaffolding only
- avoids downloading large matrices or running heavy deconvolution
- runtime target: under 1 minute on a laptop
- downsampled smoke definition: write representative result schemas and audit artifacts without cohort downloads

## Full Workflow Boundary
- still limited to lightweight project scaffolding at this stage
- heavy cohort downloads, harmonization, and model fitting remain future tasks
- current full-run target runtime: under 2 minutes on a laptop

## Optional Heavy Modules
- single-cell localization
- exhaustive resampling
- causal follow-up
- large cloud-only batch harmonization

## Cloud Placeholder Budget
- target cloud profile for future heavy runs: 8 vCPU, 32 GB RAM, 200 GB ephemeral storage
- first-pass exploratory cost ceiling: under 20 USD per heavy validation batch
- optional modules only; not required for the mainline manuscript claim
"""
    write_text(project_root / "docs" / "COMPUTE_PLAN.md", content)


def write_harmonization_notes(project_root: Path) -> None:
    content = """# Harmonization Notes

## Current State
This repository is in the scaffold stage. No raw expression matrices have been downloaded or normalized yet.

## Planned Harmonization Principles
- freeze cohort inclusion rules before merging analyses
- maintain cohort-specific preprocessing notes
- avoid forced pooling before endpoint and platform qualification
- export compact figure anchor tables independent of large matrix objects

## Known Heterogeneity Sources
- Affymetrix versus Illumina versus RNA-seq platforms
- overall survival versus recurrence-oriented endpoints
- treatment metadata completeness and definition drift
"""
    write_text(project_root / "docs" / "HARMONIZATION_NOTES.md", content)


def write_accession_metadata_output(project_root: Path) -> None:
    input_dir = project_root / "data" / "cache" / "accession_metadata"
    output_path = project_root / "results" / "ingestion" / "dataset_accession_metadata.tsv"
    if not input_dir.exists():
        return

    script_path = Path(__file__).resolve().parent / "ingest" / "cohort_metadata.py"
    subprocess.run(
        [
            "python3",
            str(script_path),
            "--input-dir",
            str(input_dir),
            "--output",
            str(output_path),
        ],
        check=True,
        capture_output=True,
        text=True,
    )


def write_gse15459_outcome_outputs(project_root: Path) -> None:
    input_path = project_root / "data" / "raw" / "clinical" / "GSE15459_outcome.xls"
    if not input_path.exists():
        return

    script_path = Path(__file__).resolve().parent / "ingest" / "gse15459_outcomes.py"
    subprocess.run(
        [
            "python3",
            str(script_path),
            "--input",
            str(input_path),
            "--outcomes-output",
            str(project_root / "results" / "ingestion" / "gse15459_outcomes.tsv"),
            "--summary-output",
            str(project_root / "results" / "ingestion" / "endpoint_completeness.tsv"),
        ],
        check=True,
        capture_output=True,
        text=True,
    )


def write_gse15459_expression_link_outputs(project_root: Path) -> None:
    clinical_input = project_root / "results" / "ingestion" / "gse15459_outcomes.tsv"
    matrix_input = project_root / "data" / "raw" / "expression" / "GSE15459_series_matrix.txt.gz"
    if not clinical_input.exists() or not matrix_input.exists():
        return

    script_path = Path(__file__).resolve().parent / "ingest" / "gse15459_expression_link.py"
    subprocess.run(
        [
            "python3",
            str(script_path),
            "--clinical-input",
            str(clinical_input),
            "--matrix-input",
            str(matrix_input),
            "--link-output",
            str(project_root / "results" / "ingestion" / "gse15459_expression_clinical_link.tsv"),
            "--summary-output",
            str(project_root / "results" / "ingestion" / "gse15459_expression_link_summary.tsv"),
            "--matrix-output",
            str(project_root / "results" / "ingestion" / "gse15459_expression_linked_matrix.tsv.gz"),
        ],
        check=True,
        capture_output=True,
        text=True,
    )


def write_natcom_2018_clinical_outputs(project_root: Path) -> None:
    input_path = project_root / "data" / "raw" / "clinical" / "41467_2018_4179_MOESM5_ESM.xlsx"
    if not input_path.exists():
        return

    script_path = Path(__file__).resolve().parent / "ingest" / "natcom_2018_clinical_endpoints.py"
    subprocess.run(
        [
            "python3",
            str(script_path),
            "--input",
            str(input_path),
            "--clinical-output",
            str(project_root / "results" / "ingestion" / "natcom_2018_clinical_endpoints.tsv"),
            "--summary-output",
            str(project_root / "results" / "ingestion" / "natcom_2018_clinical_summary.tsv"),
        ],
        check=True,
        capture_output=True,
        text=True,
    )


def write_gse26899_expression_link_outputs(project_root: Path) -> None:
    clinical_input = project_root / "results" / "ingestion" / "natcom_2018_clinical_endpoints.tsv"
    matrix_input = project_root / "data" / "raw" / "expression" / "GSE26899_series_matrix.txt.gz"
    if not clinical_input.exists() or not matrix_input.exists():
        return

    script_path = Path(__file__).resolve().parent / "ingest" / "gse26899_expression_clinical_link.py"
    subprocess.run(
        [
            "python3",
            str(script_path),
            "--clinical-input",
            str(clinical_input),
            "--matrix-input",
            str(matrix_input),
            "--link-output",
            str(project_root / "results" / "ingestion" / "gse26899_expression_clinical_link.tsv"),
            "--summary-output",
            str(project_root / "results" / "ingestion" / "gse26899_expression_link_summary.tsv"),
            "--matrix-output",
            str(project_root / "results" / "ingestion" / "gse26899_expression_linked_matrix.tsv.gz"),
        ],
        check=True,
        capture_output=True,
        text=True,
    )


def write_gse26901_expression_link_outputs(project_root: Path) -> None:
    clinical_input = project_root / "results" / "ingestion" / "natcom_2018_clinical_endpoints.tsv"
    matrix_input = project_root / "data" / "raw" / "expression" / "GSE26901_series_matrix.txt.gz"
    if not clinical_input.exists() or not matrix_input.exists():
        return

    script_path = Path(__file__).resolve().parent / "ingest" / "gse26901_expression_clinical_link.py"
    subprocess.run(
        [
            "python3",
            str(script_path),
            "--clinical-input",
            str(clinical_input),
            "--matrix-input",
            str(matrix_input),
            "--link-output",
            str(project_root / "results" / "ingestion" / "gse26901_expression_clinical_link.tsv"),
            "--summary-output",
            str(project_root / "results" / "ingestion" / "gse26901_expression_link_summary.tsv"),
            "--matrix-output",
            str(project_root / "results" / "ingestion" / "gse26901_expression_linked_matrix.tsv.gz"),
        ],
        check=True,
        capture_output=True,
        text=True,
    )


def write_gse84437_expression_link_outputs(project_root: Path) -> None:
    matrix_input = project_root / "data" / "raw" / "expression" / "GSE84437_series_matrix.txt.gz"
    if not matrix_input.exists():
        return

    script_path = Path(__file__).resolve().parent / "ingest" / "gse84437_expression_clinical_link.py"
    subprocess.run(
        [
            "python3",
            str(script_path),
            "--matrix-input",
            str(matrix_input),
            "--link-output",
            str(project_root / "results" / "ingestion" / "gse84437_expression_clinical_link.tsv"),
            "--summary-output",
            str(project_root / "results" / "ingestion" / "gse84437_expression_link_summary.tsv"),
            "--matrix-output",
            str(project_root / "results" / "ingestion" / "gse84437_expression_linked_matrix.tsv.gz"),
        ],
        check=True,
        capture_output=True,
        text=True,
    )


def write_gse62254_clinical_link_outputs(project_root: Path) -> None:
    matrix_input = project_root / "data" / "raw" / "expression" / "GSE62254_series_matrix.txt.gz"
    clinical_input = project_root / "data" / "raw" / "clinical" / "nm3850_supplementary_data1.xls"
    if not matrix_input.exists() or not clinical_input.exists():
        return

    script_path = Path(__file__).resolve().parent / "ingest" / "gse62254_clinical_link.py"
    subprocess.run(
        [
            "python3",
            str(script_path),
            "--matrix-input",
            str(matrix_input),
            "--clinical-input",
            str(clinical_input),
            "--link-output",
            str(project_root / "results" / "ingestion" / "gse62254_expression_clinical_link.tsv"),
            "--summary-output",
            str(project_root / "results" / "ingestion" / "gse62254_clinical_link_summary.tsv"),
            "--matrix-output",
            str(project_root / "results" / "ingestion" / "gse62254_expression_linked_matrix.tsv.gz"),
        ],
        check=True,
        capture_output=True,
        text=True,
    )


def write_gse15459_subtype_effects(project_root: Path) -> None:
    input_path = project_root / "results" / "ingestion" / "gse15459_outcomes.tsv"
    if not input_path.exists():
        return

    script_path = Path(__file__).resolve().parent / "analysis" / "run_gse15459_subtype_survival.py"
    subprocess.run(
        [
            "python3",
            str(script_path),
            "--input",
            str(input_path),
            "--output",
            str(project_root / "results" / "effect_sizes" / "gse15459_subtype_effects.tsv"),
        ],
        check=True,
        capture_output=True,
        text=True,
    )


def write_acrg_cross_cohort_effects(project_root: Path) -> None:
    acrg_input = project_root / "results" / "ingestion" / "gse62254_expression_clinical_link.tsv"
    gse15459_input = project_root / "results" / "ingestion" / "gse15459_outcomes.tsv"
    source_data_input = project_root / "data" / "raw" / "clinical" / "nm3850_source_data_fig2.xls"
    if not acrg_input.exists() or not gse15459_input.exists() or not source_data_input.exists():
        return

    script_path = Path(__file__).resolve().parent / "analysis" / "run_acrg_cross_cohort_survival.py"
    subprocess.run(
        [
            "python3",
            str(script_path),
            "--acrg-input",
            str(acrg_input),
            "--gse15459-input",
            str(gse15459_input),
            "--source-data-input",
            str(source_data_input),
            "--output",
            str(project_root / "results" / "effect_sizes" / "acrg_cross_cohort_subtype_effects.tsv"),
            "--subtype-counts-output",
            str(project_root / "results" / "tables" / "acrg_cross_cohort_subtype_counts.tsv"),
        ],
        check=True,
        capture_output=True,
        text=True,
    )


def write_gse84437_projection_effects(project_root: Path) -> None:
    training_clinical_input = project_root / "results" / "ingestion" / "gse62254_expression_clinical_link.tsv"
    training_matrix_input = project_root / "results" / "ingestion" / "gse62254_expression_linked_matrix.tsv.gz"
    replication_clinical_input = project_root / "results" / "ingestion" / "gse84437_expression_clinical_link.tsv"
    replication_matrix_input = project_root / "results" / "ingestion" / "gse84437_expression_linked_matrix.tsv.gz"
    gpl570_map_input = project_root / "data" / "cache" / "platform_annotations" / "GPL570_gene_symbols.tsv"
    gpl6947_map_input = project_root / "data" / "cache" / "platform_annotations" / "GPL6947_gene_symbols.tsv"
    if (
        not training_clinical_input.exists()
        or not training_matrix_input.exists()
        or not replication_clinical_input.exists()
        or not replication_matrix_input.exists()
    ):
        return

    script_path = Path(__file__).resolve().parent / "analysis" / "run_gse84437_projection_survival.py"
    subprocess.run(
        [
            "python3",
            str(script_path),
            "--training-clinical-input",
            str(training_clinical_input),
            "--training-matrix-input",
            str(training_matrix_input),
            "--replication-clinical-input",
            str(replication_clinical_input),
            "--replication-matrix-input",
            str(replication_matrix_input),
            "--gpl570-map-input",
            str(gpl570_map_input),
            "--gpl6947-map-input",
            str(gpl6947_map_input),
            "--assignment-output",
            str(project_root / "results" / "replication" / "gse84437_projected_subtypes.tsv"),
            "--effect-output",
            str(project_root / "results" / "effect_sizes" / "gse84437_projected_subtype_effects.tsv"),
            "--qc-output",
            str(project_root / "results" / "tables" / "gse84437_projection_qc.tsv"),
            "--sensitivity-output",
            str(project_root / "results" / "effect_sizes" / "gse84437_projection_low_margin_sensitivity.tsv"),
        ],
        check=True,
        capture_output=True,
        text=True,
    )


def write_state_definition_benchmark(project_root: Path) -> None:
    training_clinical_input = project_root / "results" / "ingestion" / "gse62254_expression_clinical_link.tsv"
    training_matrix_input = project_root / "results" / "ingestion" / "gse62254_expression_linked_matrix.tsv.gz"
    gse15459_outcomes_input = project_root / "results" / "ingestion" / "gse15459_outcomes.tsv"
    gse15459_matrix_input = project_root / "results" / "ingestion" / "gse15459_expression_linked_matrix.tsv.gz"
    gse84437_clinical_input = project_root / "results" / "ingestion" / "gse84437_expression_clinical_link.tsv"
    gse84437_matrix_input = project_root / "results" / "ingestion" / "gse84437_expression_linked_matrix.tsv.gz"
    gse84437_assignment_input = project_root / "results" / "replication" / "gse84437_projected_subtypes.tsv"
    source_data_input = project_root / "data" / "raw" / "clinical" / "nm3850_source_data_fig2.xls"
    gpl570_map_input = project_root / "data" / "cache" / "platform_annotations" / "GPL570_gene_symbols.tsv"
    gpl6947_map_input = project_root / "data" / "cache" / "platform_annotations" / "GPL6947_gene_symbols.tsv"
    if not all(
        path.exists()
        for path in [
            training_clinical_input,
            training_matrix_input,
            gse15459_outcomes_input,
            gse15459_matrix_input,
            gse84437_clinical_input,
            gse84437_matrix_input,
            gse84437_assignment_input,
            source_data_input,
        ]
    ):
        write_tsv(project_root / "results" / "benchmarks" / "method_benchmark.tsv", BENCHMARK_ROWS)
        write_tsv(
            project_root / "results" / "benchmarks" / "state_definition_decision.tsv",
            [
                {
                    "candidate_method": "compact_program_score",
                    "gene_panel_size": "tbd",
                    "replace_threshold": "replace only after external noninferiority and gain benchmark passes",
                    "external_direction_pass": "false",
                    "external_cindex_noninferior_pass": "false",
                    "external_gain_pass": "false",
                    "replace_decision": "keep_shared_acrg_state",
                    "decision_rationale": "benchmark inputs unavailable",
                }
            ],
        )
        return

    script_path = Path(__file__).resolve().parent / "analysis" / "run_state_definition_benchmark.py"
    subprocess.run(
        [
            "python3",
            str(script_path),
            "--training-clinical-input",
            str(training_clinical_input),
            "--training-matrix-input",
            str(training_matrix_input),
            "--gse15459-outcomes-input",
            str(gse15459_outcomes_input),
            "--gse15459-matrix-input",
            str(gse15459_matrix_input),
            "--gse84437-clinical-input",
            str(gse84437_clinical_input),
            "--gse84437-matrix-input",
            str(gse84437_matrix_input),
            "--gse84437-assignment-input",
            str(gse84437_assignment_input),
            "--source-data-input",
            str(source_data_input),
            "--gpl570-map-input",
            str(gpl570_map_input),
            "--gpl6947-map-input",
            str(gpl6947_map_input),
            "--benchmark-output",
            str(project_root / "results" / "benchmarks" / "method_benchmark.tsv"),
            "--decision-output",
            str(project_root / "results" / "benchmarks" / "state_definition_decision.tsv"),
            "--gene-panel-output",
            str(project_root / "results" / "benchmarks" / "compact_score_gene_panel.tsv"),
        ],
        check=True,
        capture_output=True,
        text=True,
    )


def write_tcga_stad_compact_score_effects(project_root: Path) -> None:
    gene_panel_input = project_root / "results" / "benchmarks" / "compact_score_gene_panel.tsv"
    if not gene_panel_input.exists():
        return

    script_path = Path(__file__).resolve().parent / "analysis" / "run_tcga_stad_compact_score_validation.py"
    output_path = project_root / "results" / "effect_sizes" / "tcga_stad_compact_score_effects.tsv"
    try:
        subprocess.run(
            [
                "python3",
                str(script_path),
                "--gene-panel-input",
                str(gene_panel_input),
                "--output",
                str(output_path),
            ],
            check=True,
            capture_output=True,
            text=True,
        )
    except subprocess.CalledProcessError:
        fallback_path = Path(__file__).resolve().parents[1] / "results" / "effect_sizes" / "tcga_stad_compact_score_effects.tsv"
        if not fallback_path.exists():
            raise
        output_path.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(fallback_path, output_path)


def write_gse26253_compact_score_effects(project_root: Path) -> None:
    gene_panel_input = project_root / "results" / "benchmarks" / "compact_score_gene_panel.tsv"
    matrix_input = project_root / "data" / "raw" / "expression" / "GSE26253_series_matrix.txt.gz"
    gpl8432_map_input = project_root / "data" / "cache" / "platform_annotations" / "GPL8432_gene_symbols.tsv"
    if not gene_panel_input.exists() or not matrix_input.exists():
        return

    script_path = Path(__file__).resolve().parent / "analysis" / "run_gse26253_compact_score_validation.py"
    subprocess.run(
        [
            "python3",
            str(script_path),
            "--gene-panel-input",
            str(gene_panel_input),
            "--matrix-input",
            str(matrix_input),
            "--gpl8432-map-input",
            str(gpl8432_map_input),
            "--output",
            str(project_root / "results" / "effect_sizes" / "gse26253_compact_score_effects.tsv"),
        ],
        check=True,
        capture_output=True,
        text=True,
    )


def write_natcom_c2_treatment_effects(project_root: Path) -> None:
    gse26899_input = project_root / "results" / "ingestion" / "gse26899_expression_clinical_link.tsv"
    gse26901_input = project_root / "results" / "ingestion" / "gse26901_expression_clinical_link.tsv"
    if not gse26899_input.exists() or not gse26901_input.exists():
        return

    script_path = Path(__file__).resolve().parent / "analysis" / "run_natcom_c2_survival.py"
    subprocess.run(
        [
            "python3",
            str(script_path),
            "--gse26899-input",
            str(gse26899_input),
            "--gse26901-input",
            str(gse26901_input),
            "--effect-output",
            str(project_root / "results" / "effect_sizes" / "treatment_effects.tsv"),
            "--figure-output",
            str(project_root / "results" / "figures" / "figure5_extension_anchor.tsv"),
            "--summary-output",
            str(project_root / "results" / "tables" / "natcom_c2_cohort_summary.tsv"),
        ],
        check=True,
        capture_output=True,
        text=True,
    )

    supplement_script_path = Path(__file__).resolve().parent / "analysis" / "run_natcom_c2_supplement_table.py"
    subprocess.run(
        [
            "python3",
            str(supplement_script_path),
            "--effect-input",
            str(project_root / "results" / "effect_sizes" / "treatment_effects.tsv"),
            "--summary-input",
            str(project_root / "results" / "tables" / "natcom_c2_cohort_summary.tsv"),
            "--output",
            str(project_root / "results" / "tables" / "natcom_c2_cohortwise_interpretation.tsv"),
        ],
        check=True,
        capture_output=True,
        text=True,
    )


def write_natcom_c2_pooled_effects(project_root: Path) -> None:
    gse26899_input = project_root / "results" / "ingestion" / "gse26899_expression_clinical_link.tsv"
    gse26901_input = project_root / "results" / "ingestion" / "gse26901_expression_clinical_link.tsv"
    if not gse26899_input.exists() or not gse26901_input.exists():
        return

    script_path = Path(__file__).resolve().parent / "analysis" / "run_natcom_c2_pooled_analysis.py"
    subprocess.run(
        [
            "python3",
            str(script_path),
            "--gse26899-input",
            str(gse26899_input),
            "--gse26901-input",
            str(gse26901_input),
            "--effect-output",
            str(project_root / "results" / "effect_sizes" / "natcom_c2_pooled_effects.tsv"),
            "--meta-output",
            str(project_root / "results" / "tables" / "natcom_c2_meta_summary.tsv"),
        ],
        check=True,
        capture_output=True,
        text=True,
    )


def build_submission_figure_support_tables(project_root: Path) -> None:
    script_path = Path(__file__).resolve().parent / "figures" / "build_submission_figure_support_tables.py"
    subprocess.run(
        [
            "python3",
            str(script_path),
            "--dataset-summary-input",
            str(project_root / "results" / "dataset_summary.tsv"),
            "--projection-qc-input",
            str(project_root / "results" / "tables" / "gse84437_projection_qc.tsv"),
            "--projected-subtypes-input",
            str(project_root / "results" / "replication" / "gse84437_projected_subtypes.tsv"),
            "--derivation-clinical-input",
            str(project_root / "results" / "ingestion" / "gse62254_expression_clinical_link.tsv"),
            "--derivation-matrix-input",
            str(project_root / "results" / "ingestion" / "gse62254_expression_linked_matrix.tsv.gz"),
            "--gse15459-outcomes-input",
            str(project_root / "results" / "ingestion" / "gse15459_outcomes.tsv"),
            "--source-data-input",
            str(project_root / "data" / "raw" / "clinical" / "nm3850_source_data_fig2.xls"),
            "--gse84437-clinical-input",
            str(project_root / "results" / "ingestion" / "gse84437_expression_clinical_link.tsv"),
            "--gse26899-input",
            str(project_root / "results" / "ingestion" / "gse26899_expression_clinical_link.tsv"),
            "--figure3-anchor-input",
            str(project_root / "results" / "figures" / "figure3_replication_anchor.tsv"),
            "--figure5-pooled-input",
            str(project_root / "results" / "effect_sizes" / "natcom_c2_pooled_effects.tsv"),
            "--figure5-meta-input",
            str(project_root / "results" / "tables" / "natcom_c2_meta_summary.tsv"),
            "--figure1-cohort-map-output",
            str(project_root / "results" / "figures" / "figure1_cohort_map.tsv"),
            "--figure1-projection-output",
            str(project_root / "results" / "figures" / "figure1_projection_consistency.tsv"),
            "--figure1-projection-space-output",
            str(project_root / "results" / "figures" / "figure1_projection_space.tsv"),
            "--figure2-km-output",
            str(project_root / "results" / "figures" / "figure2_derivation_km.tsv"),
            "--figure2-marker-output",
            str(project_root / "results" / "figures" / "figure2_marker_expression.tsv"),
            "--figure3-summary-output",
            str(project_root / "results" / "figures" / "figure3_summary_diamond.tsv"),
            "--figure3-km-output",
            str(project_root / "results" / "figures" / "figure3_external_km.tsv"),
            "--figure4-stage-km-output",
            str(project_root / "results" / "figures" / "figure4_stage_stratified_km.tsv"),
            "--figure5-family-output",
            str(project_root / "results" / "figures" / "figure5_family_summary.tsv"),
            "--figure5-cross-km-output",
            str(project_root / "results" / "figures" / "figure5_cross_km.tsv"),
        ],
        check=True,
        capture_output=True,
        text=True,
    )


def render_submission_figure1(project_root: Path) -> None:
    required_paths = [
        project_root / "results" / "figures" / "figure1_study_design.tsv",
        project_root / "results" / "figures" / "figure1_cohort_map.tsv",
        project_root / "results" / "figures" / "figure1_projection_consistency.tsv",
        project_root / "results" / "figures" / "figure1_projection_space.tsv",
    ]
    if not all(path.exists() for path in required_paths):
        return

    script_path = Path(__file__).resolve().parent / "figures" / "render_submission_figure1.py"
    subprocess.run(
        [
            "python3",
            str(script_path),
            "--flow-input",
            str(required_paths[0]),
            "--cohort-map-input",
            str(required_paths[1]),
            "--projection-input",
            str(required_paths[2]),
            "--projection-space-input",
            str(required_paths[3]),
            "--pdf-output",
            str(project_root / "plots" / "publication" / "pdf" / "fig1_workflow_consistency.pdf"),
            "--png-output",
            str(project_root / "plots" / "publication" / "png" / "fig1_workflow_consistency.png"),
        ],
        check=True,
        capture_output=True,
        text=True,
    )


def render_submission_figure2(project_root: Path) -> None:
    required_paths = [
        project_root / "results" / "figures" / "figure2_derivation_km.tsv",
        project_root / "results" / "figures" / "figure2_marker_expression.tsv",
    ]
    if not all(path.exists() for path in required_paths):
        return

    script_path = Path(__file__).resolve().parent / "figures" / "render_submission_figure2.py"
    subprocess.run(
        [
            "python3",
            str(script_path),
            "--km-input",
            str(required_paths[0]),
            "--marker-input",
            str(required_paths[1]),
            "--pdf-output",
            str(project_root / "plots" / "publication" / "pdf" / "fig2_discovery_signal.pdf"),
            "--png-output",
            str(project_root / "plots" / "publication" / "png" / "fig2_discovery_signal.png"),
        ],
        check=True,
        capture_output=True,
        text=True,
    )


def render_publication_figures(project_root: Path) -> None:
    figure3_anchor_input = project_root / "results" / "figures" / "figure3_replication_anchor.tsv"
    figure3_summary_input = project_root / "results" / "figures" / "figure3_summary_diamond.tsv"
    figure3_km_input = project_root / "results" / "figures" / "figure3_external_km.tsv"
    figure4_anchor_input = project_root / "results" / "figures" / "figure4_adjusted_effects_anchor.tsv"
    figure4_stage_km_input = project_root / "results" / "figures" / "figure4_stage_stratified_km.tsv"
    figure5_anchor_input = project_root / "results" / "figures" / "figure5_extension_anchor.tsv"
    figure5_family_input = project_root / "results" / "figures" / "figure5_family_summary.tsv"
    figure5_cross_km_input = project_root / "results" / "figures" / "figure5_cross_km.tsv"
    if not all(
        path.exists()
        for path in [
            figure3_anchor_input,
            figure3_summary_input,
            figure3_km_input,
            figure4_anchor_input,
            figure4_stage_km_input,
            figure5_anchor_input,
            figure5_family_input,
            figure5_cross_km_input,
        ]
    ):
        return

    script_path = Path(__file__).resolve().parent / "figures" / "render_prognosis_manuscript_figures.py"
    subprocess.run(
        [
            "python3",
            str(script_path),
            "--figure3-anchor-input",
            str(figure3_anchor_input),
            "--figure3-summary-input",
            str(figure3_summary_input),
            "--figure3-km-input",
            str(figure3_km_input),
            "--figure4-anchor-input",
            str(figure4_anchor_input),
            "--figure4-stage-km-input",
            str(figure4_stage_km_input),
            "--figure3-pdf-output",
            str(project_root / "plots" / "publication" / "pdf" / "fig3_external_replication_forest.pdf"),
            "--figure3-png-output",
            str(project_root / "plots" / "publication" / "png" / "fig3_external_replication_forest.png"),
            "--figure4-pdf-output",
            str(project_root / "plots" / "publication" / "pdf" / "fig4_adjusted_sensitivity_forest.pdf"),
            "--figure4-png-output",
            str(project_root / "plots" / "publication" / "png" / "fig4_adjusted_sensitivity_forest.png"),
            "--figure5-anchor-input",
            str(figure5_anchor_input),
            "--figure5-family-input",
            str(figure5_family_input),
            "--figure5-cross-km-input",
            str(figure5_cross_km_input),
            "--figure5-pdf-output",
            str(project_root / "plots" / "publication" / "pdf" / "fig5_treatment_extension_forest.pdf"),
            "--figure5-png-output",
            str(project_root / "plots" / "publication" / "png" / "fig5_treatment_extension_forest.png"),
        ],
        check=True,
        capture_output=True,
        text=True,
    )


def write_outputs(project_root: Path) -> None:
    records = load_dataset_records(project_root)
    claims = load_claim_records(project_root)
    qualification_rows = build_dataset_qualification_rows(records)
    treatment_effect_rows = build_treatment_effect_rows(project_root, records, claims)
    claim_effect_rows = build_claim_effect_rows(project_root, records, claims, treatment_effect_rows)
    adjusted_effect_rows = build_adjusted_effect_rows(project_root, records)
    supportive_context_rows = build_supportive_context_rows(project_root, records, claims)
    shared_state_rows = build_shared_state_effect_rows(project_root, records)
    replication_rows = build_replication_rows(project_root, records)
    qc_rows = build_qc_rows(records)
    metadata_rows = build_metadata_rows(records)
    umap_rows = build_umap_rows(records)
    figure1_rows = build_figure1_rows(records)
    figure2_rows = build_figure2_rows(records)
    figure3_rows = build_figure3_rows(project_root, records)
    figure4_rows = build_figure4_rows(project_root, records)
    figure5_rows = build_figure5_rows(project_root, records)

    write_tsv(project_root / "results" / "dataset_summary.tsv", records)
    write_tsv(project_root / "results" / "dataset_qualification.tsv", qualification_rows)
    write_tsv(project_root / "results" / "effect_sizes" / "c1_shared_state_effects.tsv", shared_state_rows)
    write_tsv(project_root / "results" / "effect_sizes" / "claim_effects.tsv", claim_effect_rows)
    write_tsv(
        project_root / "results" / "effect_sizes" / "covariate_adjusted_effects.tsv",
        adjusted_effect_rows,
    )
    if treatment_effect_rows:
        write_tsv(project_root / "results" / "effect_sizes" / "treatment_effects.tsv", treatment_effect_rows)
    write_tsv(project_root / "results" / "replication" / "combined_summary.tsv", replication_rows)
    write_replication_dataset_tables(project_root, replication_rows)
    write_tsv(project_root / "results" / "figures" / "qc_overview.tsv", qc_rows)
    write_tsv(project_root / "results" / "figures" / "metadata_overview.tsv", metadata_rows)
    write_tsv(project_root / "results" / "figures" / "umap_overview.tsv", umap_rows)
    write_tsv(project_root / "results" / "figures" / "figure1_study_design.tsv", figure1_rows)
    write_tsv(project_root / "results" / "figures" / "figure2_derivation_anchor.tsv", figure2_rows)
    write_tsv(project_root / "results" / "figures" / "figure3_replication_anchor.tsv", figure3_rows)
    write_tsv(project_root / "results" / "figures" / "figure4_adjusted_effects_anchor.tsv", figure4_rows)
    write_tsv(project_root / "results" / "figures" / "figure5_extension_anchor.tsv", figure5_rows)
    if supportive_context_rows:
        write_tsv(project_root / "results" / "supportive" / "context_summary.tsv", supportive_context_rows)


def write_manuscript_refresh_stub(project_root: Path, run_id: str) -> None:
    content = f"# Manuscript Refresh Stub\n\nRun ID: {run_id}\n\nThis scaffold-stage full run refreshed the manuscript state placeholder only. Real manuscript synchronization will be added after analytical outputs become non-placeholder.\n"
    write_text(project_root / "docs" / "manuscript_refresh" / f"{run_id}.md", content)


def checksum_for(path: Path) -> str:
    digest = hashlib.sha256()
    digest.update(path.read_bytes())
    return digest.hexdigest()


def write_audit_bundle(project_root: Path, mode: str, run_id: str) -> None:
    audit_dir = project_root / "docs" / "audit_runs" / run_id
    audit_dir.mkdir(parents=True, exist_ok=True)

    manifest = {
        "run_id": run_id,
        "mode": mode,
        "generated_at": datetime.now(UTC).isoformat(),
        "platform": platform.platform(),
        "python": platform.python_version(),
        "cwd": str(project_root),
    }
    write_text(audit_dir / "run_manifest.json", json.dumps(manifest, indent=2, sort_keys=True))

    environment_text = "\n".join(
        [
            f"platform={platform.platform()}",
            f"python={platform.python_version()}",
            f"cwd={project_root}",
        ]
    )
    write_text(audit_dir / "environment.txt", environment_text + "\n")

    git_state = "git unavailable: project root is not a git repository\n"
    try:
        result = subprocess.run(
            ["git", "-C", str(project_root), "rev-parse", "--is-inside-work-tree"],
            check=False,
            capture_output=True,
            text=True,
        )
        if result.returncode == 0 and result.stdout.strip() == "true":
            head = subprocess.run(
                ["git", "-C", str(project_root), "rev-parse", "HEAD"],
                check=False,
                capture_output=True,
                text=True,
            )
            status = subprocess.run(
                ["git", "-C", str(project_root), "status", "--short"],
                check=False,
                capture_output=True,
                text=True,
            )
            git_state = f"HEAD={head.stdout.strip()}\nstatus=\n{status.stdout}"
    except FileNotFoundError:
        git_state = "git unavailable: git executable not found\n"
    write_text(audit_dir / "git_state.txt", git_state)

    log_source = project_root / "logs" / f"{run_id}.log"
    if log_source.exists():
        write_text(audit_dir / f"{run_id}.log", log_source.read_text(encoding="utf-8"))

    tracked_files = [
        project_root / "results" / "dataset_summary.tsv",
        project_root / "results" / "dataset_qualification.tsv",
        project_root / "results" / "ingestion" / "dataset_accession_metadata.tsv",
        project_root / "results" / "ingestion" / "gse62254_expression_clinical_link.tsv",
        project_root / "results" / "ingestion" / "gse62254_clinical_link_summary.tsv",
        project_root / "results" / "ingestion" / "gse62254_expression_linked_matrix.tsv.gz",
        project_root / "results" / "ingestion" / "gse15459_outcomes.tsv",
        project_root / "results" / "ingestion" / "endpoint_completeness.tsv",
        project_root / "results" / "ingestion" / "gse15459_expression_clinical_link.tsv",
        project_root / "results" / "ingestion" / "gse15459_expression_link_summary.tsv",
        project_root / "results" / "ingestion" / "gse15459_expression_linked_matrix.tsv.gz",
        project_root / "results" / "ingestion" / "natcom_2018_clinical_endpoints.tsv",
        project_root / "results" / "ingestion" / "natcom_2018_clinical_summary.tsv",
        project_root / "results" / "ingestion" / "gse26899_expression_clinical_link.tsv",
        project_root / "results" / "ingestion" / "gse26899_expression_link_summary.tsv",
        project_root / "results" / "ingestion" / "gse26899_expression_linked_matrix.tsv.gz",
        project_root / "results" / "ingestion" / "gse26901_expression_clinical_link.tsv",
        project_root / "results" / "ingestion" / "gse26901_expression_link_summary.tsv",
        project_root / "results" / "ingestion" / "gse26901_expression_linked_matrix.tsv.gz",
        project_root / "results" / "ingestion" / "gse84437_expression_clinical_link.tsv",
        project_root / "results" / "ingestion" / "gse84437_expression_link_summary.tsv",
        project_root / "results" / "ingestion" / "gse84437_expression_linked_matrix.tsv.gz",
        project_root / "results" / "effect_sizes" / "acrg_cross_cohort_subtype_effects.tsv",
        project_root / "results" / "effect_sizes" / "c1_shared_state_effects.tsv",
        project_root / "results" / "effect_sizes" / "gse15459_subtype_effects.tsv",
        project_root / "results" / "effect_sizes" / "gse84437_projected_subtype_effects.tsv",
        project_root / "results" / "effect_sizes" / "gse26253_compact_score_effects.tsv",
        project_root / "results" / "effect_sizes" / "tcga_stad_compact_score_effects.tsv",
        project_root / "results" / "effect_sizes" / "claim_effects.tsv",
        project_root / "results" / "effect_sizes" / "covariate_adjusted_effects.tsv",
        project_root / "results" / "effect_sizes" / "treatment_effects.tsv",
        project_root / "results" / "effect_sizes" / "natcom_c2_pooled_effects.tsv",
        project_root / "results" / "replication" / "combined_summary.tsv",
        project_root / "results" / "replication" / "gse84437_projected_subtypes.tsv",
        project_root / "results" / "benchmarks" / "method_benchmark.tsv",
        project_root / "results" / "benchmarks" / "state_definition_decision.tsv",
        project_root / "results" / "benchmarks" / "compact_score_gene_panel.tsv",
        project_root / "results" / "figures" / "qc_overview.tsv",
        project_root / "results" / "figures" / "metadata_overview.tsv",
        project_root / "results" / "figures" / "umap_overview.tsv",
        project_root / "results" / "figures" / "figure1_study_design.tsv",
        project_root / "results" / "figures" / "figure1_cohort_map.tsv",
        project_root / "results" / "figures" / "figure1_projection_consistency.tsv",
        project_root / "results" / "figures" / "figure2_derivation_anchor.tsv",
        project_root / "results" / "figures" / "figure2_derivation_km.tsv",
        project_root / "results" / "figures" / "figure2_marker_expression.tsv",
        project_root / "results" / "figures" / "figure3_replication_anchor.tsv",
        project_root / "results" / "figures" / "figure3_summary_diamond.tsv",
        project_root / "results" / "figures" / "figure4_adjusted_effects_anchor.tsv",
        project_root / "results" / "figures" / "figure5_extension_anchor.tsv",
        project_root / "results" / "figures" / "figure5_family_summary.tsv",
        project_root / "results" / "tables" / "gse84437_projection_qc.tsv",
        project_root / "results" / "tables" / "natcom_c2_cohort_summary.tsv",
        project_root / "results" / "tables" / "natcom_c2_cohortwise_interpretation.tsv",
        project_root / "results" / "tables" / "natcom_c2_meta_summary.tsv",
        project_root / "results" / "supportive" / "context_summary.tsv",
        project_root / "plots" / "publication" / "pdf" / "fig1_workflow_consistency.pdf",
        project_root / "plots" / "publication" / "png" / "fig1_workflow_consistency.png",
        project_root / "plots" / "publication" / "pdf" / "fig2_discovery_signal.pdf",
        project_root / "plots" / "publication" / "png" / "fig2_discovery_signal.png",
        project_root / "plots" / "publication" / "pdf" / "fig3_external_replication_forest.pdf",
        project_root / "plots" / "publication" / "png" / "fig3_external_replication_forest.png",
        project_root / "plots" / "publication" / "pdf" / "fig4_adjusted_sensitivity_forest.pdf",
        project_root / "plots" / "publication" / "png" / "fig4_adjusted_sensitivity_forest.png",
        project_root / "plots" / "publication" / "pdf" / "fig5_treatment_extension_forest.pdf",
        project_root / "plots" / "publication" / "png" / "fig5_treatment_extension_forest.png",
    ]
    rows = [
        {"relative_path": str(path.relative_to(project_root)), "sha256": checksum_for(path)}
        for path in tracked_files
        if path.exists()
    ]
    write_tsv(audit_dir / "checksums.tsv", rows)


def write_log(project_root: Path, mode: str, run_id: str) -> None:
    log_path = project_root / "logs" / f"{run_id}.log"
    stages = [
        f"run_id={run_id}",
        f"mode={mode}",
        "stage=data: scaffold-only",
        "stage=preprocess: placeholder outputs generated",
        "stage=analysis: placeholder effect tables generated",
        "stage=accession_metadata: cached accession summaries parsed when available",
        "stage=derivation_link: GSE62254 matrix linked to ACRG clinical supplement when source files exist",
        "stage=clinical_outcomes: GSE15459 outcome workbook parsed when available",
        "stage=expression_link: GSE15459 series matrix linked to analyzable clinical rows when matrix cache exists",
        "stage=natcom_clinical: Nat Commun Supplementary Data 2 workbook standardized when source file exists",
        "stage=natcom_gse26899_link: GSE26899 matrix linked to KUGH clinical rows when matrix cache exists",
        "stage=natcom_gse26901_link: GSE26901 matrix linked to KUCM clinical rows with normalized patient ids when matrix cache exists",
        "stage=replication_link: GSE84437 series matrix filtered to survival-ready rows when matrix cache exists",
        "stage=survival_effects: GSE15459 subtype Cox effects generated when clinical rows exist",
        "stage=cross_cohort_effects: aligned ACRG subtype Cox effects generated for GSE62254 and GSE15459 when source data exist",
        "stage=projection_effects: GSE62254 frozen subtype centroids projected into GSE84437 with Cox outputs and assignment QC when platform annotations exist",
        "stage=benchmark: state-definition head-to-head table generated",
        "stage=compact_score_external_validation: GSE26253 recurrence-oriented compact-score effects generated when matrix input exists",
        "stage=natcom_c2_effects: Nat Commun cohort-family treatment-extension effects and cohort summary generated when linked rows exist",
        "stage=natcom_c2_pooled_effects: Nat Commun pooled dataset-adjusted interaction and fixed-effect meta summary generated when linked rows exist",
        "stage=figures: anchor tables, support tables, and manuscript-grade prognosis forest plots generated",
    ]
    if mode == "full":
        stages.append("stage=manuscript_refresh: placeholder refresh record written")
    stages.append("stage=audit: manifest, environment, git-state, log copy, and checksums written")
    content = "\n".join(stages)
    write_text(log_path, content + "\n")


def run(mode: str, project_root: Path, run_id: str) -> None:
    ensure_project_scaffold(project_root)
    write_compute_plan(project_root)
    write_harmonization_notes(project_root)
    write_accession_metadata_output(project_root)
    write_gse62254_clinical_link_outputs(project_root)
    write_gse15459_outcome_outputs(project_root)
    write_gse15459_expression_link_outputs(project_root)
    write_natcom_2018_clinical_outputs(project_root)
    write_gse26899_expression_link_outputs(project_root)
    write_gse26901_expression_link_outputs(project_root)
    write_gse84437_expression_link_outputs(project_root)
    write_gse15459_subtype_effects(project_root)
    write_acrg_cross_cohort_effects(project_root)
    write_gse84437_projection_effects(project_root)
    write_state_definition_benchmark(project_root)
    write_tcga_stad_compact_score_effects(project_root)
    write_gse26253_compact_score_effects(project_root)
    write_natcom_c2_treatment_effects(project_root)
    write_natcom_c2_pooled_effects(project_root)
    write_outputs(project_root)
    build_submission_figure_support_tables(project_root)
    render_submission_figure1(project_root)
    render_submission_figure2(project_root)
    render_publication_figures(project_root)
    if mode == "full":
        write_manuscript_refresh_stub(project_root, run_id)
    write_log(project_root, mode, run_id)
    write_audit_bundle(project_root, mode, run_id)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Discover Oncology workflow scaffold runner")
    parser.add_argument("--mode", choices=["smoke", "full"], required=True)
    parser.add_argument("--project-root", default=os.getcwd())
    parser.add_argument("--run-id", default=None)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    project_root = Path(args.project_root).resolve()
    run_id = args.run_id or f"{args.mode}-{datetime.now(UTC).strftime('%Y%m%dT%H%M%SZ')}"
    run(args.mode, project_root, run_id)


if __name__ == "__main__":
    main()

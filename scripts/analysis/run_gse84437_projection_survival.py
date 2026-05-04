from __future__ import annotations

import argparse
import csv
import gzip
import math
import re
import urllib.request
from pathlib import Path

import pandas as pd
from lifelines import CoxPHFitter


REFERENCE_SUBTYPE = "MSI"
COMPARISON_SUBTYPES = ["EMT", "MSS/TP53-", "MSS/TP53+"]
SUBTYPE_ORDER = [REFERENCE_SUBTYPE, *COMPARISON_SUBTYPES]
SCORE_COLUMNS = {
    "EMT": "score_emt",
    "MSI": "score_msi",
    "MSS/TP53-": "score_mss_tp53_minus",
    "MSS/TP53+": "score_mss_tp53_plus",
}


def format_number(value: float) -> str:
    return format(float(value), ".6g")


def parse_stage_token(raw_value: object) -> int | None:
    value = str(raw_value).strip().upper()
    if not value:
        return None
    match = re.search(r"(\d+)", value)
    if match is None:
        return None
    return int(match.group(1))


def load_probe_gene_map(path: Path) -> dict[str, str]:
    frame = pd.read_csv(path, sep="\t", dtype=str)
    frame = frame.dropna(subset=["probe_id", "gene_symbol"]).copy()
    frame["probe_id"] = frame["probe_id"].str.strip()
    frame["gene_symbol"] = frame["gene_symbol"].str.strip()
    frame = frame[(frame["probe_id"] != "") & (frame["gene_symbol"] != "")].copy()
    return dict(zip(frame["probe_id"], frame["gene_symbol"]))


def parse_gpl570_text_file(path: Path, probe_ids: set[str]) -> dict[str, str]:
    mapping: dict[str, str] = {}
    marker = "!platform_table_begin"
    with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        in_table = False
        header: list[str] | None = None
        id_idx = None
        symbol_idx = None
        for raw in handle:
            line = raw.rstrip("\n").rstrip("\r")
            if line == marker:
                in_table = True
                continue
            if not in_table:
                continue
            if header is None:
                header = line.split("\t")
                id_idx = header.index("ID")
                symbol_idx = header.index("Gene Symbol")
                continue
            if not line or line.startswith("!"):
                continue
            fields = line.split("\t")
            if len(fields) <= max(id_idx, symbol_idx):
                continue
            probe_id = fields[id_idx].strip()
            if probe_id not in probe_ids:
                continue
            symbol = fields[symbol_idx].strip()
            if symbol and symbol != "---":
                mapping[probe_id] = symbol.split(" /// ")[0].strip()
    return mapping


def parse_gpl6947_bgx_file(path: Path, probe_ids: set[str]) -> dict[str, str]:
    mapping: dict[str, str] = {}
    with gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="") as handle:
        in_probes = False
        header: list[str] | None = None
        probe_idx = None
        symbol_idx = None
        for raw in handle:
            line = raw.rstrip("\n").rstrip("\r")
            if line == "[Probes]":
                in_probes = True
                continue
            if not in_probes:
                continue
            if header is None:
                header = line.split("\t")
                probe_idx = header.index("Probe_Id")
                symbol_idx = header.index("Symbol")
                continue
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) <= max(probe_idx, symbol_idx):
                continue
            probe_id = fields[probe_idx].strip()
            if probe_id not in probe_ids:
                continue
            symbol = fields[symbol_idx].strip()
            if symbol and symbol != "---" and symbol not in {"negative", "biotin"}:
                mapping[probe_id] = symbol.split(" /// ")[0].strip()
    return mapping


def fetch_probe_gene_map(platform_accession: str, symbol_column: str, probe_ids: set[str]) -> dict[str, str]:
    url = (
        f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={platform_accession}"
        "&targ=self&form=text&view=data"
    )
    mapping: dict[str, str] = {}
    with urllib.request.urlopen(url, timeout=180) as response:
        in_table = False
        header: list[str] | None = None
        id_idx = None
        symbol_idx = None
        for raw in response:
            line = raw.decode("utf-8", "replace").strip()
            if line == "!platform_table_begin":
                in_table = True
                continue
            if not in_table:
                continue
            if header is None:
                header = line.split("\t")
                id_idx = header.index("ID")
                symbol_idx = header.index(symbol_column)
                continue
            if not line or line.startswith("!") or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) <= max(id_idx, symbol_idx):
                continue
            probe_id = fields[id_idx].strip()
            if probe_id not in probe_ids:
                continue
            symbol = fields[symbol_idx].strip()
            if not symbol or symbol == "---":
                continue
            mapping[probe_id] = symbol.split(" /// ")[0].strip()
    return mapping


def load_or_build_probe_gene_map(
    path: Path,
    platform_accession: str,
    symbol_column: str,
    probe_ids: set[str],
) -> dict[str, str]:
    if path.exists():
        mapping = load_probe_gene_map(path)
        if mapping:
            return mapping

    raw_dir = path.parent / "raw"
    if platform_accession == "GPL570":
        raw_path = raw_dir / "GPL570_data.txt"
        if raw_path.exists():
            mapping = parse_gpl570_text_file(raw_path, probe_ids)
            if mapping:
                path.parent.mkdir(parents=True, exist_ok=True)
                with path.open("w", encoding="utf-8", newline="") as handle:
                    writer = csv.DictWriter(handle, fieldnames=["probe_id", "gene_symbol"], delimiter="\t")
                    writer.writeheader()
                    for probe_id in sorted(mapping):
                        writer.writerow({"probe_id": probe_id, "gene_symbol": mapping[probe_id]})
                return mapping

    if platform_accession == "GPL6947":
        raw_path = raw_dir / "GPL6947.bgx.gz"
        if raw_path.exists():
            mapping = parse_gpl6947_bgx_file(raw_path, probe_ids)
            if mapping:
                path.parent.mkdir(parents=True, exist_ok=True)
                with path.open("w", encoding="utf-8", newline="") as handle:
                    writer = csv.DictWriter(handle, fieldnames=["probe_id", "gene_symbol"], delimiter="\t")
                    writer.writeheader()
                    for probe_id in sorted(mapping):
                        writer.writerow({"probe_id": probe_id, "gene_symbol": mapping[probe_id]})
                return mapping

    mapping = fetch_probe_gene_map(platform_accession, symbol_column, probe_ids)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["probe_id", "gene_symbol"], delimiter="\t")
        writer.writeheader()
        for probe_id in sorted(mapping):
            writer.writerow({"probe_id": probe_id, "gene_symbol": mapping[probe_id]})
    return mapping


def load_expression_matrix(path: Path) -> pd.DataFrame:
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        frame = pd.read_csv(handle, sep="\t", index_col=0)
    frame.index = frame.index.astype(str)
    frame.columns = frame.columns.astype(str)
    return frame.apply(pd.to_numeric)


def collapse_to_gene_matrix(matrix: pd.DataFrame, probe_gene_map: dict[str, str]) -> pd.DataFrame:
    probe_index = pd.Index([probe for probe in matrix.index if probe in probe_gene_map])
    collapsed = matrix.loc[probe_index].copy()
    collapsed.insert(0, "gene_symbol", [probe_gene_map[probe] for probe in probe_index])
    collapsed = collapsed.groupby("gene_symbol", sort=True).mean(numeric_only=True)
    collapsed.index = collapsed.index.astype(str)
    return collapsed


def zscore_by_gene(matrix: pd.DataFrame) -> pd.DataFrame:
    means = matrix.mean(axis=1)
    stdevs = matrix.std(axis=1, ddof=0).replace(0, float("nan"))
    standardized = matrix.sub(means, axis=0).div(stdevs, axis=0)
    standardized = standardized.dropna(axis=0, how="any")
    return standardized


def zscore_by_vector(matrix: pd.DataFrame) -> pd.DataFrame:
    means = matrix.mean(axis=0)
    stdevs = matrix.std(axis=0, ddof=0).replace(0, float("nan"))
    standardized = matrix.sub(means, axis=1).div(stdevs, axis=1)
    standardized = standardized.dropna(axis=1, how="any")
    return standardized


def build_centroids(training_matrix: pd.DataFrame, clinical_frame: pd.DataFrame) -> pd.DataFrame:
    subtype_by_sample = clinical_frame.set_index("gsm_id")["molecular_subtype_label"]
    subtype_by_sample = subtype_by_sample.loc[
        subtype_by_sample.index.isin(training_matrix.columns) & subtype_by_sample.isin(SUBTYPE_ORDER)
    ]
    aligned_matrix = training_matrix.loc[:, subtype_by_sample.index]
    standardized = zscore_by_gene(aligned_matrix)

    centroids = {}
    for subtype in SUBTYPE_ORDER:
        sample_ids = subtype_by_sample[subtype_by_sample == subtype].index
        centroids[subtype] = standardized.loc[:, sample_ids].mean(axis=1)
    return pd.DataFrame(centroids)


def project_subtypes(replication_matrix: pd.DataFrame, centroids: pd.DataFrame) -> tuple[pd.DataFrame, int]:
    shared_genes = centroids.index.intersection(replication_matrix.index)
    if len(shared_genes) == 0:
        raise ValueError("No shared genes found between derivation centroids and replication matrix")

    replication_shared = zscore_by_gene(replication_matrix.loc[shared_genes])
    shared_genes = centroids.index.intersection(replication_shared.index)
    centroid_shared = centroids.loc[shared_genes]

    sample_vectors = zscore_by_vector(replication_shared)
    centroid_vectors = zscore_by_vector(centroid_shared)
    shared_columns = sample_vectors.columns
    score_frame = sample_vectors.loc[:, shared_columns].T.dot(centroid_vectors) / len(shared_genes)
    score_frame["predicted_subtype"] = score_frame.idxmax(axis=1)
    score_frame = score_frame.rename(columns={subtype: SCORE_COLUMNS[subtype] for subtype in SUBTYPE_ORDER})
    score_frame.index.name = "gsm_id"
    return score_frame.reset_index(), len(shared_genes)


def softmax_probabilities(scores: dict[str, float]) -> dict[str, float]:
    max_score = max(scores.values())
    exponentials = {label: math.exp(value - max_score) for label, value in scores.items()}
    denominator = sum(exponentials.values())
    return {label: value / denominator for label, value in exponentials.items()}


def rank_margin(margin: float, low_cutoff: float, high_cutoff: float) -> str:
    if margin <= low_cutoff:
        return "low"
    if margin <= high_cutoff:
        return "medium"
    return "high"


def build_assignment_rows(
    score_frame: pd.DataFrame,
    replication_clinical: pd.DataFrame,
    shared_gene_count: int,
) -> list[dict[str, str]]:
    clinical_columns = replication_clinical[["gsm_id", "sample_title"]].copy()
    merged = clinical_columns.merge(score_frame, on="gsm_id", how="inner")
    score_columns = [
        SCORE_COLUMNS["EMT"],
        SCORE_COLUMNS["MSI"],
        SCORE_COLUMNS["MSS/TP53-"],
        SCORE_COLUMNS["MSS/TP53+"],
    ]
    margin_values = []
    for row in merged.to_dict(orient="records"):
        ranked = sorted(
            [
                ("EMT", float(row[SCORE_COLUMNS["EMT"]])),
                ("MSI", float(row[SCORE_COLUMNS["MSI"]])),
                ("MSS/TP53-", float(row[SCORE_COLUMNS["MSS/TP53-"]])),
                ("MSS/TP53+", float(row[SCORE_COLUMNS["MSS/TP53+"]])),
            ],
            key=lambda item: item[1],
            reverse=True,
        )
        margin_values.append(ranked[0][1] - ranked[1][1])
    low_cutoff = float(pd.Series(margin_values).quantile(1 / 3))
    high_cutoff = float(pd.Series(margin_values).quantile(2 / 3))

    rows: list[dict[str, str]] = []
    for row in merged.to_dict(orient="records"):
        predicted_subtype = row["predicted_subtype"]
        scores = {
            "EMT": float(row[SCORE_COLUMNS["EMT"]]),
            "MSI": float(row[SCORE_COLUMNS["MSI"]]),
            "MSS/TP53-": float(row[SCORE_COLUMNS["MSS/TP53-"]]),
            "MSS/TP53+": float(row[SCORE_COLUMNS["MSS/TP53+"]]),
        }
        ranked_scores = sorted(scores.items(), key=lambda item: item[1], reverse=True)
        runner_up_subtype, runner_up_score = ranked_scores[1]
        assignment_margin = ranked_scores[0][1] - runner_up_score
        probabilities = softmax_probabilities(scores)
        assignment_entropy = -sum(value * math.log(value) for value in probabilities.values()) / math.log(len(probabilities))
        rows.append(
            {
                "dataset_id": "GSE84437",
                "gsm_id": row["gsm_id"],
                "sample_title": row["sample_title"],
                "predicted_subtype": predicted_subtype,
                "assignment_score": format_number(row[SCORE_COLUMNS[predicted_subtype]]),
                "runner_up_subtype": runner_up_subtype,
                "runner_up_score": format_number(runner_up_score),
                "assignment_margin": format_number(assignment_margin),
                "assignment_margin_rank": rank_margin(assignment_margin, low_cutoff, high_cutoff),
                "assignment_entropy": format_number(assignment_entropy),
                "score_emt": format_number(row[SCORE_COLUMNS["EMT"]]),
                "score_msi": format_number(row[SCORE_COLUMNS["MSI"]]),
                "score_mss_tp53_minus": format_number(row[SCORE_COLUMNS["MSS/TP53-"]]),
                "score_mss_tp53_plus": format_number(row[SCORE_COLUMNS["MSS/TP53+"]]),
                "shared_gene_count": str(shared_gene_count),
            }
        )
    return rows


def build_qc_rows(assignment_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    frame = pd.DataFrame(assignment_rows)
    frame["assignment_margin"] = pd.to_numeric(frame["assignment_margin"])
    frame["assignment_entropy"] = pd.to_numeric(frame["assignment_entropy"])
    rows: list[dict[str, str]] = []
    for predicted_subtype, subframe in frame.groupby("predicted_subtype", sort=True):
        rows.append(
            {
                "dataset_id": "GSE84437",
                "predicted_subtype": str(predicted_subtype),
                "sample_count": str(len(subframe)),
                "median_assignment_margin": format_number(subframe["assignment_margin"].median()),
                "median_assignment_entropy": format_number(subframe["assignment_entropy"].median()),
                "low_margin_count": str(int((subframe["assignment_margin_rank"] == "low").sum())),
                "shared_gene_count": str(subframe["shared_gene_count"].iloc[0]),
            }
        )
    return rows


def build_replication_frame(replication_clinical: pd.DataFrame, score_frame: pd.DataFrame) -> pd.DataFrame:
    assignment_columns = ["gsm_id", "predicted_subtype"]
    if "assignment_margin" in score_frame.columns:
        assignment_columns.append("assignment_margin")
    if "assignment_margin_rank" in score_frame.columns:
        assignment_columns.append("assignment_margin_rank")
    merged = replication_clinical.merge(score_frame[assignment_columns], on="gsm_id", how="inner")
    merged["overall_survival_months"] = pd.to_numeric(merged["overall_survival_months"])
    merged["overall_survival_event"] = pd.to_numeric(merged["overall_survival_event"])
    merged["pt_stage_numeric"] = merged["pt_stage"].map(parse_stage_token)
    merged["pn_stage_numeric"] = merged["pn_stage"].map(parse_stage_token)
    if "assignment_margin" in merged.columns:
        merged["assignment_margin"] = pd.to_numeric(merged["assignment_margin"])
    merged["subtype"] = pd.Categorical(merged["predicted_subtype"], categories=SUBTYPE_ORDER, ordered=False)
    return merged


def fit_model(frame: pd.DataFrame, model_name: str) -> list[dict[str, str]]:
    covariate_columns: list[str] = []
    model_frame = frame[
        [
            "overall_survival_months",
            "overall_survival_event",
            "subtype",
            "pt_stage_numeric",
            "pn_stage_numeric",
        ]
    ].copy()
    subtype_dummies = pd.get_dummies(model_frame["subtype"], prefix="subtype")
    reference_column = f"subtype_{REFERENCE_SUBTYPE}"
    if reference_column in subtype_dummies.columns:
        subtype_dummies = subtype_dummies.drop(columns=[reference_column])
    covariate_columns.extend(sorted(subtype_dummies.columns))
    model_frame = pd.concat([model_frame, subtype_dummies], axis=1)

    if model_name == "cox_pt_pn_adjusted":
        model_frame = model_frame.dropna(subset=["pt_stage_numeric", "pn_stage_numeric"]).copy()
        covariate_columns.extend(["pt_stage_numeric", "pn_stage_numeric"])

    cox_frame = model_frame[["overall_survival_months", "overall_survival_event", *covariate_columns]].copy()
    fitter = CoxPHFitter()
    fitter.fit(cox_frame, duration_col="overall_survival_months", event_col="overall_survival_event")

    rows: list[dict[str, str]] = []
    for comparison in COMPARISON_SUBTYPES:
        coefficient_name = f"subtype_{comparison}"
        summary = fitter.summary.loc[coefficient_name]
        rows.append(
            {
                "dataset_id": "GSE84437",
                "outcome": "overall_survival",
                "model": model_name,
                "reference_subtype": REFERENCE_SUBTYPE,
                "comparison_subtype": comparison,
                "effect_type": "hazard_ratio",
                "effect": format_number(summary["exp(coef)"]),
                "ci_lower": format_number(summary["exp(coef) lower 95%"]),
                "ci_upper": format_number(summary["exp(coef) upper 95%"]),
                "pvalue": format_number(summary["p"]),
                "n": str(len(model_frame)),
                "events": str(int(model_frame["overall_survival_event"].sum())),
                "covariates": "pt_stage,pn_stage" if model_name == "cox_pt_pn_adjusted" else "none",
            }
        )
    return rows


def build_low_margin_sensitivity_rows(replication_frame: pd.DataFrame) -> list[dict[str, str]]:
    if "assignment_margin_rank" not in replication_frame.columns:
        raise ValueError("assignment_margin_rank column is required for low-margin sensitivity analysis")

    removed_frame = replication_frame[replication_frame["assignment_margin_rank"] == "low"].copy()
    retained_frame = replication_frame[replication_frame["assignment_margin_rank"] != "low"].copy()

    rows: list[dict[str, str]] = []
    for analysis_set, frame in [
        ("all_projected_samples", replication_frame),
        ("excluding_low_margin_assignments", retained_frame),
    ]:
        excluded_n = 0 if analysis_set == "all_projected_samples" else len(removed_frame)
        excluded_events = 0 if analysis_set == "all_projected_samples" else int(removed_frame["overall_survival_event"].sum())
        for model_name in ["cox_unadjusted", "cox_pt_pn_adjusted"]:
            for row in fit_model(frame, model_name):
                rows.append(
                    {
                        "dataset_id": row["dataset_id"],
                        "analysis_set": analysis_set,
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
                        "events": row["events"],
                        "covariates": row["covariates"],
                        "excluded_margin_group": "low" if analysis_set != "all_projected_samples" else "none",
                        "excluded_n": str(excluded_n),
                        "excluded_events": str(excluded_events),
                    }
                )
    return rows


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Project frozen GSE62254 subtypes into GSE84437 and fit survival models")
    parser.add_argument("--training-clinical-input", required=True)
    parser.add_argument("--training-matrix-input", required=True)
    parser.add_argument("--replication-clinical-input", required=True)
    parser.add_argument("--replication-matrix-input", required=True)
    parser.add_argument("--gpl570-map-input", required=True)
    parser.add_argument("--gpl6947-map-input", required=True)
    parser.add_argument("--assignment-output", required=True)
    parser.add_argument("--effect-output", required=True)
    parser.add_argument("--qc-output", required=True)
    parser.add_argument("--sensitivity-output")
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    training_clinical = pd.read_csv(Path(args.training_clinical_input), sep="\t", dtype=str)
    replication_clinical = pd.read_csv(Path(args.replication_clinical_input), sep="\t", dtype=str)
    training_probe_matrix = load_expression_matrix(Path(args.training_matrix_input))
    replication_probe_matrix = load_expression_matrix(Path(args.replication_matrix_input))

    training_matrix = collapse_to_gene_matrix(
        training_probe_matrix,
        load_or_build_probe_gene_map(
            Path(args.gpl570_map_input),
            "GPL570",
            "Gene Symbol",
            set(training_probe_matrix.index),
        ),
    )
    replication_matrix = collapse_to_gene_matrix(
        replication_probe_matrix,
        load_or_build_probe_gene_map(
            Path(args.gpl6947_map_input),
            "GPL6947",
            "Symbol",
            set(replication_probe_matrix.index),
        ),
    )

    centroids = build_centroids(training_matrix, training_clinical)
    score_frame, shared_gene_count = project_subtypes(replication_matrix, centroids)
    assignment_rows = build_assignment_rows(score_frame, replication_clinical, shared_gene_count)
    assignment_frame = pd.DataFrame(assignment_rows)
    replication_frame = build_replication_frame(replication_clinical, assignment_frame)
    effect_rows = [
        *fit_model(replication_frame, "cox_unadjusted"),
        *fit_model(replication_frame, "cox_pt_pn_adjusted"),
    ]

    write_tsv(Path(args.assignment_output), assignment_rows)
    write_tsv(Path(args.effect_output), effect_rows)
    write_tsv(Path(args.qc_output), build_qc_rows(assignment_rows))
    if args.sensitivity_output:
        write_tsv(Path(args.sensitivity_output), build_low_margin_sensitivity_rows(replication_frame))


if __name__ == "__main__":
    main()
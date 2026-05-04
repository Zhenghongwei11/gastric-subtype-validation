from __future__ import annotations

import argparse
import csv
import gzip
import math
import re
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
    roman_map = {"I": 1, "II": 2, "III": 3, "IV": 4}
    for token in ("IV", "III", "II", "I"):
        if token in value:
            return roman_map[token]
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


def load_or_build_gpl6947_gene_map(path: Path, probe_ids: set[str]) -> dict[str, str]:
    if path.exists():
        mapping = load_probe_gene_map(path)
        if mapping:
            return mapping

    raw_dir = path.parent / "raw"
    raw_path = raw_dir / "GPL6947.bgx.gz"
    if not raw_path.exists():
        raise FileNotFoundError(
            "GPL6947 annotation cache missing. Expected a BGX cache at "
            f"{raw_path}. Run the existing workflow once to populate the cache."
        )
    mapping = parse_gpl6947_bgx_file(raw_path, probe_ids)
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
    if probe_index.empty:
        raise ValueError("No overlapping probes found between matrix and GPL gene-symbol map")
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
    score_frame = sample_vectors.T.dot(centroid_vectors) / len(shared_genes)
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
    dataset_id: str,
    score_frame: pd.DataFrame,
    replication_clinical: pd.DataFrame,
    shared_gene_count: int,
) -> list[dict[str, str]]:
    clinical_columns = replication_clinical[["gsm_id", "matrix_sample_title"]].copy()
    merged = clinical_columns.merge(score_frame, on="gsm_id", how="inner")
    margin_values: list[float] = []
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
        assignment_entropy = -sum(value * math.log(value) for value in probabilities.values()) / math.log(
            len(probabilities)
        )
        rows.append(
            {
                "dataset_id": dataset_id,
                "gsm_id": row["gsm_id"],
                "sample_title": row["matrix_sample_title"],
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


def build_qc_rows(dataset_id: str, assignment_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    frame = pd.DataFrame(assignment_rows)
    frame["assignment_margin"] = pd.to_numeric(frame["assignment_margin"])
    frame["assignment_entropy"] = pd.to_numeric(frame["assignment_entropy"])
    rows: list[dict[str, str]] = []
    for predicted_subtype, subframe in frame.groupby("predicted_subtype", sort=True):
        rows.append(
            {
                "dataset_id": dataset_id,
                "predicted_subtype": str(predicted_subtype),
                "sample_count": str(len(subframe)),
                "median_assignment_margin": format_number(subframe["assignment_margin"].median()),
                "median_assignment_entropy": format_number(subframe["assignment_entropy"].median()),
                "low_margin_count": str(int((subframe["assignment_margin_rank"] == "low").sum())),
                "shared_gene_count": str(subframe["shared_gene_count"].iloc[0]),
            }
        )
    return rows


def build_replication_frame(replication_clinical: pd.DataFrame, assignment_frame: pd.DataFrame) -> pd.DataFrame:
    merged = replication_clinical.merge(
        assignment_frame[["gsm_id", "predicted_subtype", "assignment_margin_rank"]],
        on="gsm_id",
        how="inner",
    )
    merged["overall_survival_months"] = pd.to_numeric(merged["overall_survival_months"])
    merged["overall_survival_event"] = pd.to_numeric(merged["overall_survival_event"])
    merged["stage_numeric"] = merged["ajcc_stage"].map(parse_stage_token)
    merged["subtype"] = pd.Categorical(merged["predicted_subtype"], categories=SUBTYPE_ORDER, ordered=False)
    return merged


def fit_model(dataset_id: str, frame: pd.DataFrame, model_name: str) -> list[dict[str, str]]:
    covariate_columns: list[str] = []
    model_frame = frame[["overall_survival_months", "overall_survival_event", "subtype", "stage_numeric"]].copy()
    subtype_dummies = pd.get_dummies(model_frame["subtype"], prefix="subtype")
    reference_column = f"subtype_{REFERENCE_SUBTYPE}"
    if reference_column in subtype_dummies.columns:
        subtype_dummies = subtype_dummies.drop(columns=[reference_column])
    covariate_columns.extend(sorted(subtype_dummies.columns))
    model_frame = pd.concat([model_frame, subtype_dummies], axis=1)

    covariates = "none"
    if model_name == "cox_stage_adjusted":
        model_frame = model_frame.dropna(subset=["stage_numeric"]).copy()
        model_frame["stage_numeric"] = pd.to_numeric(model_frame["stage_numeric"])
        covariate_columns.append("stage_numeric")
        covariates = "stage"

    cox_frame = model_frame[["overall_survival_months", "overall_survival_event", *covariate_columns]].copy()
    fitter = CoxPHFitter()
    fitter.fit(cox_frame, duration_col="overall_survival_months", event_col="overall_survival_event")

    rows: list[dict[str, str]] = []
    for comparison in COMPARISON_SUBTYPES:
        coefficient_name = f"subtype_{comparison}"
        summary = fitter.summary.loc[coefficient_name]
        rows.append(
            {
                "dataset_id": dataset_id,
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
                "covariates": covariates,
            }
        )
    return rows


def subtype_counts(dataset_id: str, frame: pd.DataFrame, model_name: str) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for subtype, subframe in frame.dropna(subset=["subtype"]).groupby("subtype", observed=True):
        rows.append(
            {
                "dataset_id": dataset_id,
                "outcome": "overall_survival",
                "model": model_name,
                "subtype": str(subtype),
                "n": str(len(subframe)),
                "events": str(int(pd.to_numeric(subframe["overall_survival_event"]).sum())),
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
    parser = argparse.ArgumentParser(
        description="Transfer frozen GSE62254 ACRG centroids into Nat Commun cohort family (GSE26899 and GSE26901)"
    )
    parser.add_argument("--training-clinical-input", required=True)
    parser.add_argument("--training-matrix-input", required=True)
    parser.add_argument("--gpl570-map-input", required=True)
    parser.add_argument("--gpl6947-map-input", required=True)
    parser.add_argument("--gse26899-clinical-input", required=True)
    parser.add_argument("--gse26899-matrix-input", required=True)
    parser.add_argument("--gse26901-clinical-input", required=True)
    parser.add_argument("--gse26901-matrix-input", required=True)
    parser.add_argument("--assignment-output", required=True)
    parser.add_argument("--effect-output", required=True)
    parser.add_argument("--qc-output", required=True)
    parser.add_argument("--subtype-counts-output", required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    training_clinical = pd.read_csv(Path(args.training_clinical_input), sep="\t", dtype=str)
    training_probe_matrix = load_expression_matrix(Path(args.training_matrix_input))
    training_matrix = collapse_to_gene_matrix(
        training_probe_matrix,
        load_probe_gene_map(Path(args.gpl570_map_input)),
    )
    centroids = build_centroids(training_matrix, training_clinical)

    gpl6947_gene_map = load_or_build_gpl6947_gene_map(
        Path(args.gpl6947_map_input),
        probe_ids=set(load_expression_matrix(Path(args.gse26899_matrix_input)).index),
    )

    assignment_rows: list[dict[str, str]] = []
    effect_rows: list[dict[str, str]] = []
    qc_rows: list[dict[str, str]] = []
    counts_rows: list[dict[str, str]] = []

    for dataset_id, clinical_path, matrix_path in [
        ("GSE26899", Path(args.gse26899_clinical_input), Path(args.gse26899_matrix_input)),
        ("GSE26901", Path(args.gse26901_clinical_input), Path(args.gse26901_matrix_input)),
    ]:
        replication_clinical = pd.read_csv(clinical_path, sep="\t", dtype=str)
        replication_probe_matrix = load_expression_matrix(matrix_path)
        replication_matrix = collapse_to_gene_matrix(replication_probe_matrix, gpl6947_gene_map)
        score_frame, shared_gene_count = project_subtypes(replication_matrix, centroids)
        dataset_assignment_rows = build_assignment_rows(dataset_id, score_frame, replication_clinical, shared_gene_count)
        assignment_rows.extend(dataset_assignment_rows)
        assignment_frame = pd.DataFrame(dataset_assignment_rows)
        replication_frame = build_replication_frame(replication_clinical, assignment_frame)

        effect_rows.extend(fit_model(dataset_id, replication_frame, "cox_unadjusted"))
        effect_rows.extend(fit_model(dataset_id, replication_frame, "cox_stage_adjusted"))
        qc_rows.extend(build_qc_rows(dataset_id, dataset_assignment_rows))
        counts_rows.extend(subtype_counts(dataset_id, replication_frame, "cox_unadjusted"))
        counts_rows.extend(subtype_counts(dataset_id, replication_frame.dropna(subset=["stage_numeric"]), "cox_stage_adjusted"))

    write_tsv(Path(args.assignment_output), assignment_rows)
    write_tsv(Path(args.effect_output), effect_rows)
    write_tsv(Path(args.qc_output), qc_rows)
    write_tsv(Path(args.subtype_counts_output), counts_rows)


if __name__ == "__main__":
    main()

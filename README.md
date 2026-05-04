# Gastric subtype validation — reproducibility package

This repository contains a public reproducibility package for a gastric adenocarcinoma transcriptomic subtype validation analysis.

## What is included

- Analysis scripts used to generate the main result tables and publication figures.
- Result tables under `results/` that support the reported effect estimates.
- Publication figures under `plots/publication/` (PDF and PNG).
- A dataset registry under `data/manifest.tsv`.

## What is intentionally NOT included

To avoid inappropriate disclosure of journal submission materials or internal writing workflow artifacts, this repository does **not** include:

- Manuscript text, cover letters, declarations, or any journal upload artifacts.
- Internal writing/review checklists or citation-audit logs.
- Local tokens, credentials, or other secrets.
- Large raw expression matrices; these can be re-obtained from public accessions listed in `data/manifest.tsv`.

## Environment

- Python 3.10+ is recommended.

```bash
python3 -m pip install -r requirements.txt
```

## Reproduce figures and tables

Run the analysis scripts under `scripts/analysis/` and figure scripts under `scripts/figures/`. See the script headers for inputs/outputs.

## Citation

Zenodo: `10.5281/zenodo.20025337` (v1.0.2).

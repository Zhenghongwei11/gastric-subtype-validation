# Gastric subtype validation — reproducibility package

This repository contains a public reproducibility package for a gastric adenocarcinoma transcriptomic subtype validation analysis.

## What is included

- Analysis scripts used to generate the main result tables and publication figures.
- Result tables under `results/` that support the reported effect estimates.
- Publication figures under `plots/publication/` (PDF and PNG).
- A dataset registry under `data/manifest.tsv`.

## Scope

- This package provides analysis code, derived result tables, and publication figures to support reproducibility.
- Primary inputs are public datasets listed in `data/manifest.tsv`.
- Raw expression matrices are not mirrored here; please obtain them from the original public sources.

## Environment

- Python 3.10+ is recommended.

```bash
python3 -m pip install -r requirements.txt
```

## Reproduce figures and tables

Run the analysis scripts under `scripts/analysis/` and figure scripts under `scripts/figures/`. See the script headers for inputs/outputs.

## Citation

Zenodo: `10.5281/zenodo.20018751` (v1.0.1).

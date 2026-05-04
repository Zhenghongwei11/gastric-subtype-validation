# Gastric subtype validation

This repository contains the reproducibility package for a public-data gastric adenocarcinoma study evaluating whether published transcriptomic subtype states retain overall-survival associations across independent cohorts without cohort-specific refitting.

## Study focus

The primary analysis tests a published ACRG subtype definition across one derivation cohort and two external validation cohorts. Overall survival is the primary endpoint. Treatment-by-subtype analyses are included only as exploratory analyses because public treatment metadata are heterogeneous.

## Repository contents

- `scripts/`: analysis, ingestion, figure, and workflow entry points.
- `results/`: source tables, effect-size tables, cohort summaries, and figure-source tables used to support the reported analyses.
- `plots/publication/`: publication figure exports in PNG and PDF format.
- `supplementary_tables/`: compact supplementary tables prepared for journal review.
- `data/manifest.tsv`: public dataset registry and source-accession manifest.
- `logs/`: lightweight run traces from local workflow execution.

Full manuscript text, cover letters, journal-only upload files, OpenSpec governance files, local tokens, raw public expression matrices, and internal AI/review scaffolding are intentionally excluded from this public package.

## Environment

Python 3.10 or newer is recommended.

```bash
python3 -m pip install -r requirements.txt
```

## Reproduce core workflow

Run a lightweight smoke workflow:

```bash
bash scripts/run_smoke.sh
```

Run the fuller local workflow when source matrices and public supplements are available or cached:

```bash
bash scripts/run_full.sh
```

Large raw public matrices are not stored in this GitHub repository. They can be obtained from the public accessions listed in `data/manifest.tsv`.

## Key outputs

- Main effect-size anchors: `results/effect_sizes/`
- External subtype projection anchors: `results/replication/`
- Figure-source tables: `results/figures/`
- Publication figures: `plots/publication/`

## Citation

If using this code or result package, cite the archived GitHub/Zenodo release once the Zenodo DOI is minted from the release record.

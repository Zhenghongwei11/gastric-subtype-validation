# Compute Plan

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

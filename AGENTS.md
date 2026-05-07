# Codex Guidance

This repository analyzes neutral karyotype evolution by fitting chromosome mis-segregation and whole-genome-doubling parameters to observed endpoint karyotypes, then checking fitted behavior with forward simulations.

This branch is a local-running variant. Preserve laptop-friendly execution paths and resource settings, but do not treat reduced local settings as final scientific defaults unless the user explicitly says so.

## Canonical Path

- Treat `updated_workflow.Rmd` as the current analysis document.
- Reusable logic lives under `R/`: `db_utils.R`, `estimate_pmis.R`, `plot_utils.R`, and `workflow_utils.R`.
- Expensive/cached stages live under `scripts/`: `fit_group_markov.R`, `run_forward_validation.R`, and `run_forward_validation_all.R`.
- `ksim2.cpp` is compiled from R through `Rcpp::sourceCpp()`.
- This branch should pull improvements from `main`, but it is not intended to be merged back into `main` wholesale.

## Data Contracts

- Core cached inputs live under `core_data/`.
- Heavy generated outputs live under `results/`.
- Lightweight forward-validation summaries live under `result_summaries/`.
- Karyotypes are expected to use chromosomes 1:22 after marker-chromosome exclusion.
- Copy-number states are modeled as bounded integer states; check `R/estimate_pmis.R` before changing limits.
- `grouped_intervals.Rds` is the main input for grouped fitting and forward validation.
- `group_fit_df.Rds` stores grouped Markov grid-search results.
- `p_mis` / `pmis` is the chromosome mis-segregation parameter; `p_wgd` / `pwgd` is the whole-genome-doubling parameter.

## Commands

```bash
Rscript scripts/fit_group_markov.R --grouped_intervals_path=core_data/grouped_intervals.Rds --output_path=results/group_fit_df.Rds
Rscript scripts/run_forward_validation.R --grouped_intervals_path=core_data/grouped_intervals.Rds --group_fit_path=results/group_fit_df.Rds --forward_group_id=<group_id> --output_path=results/forward_validation_<group_id>.Rds --n_reps=30 --n_null_pairs=100 --n_cores=4
Rscript scripts/run_forward_validation_all.R --grouped_intervals_path=core_data/grouped_intervals.Rds --group_fit_path=results/group_fit_df.Rds --output_dir=results --summary_output_dir=result_summaries --n_reps=30 --n_null_pairs=100 --n_cores=4
```

Forward validation defaults to `--distance_metric=chrom_weighted_wasserstein`; use `--distance_metric=wasserstein` to reproduce the older unweighted endpoint-distance check.

For local laptop smoke checks, reduce settings explicitly, for example `--n_reps=3 --n_null_pairs=10 --n_cores=1 --bottleneck_size=200 --expansion_factor=8`. Run the R Markdown workflow manually in RStudio or another R Markdown environment. There is no formal test, lint, build, or package manifest yet.

## Dependencies

Observed R dependencies include `Rcpp`, `igraph`, `ggplot2`, `parallel`, `transport`, `reshape2`, `DBI`, `RMariaDB`, `tibble`, `stringr`, `data.table`, and `optparse`.

Database-backed workflow steps expect `db_creds.txt` with `HOST`, `DBNAME`, `USER`, and `PASSWORD`; this file is gitignored and must not be committed.

## Change Safety

- Do not silently change model constants, copy-number bounds, grid spacing for `pmis`/`pwgd`, WGD assumptions, bottleneck/expansion settings, null replicate counts, or distance metrics.
- Treat `.Rds` files and generated outputs as scientific artifacts. Ask before regenerating, deleting, moving, or committing new large outputs.
- Reduced local-running settings must be labeled as smoke/debug settings.
- Preserve user changes in dirty files. This repo often contains generated artifacts and exploratory edits.

## Known Fragile Areas

- The local-running branch can drift from `main`; resolve syncs intentionally and document branch-only resource changes.
- Forward validation is still described as provisional model checking, not a final statistical endpoint.
- There are no automated tests for simulator invariants, data schemas, lineage mapping, or posterior-predictive calculations.

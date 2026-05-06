# Codex Guidance

This repository analyzes neutral karyotype evolution by fitting chromosome mis-segregation rates to observed karyotype distributions, then comparing observed endpoints to simulated null populations.

## Canonical Path

- Treat `updated_workflow.Rmd` as the likely current workflow and `workflow.Rmd` as older exploratory context. If this is wrong, ask before changing behavior.
- The main local batch entry point is `R/run_local.R`, which reads `fit_objects.Rds` and writes `data/fits/*.Rds` plus `data/pops/*.Rds`.
- The one-object/HPC entry point is `R/analyse_lineage.R`, called by `fit_jobs.sub`.
- `R/estimate_pmis.R` compiles `ksim2.cpp` with `Rcpp::sourceCpp()` and exposes the helper functions used by both runners.

## Data Contracts

- Karyotypes are expected to be cell-by-22 matrices after excluding the marker chromosome.
- Copy-number states are expected to be integers in `0:8`.
- `fit_objects.Rds` entries contain `id_start`, `id_end`, `replicates`, `K0`, `KT`, `KT_list`, `H0`, `HT`, and `delta_pass`; completed fits add `resdf` and `test_null_res`.
- `K0` is the ancestor/MRCA karyotype matrix, `KT` is the merged endpoint matrix, and `KT_list` keeps per-replicate endpoint matrices.
- `H0` and `HT` are 9-by-22 copy-number histograms.
- `p_misseg`, `pmis`, and `cin_rate` refer to the per-chromosome mis-segregation probability unless a future maintainer clarifies otherwise.

## Commands

```bash
Rscript R/run_local.R
Rscript R/analyse_lineage.R 1 16
sbatch fit_jobs.sub
python3 make_neutrality_plots_schema_flexible.py --csv data/neutral_probability.csv --out_dir neutrality_plots_out --bundle_zip neutrality_plots_bundle.zip
```

Run the R Markdown workflows manually in RStudio or another R Markdown environment. There is no formal test, lint, build, or package manifest yet.

## Dependencies

Observed R dependencies include `Rcpp`, `igraph`, `ggplot2`, `ggrepel`, `parallel`, `transport`, `reshape2`, `DBI`, `RMariaDB`, `tibble`, `stringr`, `pbapply`, and `data.table`.

The plotting script needs Python packages `numpy`, `pandas`, and `matplotlib`.

Database-backed workflow steps expect `db_creds.txt` with `HOST`, `DBNAME`, `USER`, and `PASSWORD`; this file is gitignored and must not be committed.

## Change Safety

- Do not silently change model constants: 22 chromosomes, states `0:8`, `delta_pass * 5`, `dt = 0.1`, `rate = 1.0`, `cull_keep = 1/(2^5)`, population caps, or grid spacing for `pmis`.
- Treat tracked `.Rds` files and `data/` outputs as scientific artifacts. Ask before regenerating, deleting, moving, or committing new large outputs.
- Keep root `db_utils.R` and `R/db_utils.R` in mind; they are currently duplicate-looking files and may drift.
- Preserve user changes in dirty files. This repo often contains generated artifacts and exploratory edits.

## Known Fragile Areas

- `R/estimate_pmis.R` has a likely `get_pop(nsteps)` / `n_steps` naming bug; fix intentionally with a small smoke check, not as drive-by cleanup.
- `updated_workflow.Rmd` has a hard-coded local `root.dir`.
- `R/db_utils.R` builds SQL with string interpolation.
- `fit_jobs.sub` may have an array size that does not match the current `fit_objects.Rds`.
- There are no automated tests for simulator invariants, data schemas, lineage mapping, or p-value calculations.


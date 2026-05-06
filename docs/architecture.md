# Architecture Notes

## Purpose

NeutralKaryotypes estimates whether observed karyotype changes across passaged cell-line samples are consistent with a neutral chromosome mis-segregation model. The pipeline estimates a per-chromosome mis-segregation probability, generates simulated endpoint populations under that fitted rate, and compares simulated and observed karyotype distributions with Wasserstein distances.

## Canonical Workflow

The current workflow appears to be `updated_workflow.Rmd`; `workflow.Rmd` appears to be older exploratory analysis. This should be treated as a question if the distinction matters for a change.

The high-level flow is:

1. Read karyotyped sample IDs from `karyotyped_samples.txt`.
2. Query CLONEID-style MariaDB tables through `R/db_utils.R`, using `db_creds.txt`.
3. Map each endpoint sample to a proxy MRCA/ancestor by passaging ancestry.
4. Fetch karyotype vectors, drop the marker chromosome, and keep integer-resolved 22-chromosome karyotypes.
5. Merge technical or biological replicates into condition-level fit objects.
6. Save `fit_objects.Rds`.
7. Run `R/run_local.R` locally or `R/analyse_lineage.R` through `fit_jobs.sub` on SLURM.
8. Save completed fit summaries in `data/fits/` and simulated null populations in `data/pops/`.
9. Collect per-replicate neutrality p-values into `data/neutral_probability.csv`.
10. Optionally generate plot bundles with `make_neutrality_plots_schema_flexible.py`.

## Core Data Objects

- `karyotypes.Rds`: List with observed karyotype vectors and passaging metadata used by the older workflow.
- `karyotypes_mod.Rds`: Modified karyotype/passaging object from exploratory lineage work.
- `fit_objects.Rds`: Main batch input for simulations. Each entry usually contains:
  - `id_start`: ancestor/MRCA sample ID.
  - `id_end`: merged endpoint or condition ID.
  - `replicates`: endpoint IDs grouped into this condition.
  - `K0`: ancestor cell-by-22 karyotype matrix.
  - `KT`: merged endpoint cell-by-22 karyotype matrix.
  - `KT_list`: per-replicate endpoint matrices.
  - `H0`, `HT`: copy-number histograms with rows for states `0:8` and columns for 22 chromosomes.
  - `delta_pass`: passage difference used to derive simulated time.
- `data/fits/*.Rds`: Completed fit object with likelihood grid `resdf` and `test_null_res`.
- `data/pops/*.Rds`: Lists of simulated endpoint population histograms, keyed by dot-separated karyotype strings.
- `data/neutral_probability.csv`: Flat results table with condition, replicate ID, ancestor, p-value, fitted CIN rate, and passage delta.

## Main Modules

- `ksim2.cpp`: Current Rcpp simulator. It stores a fixed-size population array, runs stochastic cell divisions, applies balanced chromosome mis-segregation, removes invalid daughters, culls at population caps, and records named karyotype-count histograms.
- `R/estimate_pmis.R`: R helper layer for compiling `ksim2.cpp`, converting between matrices and histograms, generating synthetic starting populations, computing negative log likelihoods, and returning simulated endpoint populations.
- `R/run_local.R`: Processes all entries in `fit_objects.Rds`. It estimates `pmis`, generates null populations, computes Wasserstein tests for merged and individual endpoints, and writes outputs.
- `R/analyse_lineage.R`: Similar one-index runner intended for SLURM array execution.
- `R/db_utils.R`: Database credential loading, karyotype extraction, and lineage helper functions.
- `make_neutrality_plots_schema_flexible.py`: Report/plot generator from `data/neutral_probability.csv`.

## Expected Commands

```bash
# Run all pending fit objects locally.
Rscript R/run_local.R

# Run a single fit object by index with a chosen core count.
Rscript R/analyse_lineage.R 1 16

# Submit the configured HPC array.
sbatch fit_jobs.sub

# Generate plots from collected p-values, once Python dependencies exist.
python3 make_neutrality_plots_schema_flexible.py --csv data/neutral_probability.csv --out_dir neutrality_plots_out --bundle_zip neutrality_plots_bundle.zip
```

`Rcpp::sourceCpp("ksim2.cpp")` is called from `R/estimate_pmis.R`, so simulator compilation happens as a side effect of sourcing that file.

## Known Fragile Areas

- `R/estimate_pmis.R` defines `get_pop(p_mis, K0, nsteps, ...)` but passes `n_steps` into `run_karyotype_neutral()`. This looks like a naming bug because callers pass `nsteps = n_steps` and rely on a global `n_steps`.
- `updated_workflow.Rmd` hard-codes a local absolute `root.dir`.
- `db_utils.R` exists both at repo root and under `R/`; the canonical copy is not documented.
- SQL query construction in `R/db_utils.R` interpolates IDs into a string.
- `fit_jobs.sub` declares a fixed SLURM array range that may not match the number of fit objects.
- There is no automated test suite or schema validation for RDS objects.
- The current Python plotting script depends on `pandas` and `matplotlib`; these may not be present in a bare environment.

## Assumptions Not To Change Silently

- Karyotypes are represented as 22 chromosomes after marker-chromosome exclusion.
- Copy-number states are modeled as integers `0:8`.
- Simulation time is derived as `delta_pass * 5` doublings, then converted through `log(2)` growth time with `dt = 0.1`.
- The simulator uses `rate = 1.0` in the R helper functions.
- Population caps, `cull_keep = 1/(2^5)`, number of null replicates, likelihood grid spacing, and `record_every` affect scientific outputs.
- `p_misseg`, `pmis`, and `cin_rate` are treated as the same fitted per-chromosome mis-segregation probability in current outputs.

## Open Questions

- Is `updated_workflow.Rmd` officially canonical, or should `workflow.Rmd` remain the primary documented workflow?
- Should `ksim.cpp` be kept as a historical simulator, a comparison implementation, or removed in a later cleanup?
- Should generated `data/` outputs be reproducible artifacts outside git, or committed result snapshots?
- Should database credentials eventually move from `db_creds.txt` to environment variables or a standard config file?
- What null replicate count is expected for final results: local `100`, SLURM `150`, or another value?


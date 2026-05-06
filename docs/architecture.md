# Architecture Notes

## Purpose

NeutralKaryotypes fits and checks a neutral model of karyotype evolution. The current pipeline estimates grouped chromosome mis-segregation (`p_mis`/`pmis`) and whole-genome-doubling (`p_wgd`/`pwgd`) parameters from observed endpoint karyotypes, then compares fitted behavior to forward simulations using posterior-predictive Wasserstein checks.

This branch is a local-running variant. Its purpose is to keep the workflow usable on a laptop with limited cores and memory while still incorporating structural improvements from `main`.

## Canonical Workflow

`updated_workflow.Rmd` is the current analysis document. It is intended to stay readable and mostly orchestration-focused. Reusable logic lives under `R/`, and expensive stages are run through scripts under `scripts/` with cached outputs.

The high-level flow is:

1. Load cached karyotyped sample IDs from `core_data/karyotyped_samples.txt`.
2. Load or build collapsed passaging data and media annotations in `core_data/`.
3. Build connected passage trees for karyotyped samples.
4. Assign coarse media-derived conditions to passages.
5. Build simulation intervals and grouped intervals.
6. Load or run grouped Markov grid fitting.
7. Inspect ploidy and marginal copy-number diagnostics.
8. Run or load forward-validation results.
9. Summarize cached posterior-predictive checks from `result_summaries/`.

## Core Data Objects

- `core_data/karyotyped_samples.txt`: sample IDs used as the observed karyotype set.
- `core_data/karyotypes.Rds`: cached endpoint karyotype vectors.
- `core_data/db_col.Rds`: collapsed passage-level ancestry graph.
- `core_data/media_raw.Rds`: cached media table for condition annotation.
- `core_data/simulation_intervals.Rds`: ancestor-descendant intervals derived from connected trees and karyotypes.
- `core_data/grouped_intervals.Rds`: grouped interval object used by Markov fitting and forward validation.
- `results/group_fit_df.Rds`: grouped grid-search output from `scripts/fit_group_markov.R`.
- `results/forward_validation_<group>.Rds`: full forward-validation result for one group.
- `result_summaries/forward_validation_<group>_summary.Rds`: compact sidecar summary for reporting.

Older branch artifacts such as root-level `fit_objects.Rds`, `karyotypes.Rds`, `karyotypes_mod.Rds`, `data/fits/`, and `data/pops/` are not the current merged workflow path unless a maintainer explicitly restores them for a branch-specific reason.

## Main Modules

- `updated_workflow.Rmd`: current end-to-end analysis notebook.
- `R/db_utils.R`: database loading, passaging collapse, and karyotype extraction helpers.
- `R/estimate_pmis.R`: Markov fitting, simulator bridge, and forward-simulation helpers.
- `R/plot_utils.R`: passage-tree plotting helpers.
- `R/workflow_utils.R`: connected-tree construction, interval grouping, diagnostics, and validation summaries.
- `ksim2.cpp`: Rcpp simulation engine compiled from R.
- `scripts/fit_group_markov.R`: grouped Markov grid fitting.
- `scripts/run_forward_validation.R`: forward validation for one group.
- `scripts/run_forward_validation_all.R`: forward validation for every best-fit group.

## Expected Commands

Grouped Markov fitting:

```bash
Rscript scripts/fit_group_markov.R \
  --grouped_intervals_path=core_data/grouped_intervals.Rds \
  --output_path=results/group_fit_df.Rds
```

Forward validation for one group:

```bash
Rscript scripts/run_forward_validation.R \
  --grouped_intervals_path=core_data/grouped_intervals.Rds \
  --group_fit_path=results/group_fit_df.Rds \
  --forward_group_id=<group_id> \
  --output_path=results/forward_validation_<group_id>.Rds \
  --n_reps=30 \
  --n_null_pairs=100 \
  --n_cores=4
```

Local laptop smoke check:

```bash
Rscript scripts/run_forward_validation.R \
  --grouped_intervals_path=core_data/grouped_intervals.Rds \
  --group_fit_path=results/group_fit_df.Rds \
  --forward_group_id=<group_id> \
  --output_path=results/local_smoke_forward_validation_<group_id>.Rds \
  --n_reps=3 \
  --n_null_pairs=10 \
  --n_cores=1 \
  --bottleneck_size=200 \
  --expansion_factor=8
```

Forward validation for all groups:

```bash
Rscript scripts/run_forward_validation_all.R \
  --grouped_intervals_path=core_data/grouped_intervals.Rds \
  --group_fit_path=results/group_fit_df.Rds \
  --output_dir=results \
  --summary_output_dir=result_summaries \
  --n_reps=30 \
  --n_null_pairs=100 \
  --n_cores=4
```

## Known Fragile Areas

- The local-running branch can drift from `main`; syncs should preserve laptop execution without silently changing scientific defaults.
- Forward validation is still a model-checking layer, not a settled final statistical endpoint.
- Ploidy diagnostics are secondary because grouped fitting targets marginal chromosome copy-number distributions, not full joint karyotypes.
- Generated `.Rds` files can be large and scientifically meaningful; do not regenerate or commit them casually.
- There is no formal automated test suite for simulator invariants, data schemas, connected-tree construction, or posterior-predictive summaries.
- Database-backed cache regeneration depends on `db_creds.txt`, which must remain untracked.

## Assumptions Not To Change Silently

- Karyotypes are modeled across chromosomes 1:22 after marker-chromosome exclusion.
- Copy-number state bounds are part of the model implementation.
- `pmis`/`p_mis` and `pwgd`/`p_wgd` grid definitions affect scientific conclusions.
- Forward-validation settings such as bottleneck size, expansion factor, replicate count, null pair count, core count, and recording frequency affect runtime and output interpretation.
- Reduced local settings are for smoke checks and debugging unless explicitly promoted by a maintainer.
- Wasserstein distance is the current posterior-predictive comparison metric.

## Open Questions

- Should local-running resource settings become a named preset or wrapper script instead of examples in documentation?
- Which generated outputs, if any, should be tracked on this local-running branch?
- Should this branch keep a local-only runner, or should all local runs use the `scripts/` entry points from `main`?
- Which model should be treated as the primary reference when the marginal Markov approximation and whole-cell forward simulator disagree?


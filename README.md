# Neutral Karyotype Evolution Simulator

This repository contains the current analysis pipeline for fitting and checking a neutral model of karyotype evolution under chromosome mis-segregation and whole-genome doubling. The main goal is to estimate grouped `p_mis` and `p_wgd` values from observed endpoint karyotypes, then compare those fitted parameters to forward simulations.

## Project overview

The project has three main layers:

1. `updated_workflow.Rmd`
   This is the main analysis document. It loads cached inputs, reconstructs connected passage trees, builds simulation intervals, loads grouped Markov fits, inspects diagnostics, and plots cached forward-validation results.
2. `R/`
   This folder contains the reusable analysis code:
   `R/db_utils.R` for data loading and preprocessing,
   `R/estimate_pmis.R` for Markov fitting and forward simulation helpers,
   `R/plot_utils.R` for tree plots,
   `R/workflow_utils.R` for workflow-level diagnostics and summaries.
3. `scripts/`
   This folder contains the expensive stages that are meant to run as cached jobs rather than inline in the notebook:
   `scripts/fit_group_markov.R`,
   `scripts/run_forward_validation.R`,
   `scripts/run_forward_validation_all.R`.

The simulator itself is implemented in [ksim2.cpp](/C:/Users/4473331/Documents/projects/021_ltee_models/NeutralKaryotypes/ksim2.cpp) and compiled from R with `Rcpp`.

## Repository data products

### `core_data/`

`core_data/` contains the core inputs and cached intermediate objects used by the workflow.

Important files:

- `karyotyped_samples.txt`
  The sample ids used as the observed karyotype set.
- `karyotypes.Rds`
  Cached endpoint karyotype vectors pulled from the database.
- `db_col.Rds`
  The collapsed passage-level ancestry graph derived from the raw passaging table.
- `media_raw.Rds`
  Cached media table used for condition annotation.
- `simulation_intervals.Rds`
  The interval objects built from connected passage trees and attached endpoint karyotypes.
- `grouped_intervals.Rds`
  The simulation intervals grouped by connected tree and resolved condition. This is the main input to grouped Markov fitting and forward validation.

### `results/`

`results/` contains heavier cached outputs.

Important files:

- `group_fit_df.Rds`
  The grouped grid-search output from the Markov fitting step.
- `forward_validation_<group>.Rds`
  Full forward-validation outputs for individual groups. These files contain the simulation bank and are relatively large.

### `result_summaries/`

`result_summaries/` contains lightweight summary sidecars for the forward-validation runs. These are intended for reporting and for the final summary chunk in the workflow.

Important files:

- `forward_validation_<group>_summary.Rds`
  Compact summaries extracted from the corresponding full forward-validation result.

Each summary file contains:

- `forward_group_id`
  The grouped interval id being validated.
- `fit_row`
  The selected best-fit parameter row used for the forward simulation.
- `settings`
  The simulation and posterior-predictive settings used to generate the result.
- `interval_summary`
  The main interval-level validation table.
- `replicate_summary`
  A higher-detail table summarizing the simulation bank for each interval.

The most useful columns in `interval_summary` are:

- `posterior_predictive_p`
  Tail probability from the posterior-predictive Wasserstein check. Smaller values indicate the observed endpoint looks less typical under the simulated null.
- `test_wasserstein`
  The median Wasserstein distance between the observed matched-size endpoint sample and the matched-size simulated endpoint samples.
- `null_wasserstein_median`
  The median Wasserstein distance between independently simulated matched-size endpoint samples.
- `null_wasserstein_q25`, `null_wasserstein_q75`
  A compact reference interval for the simulated null distances.

Interpretation:
If `test_wasserstein` is substantially larger than the null summary and `posterior_predictive_p` is small, the fitted forward model is not reproducing that endpoint interval well. If the observed statistic is similar to the simulated null, the forward behavior is more consistent with the fitted model.

## Current status

The project is now in a reasonably stable operational state.

- The workflow runs end-to-end from `updated_workflow.Rmd`.
- Expensive steps have been moved into scripts and cached to disk.
- Forward validation is now based on a posterior-predictive Wasserstein check rather than the earlier permutation-based check.
- Forward simulation is synchronized so one timestep corresponds to one division attempt per cell.
- The workflow now shows both ploidy diagnostics and marginal chromosome diagnostics, which better reflect the actual Markov fitting target.

What is still provisional:

- The forward-validation framework is useful, but it should still be treated as a model-checking step rather than a final statistical endpoint.
- The ploidy diagnostic is intentionally secondary to the marginal chromosome diagnostic, because the Markov fit is to marginals rather than full joint karyotypes.
- There is still real structural mismatch between the marginal Markov approximation and the whole-cell forward simulator.

## How to run

### Notebook

Open `updated_workflow.Rmd` in RStudio or another R Markdown environment and run the chunks sequentially.

### Grouped Markov fitting

```powershell
Rscript scripts/fit_group_markov.R --grouped_intervals_path=core_data/grouped_intervals.Rds --output_path=results/group_fit_df.Rds
```

### Forward validation for one group

```powershell
Rscript scripts/run_forward_validation.R --grouped_intervals_path=core_data/grouped_intervals.Rds --group_fit_path=results/group_fit_df.Rds --forward_group_id=<group_id> --output_path=results/forward_validation_<group_id>.Rds --n_reps=30 --n_null_pairs=100 --n_cores=4
```

### Forward validation for all groups

```powershell
Rscript scripts/run_forward_validation_all.R --grouped_intervals_path=core_data/grouped_intervals.Rds --group_fit_path=results/group_fit_df.Rds --output_dir=results --summary_output_dir=result_summaries --n_reps=30 --n_null_pairs=100 --n_cores=4
```

## Suggested next steps

- Parallelize the grouped Markov fitting path in `scripts/fit_group_markov.R`, since `estimate_group_pmis_grid()` is still one of the slowest stages.
- Add marginal chromosome diagnostic panels for more than two chromosomes, or make the selected chromosomes configurable from the workflow.
- Decide whether the whole-cell forward model or the marginal Markov approximation is the primary reference model, then tighten the mismatch between WGD and survival assumptions accordingly.
- Add a compact project-level report that aggregates all `result_summaries/` files into one stable summary artifact rather than rebuilding the table in the notebook each time.

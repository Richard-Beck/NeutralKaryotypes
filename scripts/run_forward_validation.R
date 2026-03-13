# Script helpers for forward validation. Safe to source from the workflow.

if (!requireNamespace("optparse", quietly = TRUE)) {
  stop("scripts/run_forward_validation.R requires the 'optparse' package.")
}

# Build a lightweight summary object from one forward-validation result.
# Input: full result list returned by `run_forward_validation_main()`.
# Output: compact list containing settings, fit row, and interval summaries only.
forward_validation_summary_object <- function(forward_result) {
  list(
    forward_group_id = forward_result$forward_group_id,
    fit_row = forward_result$fit_row,
    settings = forward_result$settings,
    interval_summary = forward_result$forward_validation$interval_summary,
    replicate_summary = forward_result$forward_validation$replicate_summary
  )
}

# Derive a sidecar summary path from a full forward-validation output path.
# Input: `output_path` for the full forward-validation RDS.
# Output: character path for the summary RDS.
forward_validation_summary_path <- function(output_path) {
  file.path(
    "result_summaries",
    basename(sub("[.]Rds$", "_summary.Rds", output_path))
  )
}

# Choose a conservative default worker count for local runs.
# Input: none.
# Output: integer core count for parallel forward simulation.
default_forward_n_cores <- function() {
  max(1L, parallel::detectCores(logical = FALSE) - 1L)
}

# Run forward simulations across interval/replicate tasks, optionally in parallel.
# Inputs: interval list, fit row, simulation settings, and worker count.
# Output: nested list matching `simulate_group_forward()`.
run_parallel_forward_runs <- function(group_intervals,
                                      fit_row,
                                      bottleneck_size,
                                      expansion_factor,
                                      n_reps,
                                      record_every,
                                      seed,
                                      n_cores = 1L) {
  n_cores <- max(1L, as.integer(n_cores))
  if (n_cores <= 1L) {
    return(simulate_group_forward(
      group_intervals,
      fit_row = fit_row,
      bottleneck_size = bottleneck_size,
      expansion_factor = expansion_factor,
      n_reps = n_reps,
      seed = seed,
      record_every = record_every
    ))
  }

  tasks <- expand.grid(
    interval_idx = seq_along(group_intervals),
    rep_idx = seq_len(n_reps),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  cl <- parallel::makeCluster(n_cores)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  parallel::clusterCall(cl, setwd, getwd())
  parallel::clusterEvalQ(cl, {
    source("R/estimate_pmis.R")
    NULL
  })
  parallel::clusterExport(
    cl,
    varlist = c("group_intervals", "fit_row", "bottleneck_size", "expansion_factor", "record_every", "seed", "tasks"),
    envir = environment()
  )

  task_results <- parallel::parLapply(cl, seq_len(nrow(tasks)), function(task_idx) {
    task <- tasks[task_idx, ]
    simulate_interval_forward(
      interval = group_intervals[[task$interval_idx]],
      p_mis = fit_row$pmis,
      p_wgd = fit_row$pwgd,
      bottleneck_size = bottleneck_size,
      expansion_factor = expansion_factor,
      seed = seed + task$rep_idx + 1000L * task$interval_idx,
      record_every = record_every
    )
  })

  split_results <- split(task_results, tasks$interval_idx)
  lapply(split_results, function(interval_runs) {
    names(interval_runs) <- NULL
    interval_runs
  })
}

# Run forward simulations plus posterior-predictive validation and save the result.
# Inputs: grouped intervals path, group-fit path, forward group id, output path, and simulation settings.
# Output: list containing the chosen fit row, forward runs, and forward validation; also written to `output_path`.
run_forward_validation_main <- function(grouped_intervals_path = "core_data/grouped_intervals.Rds",
                                        group_fit_path = "results/group_fit_df.Rds",
                                        forward_group_id,
                                        output_path = "results/forward_validation.Rds",
                                        summary_output_path = forward_validation_summary_path(output_path),
                                        bottleneck_size = 1000L,
                                        expansion_factor = 32L,
                                        n_reps = 30L,
                                        record_every = 50L,
                                        n_null_pairs = 100L,
                                        n_cores = default_forward_n_cores(),
                                        seed = 1L) {
  if (missing(forward_group_id) || !nzchar(forward_group_id)) {
    stop("run_forward_validation_main() requires a non-empty `forward_group_id`.")
  }

  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(summary_output_path), recursive = TRUE, showWarnings = FALSE)

  source("R/estimate_pmis.R")
  source("R/workflow_utils.R")

  grouped_intervals <- readRDS(grouped_intervals_path)
  group_fit_df <- readRDS(group_fit_path)
  best_group_pmis <- select_best_group_fits(group_fit_df)
  fit_row <- best_group_pmis[best_group_pmis$group_id == forward_group_id, , drop = FALSE]
  if (!nrow(fit_row)) {
    stop("Could not find forward_group_id in best group fits: ", forward_group_id)
  }
  fit_row <- fit_row[1, , drop = FALSE]

  forward_runs <- run_parallel_forward_runs(
    grouped_intervals[[forward_group_id]],
    fit_row = fit_row,
    bottleneck_size = as.integer(bottleneck_size),
    expansion_factor = as.integer(expansion_factor),
    n_reps = as.integer(n_reps),
    record_every = as.integer(record_every),
    seed = as.integer(seed),
    n_cores = as.integer(n_cores)
  )

  forward_validation <- summarize_forward_validation(
    forward_runs = forward_runs,
    group_intervals = grouped_intervals[[forward_group_id]],
    group_id = forward_group_id,
    n_null_pairs = as.integer(n_null_pairs),
    seed_base = as.integer(seed) * 10000L
  )

  out <- list(
    forward_group_id = forward_group_id,
    fit_row = fit_row,
    settings = list(
      bottleneck_size = as.integer(bottleneck_size),
      expansion_factor = as.integer(expansion_factor),
      n_reps = as.integer(n_reps),
      record_every = as.integer(record_every),
      n_null_pairs = as.integer(n_null_pairs),
      n_cores = as.integer(n_cores),
      seed = as.integer(seed)
    ),
    forward_runs = forward_runs,
    forward_validation = forward_validation
  )

  saveRDS(out, output_path)
  saveRDS(forward_validation_summary_object(out), summary_output_path)
  message("Saved forward validation to ", output_path)
  message("Saved forward validation summary to ", summary_output_path)
  out
}

# Load cached forward-validation results for the workflow, computing them if absent.
# Inputs: forward group id, grouped interval path, fit path, output path, and simulation settings.
# Output: cached forward-validation result list.
load_forward_validation_result <- function(forward_group_id,
                                           grouped_intervals_path = "core_data/grouped_intervals.Rds",
                                           group_fit_path = "results/group_fit_df.Rds",
                                           output_path,
                                           summary_output_path = forward_validation_summary_path(output_path),
                                           bottleneck_size = 1000L,
                                           expansion_factor = 32L,
                                           n_reps = 30L,
                                           record_every = 50L,
                                           n_null_pairs = 100L,
                                           n_cores = default_forward_n_cores(),
                                           seed = 1L) {
  if (missing(forward_group_id) || !nzchar(forward_group_id)) {
    stop("load_forward_validation_result() requires a non-empty `forward_group_id`.")
  }
  if (missing(output_path) || !nzchar(output_path)) {
    stop("load_forward_validation_result() requires `output_path`.")
  }

  if (!file.exists(output_path)) {
    run_forward_validation_main(
      grouped_intervals_path = grouped_intervals_path,
      group_fit_path = group_fit_path,
      forward_group_id = forward_group_id,
      output_path = output_path,
      summary_output_path = summary_output_path,
      bottleneck_size = bottleneck_size,
      expansion_factor = expansion_factor,
      n_reps = n_reps,
      record_every = record_every,
      n_null_pairs = n_null_pairs,
      n_cores = n_cores,
      seed = seed
    )
  }

  readRDS(output_path)
}

if (sys.nframe() == 0L) {
  option_list <- list(
    optparse::make_option("--grouped_intervals_path", type = "character", default = "core_data/grouped_intervals.Rds"),
    optparse::make_option("--group_fit_path", type = "character", default = "results/group_fit_df.Rds"),
    optparse::make_option("--forward_group_id", type = "character", default = NULL),
    optparse::make_option("--output_path", type = "character", default = "results/forward_validation.Rds"),
    optparse::make_option("--summary_output_path", type = "character", default = NULL),
    optparse::make_option("--bottleneck_size", type = "integer", default = 1000L),
    optparse::make_option("--expansion_factor", type = "integer", default = 32L),
    optparse::make_option("--n_reps", type = "integer", default = 30L),
    optparse::make_option("--record_every", type = "integer", default = 50L),
    optparse::make_option("--n_null_pairs", type = "integer", default = 100L),
    optparse::make_option("--n_cores", type = "integer", default = default_forward_n_cores()),
    optparse::make_option("--seed", type = "integer", default = 1L)
  )
  opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))

  if (is.null(opt$forward_group_id) || !nzchar(opt$forward_group_id)) {
    stop("scripts/run_forward_validation.R requires --forward_group_id=<group id>.")
  }

  summary_output_path <- opt$summary_output_path
  if (is.null(summary_output_path) || !nzchar(summary_output_path)) {
    summary_output_path <- forward_validation_summary_path(opt$output_path)
  }

  run_forward_validation_main(
    grouped_intervals_path = opt$grouped_intervals_path,
    group_fit_path = opt$group_fit_path,
    forward_group_id = opt$forward_group_id,
    output_path = opt$output_path,
    summary_output_path = summary_output_path,
    bottleneck_size = opt$bottleneck_size,
    expansion_factor = opt$expansion_factor,
    n_reps = opt$n_reps,
    record_every = opt$record_every,
    n_null_pairs = opt$n_null_pairs,
    n_cores = opt$n_cores,
    seed = opt$seed
  )
}

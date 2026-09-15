# Run from the NeutralKaryotypes root. Tiny smoke fixtures, not scientific runs.
source("R/estimate_pmis.R")
source("R/forward_checkpoint_utils.R")
source("scripts/run_forward_validation.R")

expect_error <- function(expr, pattern) {
  caught <- tryCatch({ force(expr); NULL }, error = function(e) conditionMessage(e))
  if (is.null(caught) || !grepl(pattern, caught, fixed = TRUE)) {
    stop("Expected error containing '", pattern, "', received: ", caught)
  }
}

fixture <- list(
  start = "reviewed_observation_start", end = "reviewed_observation_end",
  start_karyotype = matrix(c(2, 2, 2, 3), ncol = 2L),
  end_karyotype = matrix(c(2, 3, 2, 3), ncol = 2L),
  n_divisions = 2L, conditions = c(control = 1L, glucose = 1L),
  physical_history = c("reviewed_cycle_1", "reviewed_cycle_2"),
  ontology_generation = "explicit_test_fixture"
)
intervals <- list(fixture, fixture)
intervals[[2L]]$end <- "reviewed_observation_end_2"
fit <- data.frame(group_id = "synthetic_smoke", pmis = 0.02, pwgd = 0, n_members = 2L)
test_root <- tempfile("neutral-forward-checkpoint-smoke-")
dir.create(test_root)
# Every output is confined to this newly created tiny test directory.
on.exit_cleanup <- function() unlink(test_root, recursive = TRUE)

settings <- list(group_intervals = intervals, fit_row = fit, bottleneck_size = 4L,
                 expansion_factor = 2L, n_reps = 3L, record_every = 1L, seed = 17L)
run_checkpoints <- function(path, ...) {
  do.call(run_checkpointed_forward_runs, c(settings, list(checkpoint_dir = path), list(...)))
}

serial_path <- file.path(test_root, "serial")
baseline <- do.call(simulate_group_forward, settings)
serial <- run_checkpoints(serial_path, n_cores = 1L, batch_size = 2L)
stopifnot(all(!attr(serial, "checkpoint_receipts")$reused))
attr(serial, "checkpoint_receipts") <- NULL
stopifnot(identical(serial, baseline))
message("PASS: native serial simulator parity, seed sequence, and lossless checkpoint readback.")

original_simulate <- simulate_interval_forward
simulate_interval_forward <- function(...) stop("A validated task was unexpectedly simulated again.")
resumed <- run_checkpoints(serial_path, n_cores = 1L, batch_size = 1L)
stopifnot(all(attr(resumed, "checkpoint_receipts")$reused))
attr(resumed, "checkpoint_receipts") <- NULL
stopifnot(identical(resumed, baseline))
simulate_interval_forward <- original_simulate
message("PASS: full resume performs zero simulations, despite changed batch size.")

subset_path <- file.path(test_root, "subset")
partial <- run_checkpoints(subset_path, interval_indices = 2L, replicate_indices = 1L, load_results = FALSE)
stopifnot(!partial$complete, nrow(partial$receipts) == 1L,
          partial$receipts$seed == 2018L,
          identical(load_forward_checkpoint_runs(subset_path, 2L, 1L), list(baseline[[2L]][1L])))
completed <- run_checkpoints(subset_path, n_cores = 2L, batch_size = 2L)
stopifnot(sum(attr(completed, "checkpoint_receipts")$reused) == 1L)
attr(completed, "checkpoint_receipts") <- NULL
stopifnot(identical(completed, baseline))
message("PASS: partial scheduling preserves native interval indices and resumes with parallel workers.")

failure_path <- file.path(test_root, "failure")
simulate_interval_forward <- function(..., seed) {
  if (seed == 1019L) stop("Injected later-task failure")
  original_simulate(..., seed = seed)
}
expect_error(run_checkpoints(failure_path, n_cores = 1L, batch_size = 3L), "Injected later-task failure")
stopifnot(nrow(forward_checkpoint_receipts(failure_path)) == 1L)
simulate_interval_forward <- original_simulate
repaired <- run_checkpoints(failure_path, n_cores = 1L, load_results = FALSE)
stopifnot(repaired$complete, sum(repaired$receipts$reused) == 1L)
message("PASS: each completed task survives a later task failure within its batch.")

changed <- settings
changed$group_intervals[[1L]]$end_karyotype[1L, 1L] <- 4
expect_error(do.call(run_checkpointed_forward_runs, c(changed, list(checkpoint_dir = serial_path))), "identity mismatch")
changed <- settings
changed$group_intervals[[1L]]$physical_history[1L] <- "different_reviewed_cycle"
expect_error(do.call(run_checkpointed_forward_runs, c(changed, list(checkpoint_dir = serial_path))), "identity mismatch")
changed <- settings
changed$fit_row$pmis <- 0.03
expect_error(do.call(run_checkpointed_forward_runs, c(changed, list(checkpoint_dir = serial_path))), "identity mismatch")
changed <- settings
changed$seed <- 18L
expect_error(do.call(run_checkpointed_forward_runs, c(changed, list(checkpoint_dir = serial_path))), "identity mismatch")
changed <- settings
changed$record_every <- 2L
expect_error(do.call(run_checkpointed_forward_runs, c(changed, list(checkpoint_dir = serial_path))), "identity mismatch")
message("PASS: matrices, full physical history, fit parameters, seed, and simulation settings fail closed on changes.")

task_path <- forward_checkpoint_task_path(failure_path, 1L, 1L)
task <- readRDS(task_path)
task$result$final_matrix[1L, 1L] <- task$result$final_matrix[1L, 1L] + 1
saveRDS(task, task_path)
expect_error(run_checkpoints(failure_path, load_results = FALSE), "payload checksum mismatch")
message("PASS: corrupted completed results are not silently reused or regenerated.")

# Test explicit output routing without transport or any scientific calculation.
# These local stubs operate only on a copied entrypoint environment, never owners.
route_env <- new.env(parent = globalenv())
route_env$source <- function(...) invisible(NULL)
route_env$select_best_group_fits <- function(x) x
route_env$run_parallel_forward_runs <- function(...) list()
route_env$summarize_forward_validation <- function(...) list(interval_summary = data.frame(), replicate_summary = data.frame())
route_env$aggregate_called <- FALSE
route_env$write_neutral_probability_csv <- function(...) route_env$aggregate_called <- TRUE
route_runner <- run_forward_validation_main
environment(route_runner) <- route_env
saveRDS(list(synthetic_smoke = intervals), file.path(test_root, "groups.Rds"))
saveRDS(fit, file.path(test_root, "fit.Rds"))
old_wd <- getwd()
setwd(test_root)
invisible(route_runner(grouped_intervals_path = "groups.Rds", group_fit_path = "fit.Rds",
             forward_group_id = "synthetic_smoke", output_path = "routed.Rds",
             summary_output_path = "result_summaries/explicit_summary.Rds"))
stopifnot(!route_env$aggregate_called)
invisible(route_runner(grouped_intervals_path = "groups.Rds", group_fit_path = "fit.Rds",
             forward_group_id = "synthetic_smoke", output_path = "default.Rds"))
stopifnot(route_env$aggregate_called)
setwd(old_wd)
message("PASS: explicit summary paths never rewrite the canonical aggregate by default; omitted paths preserve legacy behavior.")

on.exit_cleanup()
message("All forward checkpoint smoke tests passed.")

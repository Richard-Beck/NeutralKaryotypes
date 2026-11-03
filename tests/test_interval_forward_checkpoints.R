# Tiny synthetic contracts only; never read or regenerate scientific outputs.
source("R/estimate_pmis.R")
source("R/interval_fit_utils.R")
source("R/workflow_utils.R")
source("R/forward_checkpoint_utils.R")
source("scripts/run_forward_validation.R")
source("scripts/run_forward_validation_all.R")

expect_error <- function(expr, pattern) {
  result <- tryCatch({ force(expr); NULL }, error = function(e) conditionMessage(e))
  if (is.null(result) || !grepl(pattern, result, fixed = TRUE)) stop("Expected '", pattern, "', received: ", result)
}

test_forward_default_routing <- function() {
  native_root <- getwd()
  route_root <- tempfile("interval-forward-routing-")
  dir.create(route_root)
  on.exit({ setwd(native_root); unlink(route_root, recursive = TRUE) }, add = TRUE)
  new_mode <- "interval_shared_wgd_v1"
  expected <- list(group_fit_path = "results/interval_fit_grid.Rds",
                   output_dir = "results/interval_forward_validation",
                   output_path = "results/interval_forward_validation/forward_validation.Rds",
                   summary_output_dir = "result_summaries/interval_forward_validation",
                   aggregate_output_path = "data/neutral_interval_diagnostics.csv")
  stopifnot(identical(forward_validation_default_paths(new_mode), expected),
            identical(forward_validation_summary_path("custom/example.Rds", new_mode),
                      "result_summaries/interval_forward_validation/example_summary.Rds"),
            identical(forward_validation_summary_path("custom/example.Rds"), "result_summaries/example_summary.Rds"))

  # Exercise the real entrypoint routing with isolated data/simulation stubs.
  # Sentinels prove that omitted interval paths never rewrite pooled artifacts.
  route <- new.env(parent = globalenv())
  route$source <- function(...) invisible(NULL)
  route$reads <- character()
  route$simulation_args <- NULL
  route$aggregate_args <- NULL
  route$readRDS <- function(path) {
    route$reads <- c(route$reads, path)
    if (path == "core_data/grouped_intervals.Rds") return(list(group_a = list(list())))
    if (path %in% c("results/group_fit_df.Rds", "results/interval_fit_grid.Rds", "custom/input_fit.Rds")) {
      return(data.frame(group_id = "group_a", pmis = 0.01, pwgd = 0, negll = 1,
                        estimability_status = "unique_grid_optimum"))
    }
    base::readRDS(path)
  }
  route$select_best_group_fits <- function(x) x
  route$select_forward_fits <- function(x, ...) x
  route$nk_interval_roster <- function(...) data.frame(physical_interval_id = "test_physical_interval")
  route$forward_rng_roster <- function(...) data.frame(seed = 1L)
  route$run_parallel_forward_runs <- function(...) {
    route$simulation_args <- list(...)
    list(list())
  }
  route$summarize_forward_validation <- function(...) {
    list(interval_summary = data.frame(), replicate_summary = data.frame())
  }
  route$forward_checkpoint_atomic_save <- function(object, path, ...) saveRDS(object, path)
  route$write_neutral_probability_csv <- function(...) {
    route$aggregate_args <- list(...)
    dir.create(dirname(route$aggregate_args$output_path), recursive = TRUE, showWarnings = FALSE)
    writeLines("synthetic routing check", route$aggregate_args$output_path)
  }
  single <- run_forward_validation_main
  environment(single) <- route
  route$run_forward_validation_main <- single
  loader <- load_forward_validation_result
  environment(loader) <- route
  all_groups <- run_forward_validation_all_main
  environment(all_groups) <- route
  setwd(route_root)
  legacy_paths <- c("results/group_fit_df.Rds", "results/forward_validation.Rds", "results/forward_validation_group_a.Rds",
                    "result_summaries/forward_validation_summary.Rds", "result_summaries/forward_validation_group_a_summary.Rds",
                    "data/neutral_probability.csv")
  for (path in legacy_paths) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    writeLines("frozen pooled sentinel", path)
  }
  legacy_hashes <- tools::md5sum(legacy_paths)
  invisible(single(forward_group_id = "group_a", fit_mode = new_mode, n_cores = 1L))
  stopifnot(expected$group_fit_path %in% route$reads, !"results/group_fit_df.Rds" %in% route$reads,
            file.exists(expected$output_path),
            file.exists(file.path(expected$summary_output_dir, "forward_validation_summary.Rds")),
            identical(route$simulation_args$checkpoint_dir, paste0(expected$output_path, ".checkpoints")),
            is.null(route$aggregate_args))
  invisible(single(forward_group_id = "group_a", fit_mode = new_mode, n_cores = 1L, update_aggregate = TRUE))
  stopifnot(identical(route$aggregate_args$output_path, expected$aggregate_output_path))
  invisible(loader(forward_group_id = "group_a", fit_mode = new_mode, n_cores = 1L))
  expect_error(loader(forward_group_id = "group_a", n_cores = 1L), "requires `output_path`")
  invisible(single(forward_group_id = "group_a", fit_mode = new_mode, n_cores = 1L,
                   group_fit_path = "custom/input_fit.Rds", output_path = "custom/output.Rds",
                   summary_output_path = "custom/summary.Rds", checkpoint_dir = "custom/checkpoints"))
  invisible(loader(forward_group_id = "group_a", fit_mode = new_mode, n_cores = 1L,
                   group_fit_path = "custom/input_fit.Rds", output_path = "custom/cache.Rds",
                   summary_output_path = "custom/cache_summary.Rds", checkpoint_dir = "custom/cache_checkpoints"))
  stopifnot(file.exists("custom/output.Rds"), file.exists("custom/summary.Rds"),
            file.exists("custom/cache.Rds"), file.exists("custom/cache_summary.Rds"),
            identical(route$simulation_args$checkpoint_dir, "custom/cache_checkpoints"))
  default_all <- all_groups(fit_mode = new_mode, n_cores = 1L)
  stopifnot(identical(default_all, file.path(expected$output_dir, "forward_validation_group_a.Rds")),
            file.exists(file.path(expected$summary_output_dir, "forward_validation_group_a_summary.Rds")),
            identical(route$aggregate_args$output_path, expected$aggregate_output_path),
            identical(route$simulation_args$checkpoint_dir, file.path(expected$output_dir, "forward_checkpoints_interval_v1")))
  invisible(all_groups(fit_mode = new_mode, n_cores = 1L, group_fit_path = "custom/input_fit.Rds",
                       output_dir = "custom/all", summary_output_dir = "custom/all_summaries",
                       aggregate_output_path = "custom/diagnostics.csv"))
  stopifnot(file.exists("custom/all/forward_validation_group_a.Rds"),
            file.exists("custom/all_summaries/forward_validation_group_a_summary.Rds"),
            identical(route$aggregate_args$output_path, "custom/diagnostics.csv"),
            identical(tools::md5sum(legacy_paths), legacy_hashes))
  # Explicitly routed summaries in the legacy directory cannot cause an
  # implicit write to the legacy aggregate, even when update is requested.
  route$aggregate_args <- NULL
  invisible(single(forward_group_id = "group_a", fit_mode = new_mode, n_cores = 1L,
                   summary_output_path = "result_summaries/explicit_interval_summary.Rds", update_aggregate = TRUE))
  stopifnot(is.null(route$aggregate_args), identical(tools::md5sum(legacy_paths), legacy_hashes))
  invisible(single(forward_group_id = "group_a", n_cores = 1L))
  stopifnot(identical(route$aggregate_args$output_path, "data/neutral_probability.csv"),
            "results/group_fit_df.Rds" %in% route$reads)
  pooled_all <- all_groups(n_cores = 1L)
  stopifnot(identical(pooled_all, "results/forward_validation_group_a.Rds"),
            identical(route$aggregate_args$summary_dir, "result_summaries"),
            identical(route$aggregate_args$output_path, "data/neutral_probability.csv"))

  writer_env <- new.env(parent = globalenv())
  writer_env$write_interval_neutral_probability_csv <- function(summary_dir, grouped_intervals_path, output_path, summary_pattern, fit_mode) {
    list(summary_dir = summary_dir, output_path = output_path)
  }
  writer <- write_neutral_probability_csv
  environment(writer) <- writer_env
  stopifnot(identical(writer(fit_mode = new_mode), list(summary_dir = expected$summary_output_dir, output_path = expected$aggregate_output_path)),
            identical(writer(summary_dir = "custom/summaries", output_path = "custom/export.csv", fit_mode = new_mode),
                      list(summary_dir = "custom/summaries", output_path = "custom/export.csv")))

  # Evaluate the actual CLI blocks with supplied parser arguments and captured
  # entrypoint calls. No subprocess or scientific work is needed for this test.
  capture_cli <- function(path, args, entrypoint) {
    expression <- tail(parse(file.path(native_root, path)), 1L)[[1L]][[3L]]
    parser_line <- which(vapply(as.list(expression), function(x) {
      is.call(x) && identical(x[[1L]], as.name("<-")) && identical(x[[2L]], as.name("opt"))
    }, logical(1)))
    stopifnot(length(parser_line) == 1L)
    expression[[parser_line]][[3L]]$args <- as.name("cli_args")
    cli <- new.env(parent = globalenv())
    cli$cli_args <- args
    cli[[entrypoint]] <- function(...) list(...)
    eval(expression, cli)
  }
  for (mode in c("pooled_group_v1", new_mode, "interval_independent_v1")) {
    defaults <- forward_validation_default_paths(mode)
    single_cli <- capture_cli("scripts/run_forward_validation.R", c(paste0("--fit_mode=", mode), "--forward_group_id=group_a", "--n_cores=1"), "run_forward_validation_main")
    all_cli <- capture_cli("scripts/run_forward_validation_all.R", c(paste0("--fit_mode=", mode), "--n_cores=1"), "run_forward_validation_all_main")
    stopifnot(identical(single_cli$group_fit_path, defaults$group_fit_path), identical(single_cli$output_path, defaults$output_path),
              identical(single_cli$summary_output_path, forward_validation_summary_path(defaults$output_path, mode)),
              identical(single_cli$update_aggregate, mode == "pooled_group_v1"),
              !"calibration_status" %in% names(single_cli),
              identical(all_cli$group_fit_path, defaults$group_fit_path), identical(all_cli$output_dir, defaults$output_dir),
              identical(all_cli$summary_output_dir, defaults$summary_output_dir), identical(all_cli$aggregate_output_path, defaults$aggregate_output_path))
  }
  single_cli <- capture_cli("scripts/run_forward_validation.R", c(paste0("--fit_mode=", new_mode), "--forward_group_id=group_a", "--group_fit_path=custom/fit.Rds", "--output_path=custom/run.Rds", "--summary_output_path=custom/summary.Rds", "--checkpoint_dir=custom/checkpoints"), "run_forward_validation_main")
  all_cli <- capture_cli("scripts/run_forward_validation_all.R", c(paste0("--fit_mode=", new_mode), "--group_fit_path=custom/fit.Rds", "--output_dir=custom/results", "--summary_output_dir=custom/summaries", "--aggregate_output_path=custom/diagnostics.csv"), "run_forward_validation_all_main")
  stopifnot(identical(single_cli$group_fit_path, "custom/fit.Rds"), identical(single_cli$output_path, "custom/run.Rds"),
            identical(single_cli$summary_output_path, "custom/summary.Rds"), identical(single_cli$checkpoint_dir, "custom/checkpoints"),
            identical(all_cli$group_fit_path, "custom/fit.Rds"), identical(all_cli$output_dir, "custom/results"),
            identical(all_cli$summary_output_dir, "custom/summaries"), identical(all_cli$aggregate_output_path, "custom/diagnostics.csv"))
  # Capture the approved production arguments without launching simulations.
  production <- c("--fit_mode=interval_independent_v1", "--bottleneck_size=3000", "--expansion_factor=32",
                  "--record_every=50", "--n_reps=200", "--n_null_pairs=1000", "--distance_metric=wasserstein")
  single_production <- capture_cli("scripts/run_forward_validation.R", c(production, "--forward_group_id=group_a"), "run_forward_validation_main")
  all_production <- capture_cli("scripts/run_forward_validation_all.R", production, "run_forward_validation_all_main")
  fields <- c("bottleneck_size", "expansion_factor", "record_every", "n_reps", "n_null_pairs", "distance_metric")
  expected_settings <- list(bottleneck_size = 3000L, expansion_factor = 32L, record_every = 50L,
                            n_reps = 200L, n_null_pairs = 1000L, distance_metric = "wasserstein")
  stopifnot(identical(single_production[fields], expected_settings), identical(all_production[fields], expected_settings))
  message("PASS: interval function/loader/CLI/writer defaults preserve pooled artifacts; explicit paths and pooled defaults remain unchanged.")
}

test_forward_default_routing()
if ("--routing-only" %in% commandArgs(trailingOnly = TRUE)) {
  message("Focused interval-forward routing tests passed; simulator tests were not repeated.")
  quit(status = 0L)
}

fixture <- function(end, end_value = 2) list(
  connected_tree = "synthetic_exact_tree", start = "physical_start", end = end,
  start_ids = "physical_start_sample", end_ids = paste0(end, "_sample"),
  start_alias = NULL, n_collapsed_passages = 0L, n_divisions = 2L, elapsed_passages = 1L,
  start_karyotype = matrix(c(2, 2, 2, 3), ncol = 2L),
  end_karyotype = matrix(c(end_value, 2, 2, 3), ncol = 2L),
  conditions = c(control = 1), fixture_scope = "synthetic_smoke_only")

test_independent_forward_contract <- function() {
  test_root <- tempfile("independent-forward-contract-")
  dir.create(test_root)
  original_nll <- get_interval_nll_markov
  original_simulator <- simulate_interval_forward
  on.exit({
    assign("get_interval_nll_markov", original_nll, envir = globalenv())
    assign("simulate_interval_forward", original_simulator, envir = globalenv())
    unlink(test_root, recursive = TRUE)
  }, add = TRUE)
  assign("get_interval_nll_markov", function(p_mis, p_wgd, K0, KT, n_div, states = 0:8, eps = 1e-9) {
    if (KT[1L, 1L] == 5) return(Inf)
    target_p <- if (KT[1L, 1L] == 2) 0.001 else 0.04
    target_q <- if (KT[1L, 1L] == 2) 0 else 0.02
    1e5 * ((p_mis - target_p)^2 + (p_wgd - target_q)^2)
  }, envir = globalenv())
  assign("simulate_interval_forward", function(interval, p_mis, p_wgd, bottleneck_size,
                                                expansion_factor, record_every, divide_all_each_step, seed) {
    set.seed(seed)
    list(interval = interval, p_mis = p_mis, p_wgd = p_wgd,
         final_matrix = interval$start_karyotype[sample(seq_len(nrow(interval$start_karyotype)), 4L, replace = TRUE), , drop = FALSE])
  }, envir = globalenv())
  mode <- "interval_independent_v1"
  groups <- list(group_a = list(fixture("independent_end_a"), fixture("independent_end_b", 3)))
  fit <- function(grouped, fit_mode = mode) nk_select_interval_fits(nk_interval_fit_grid(
    grouped, c(0.001, 0.04), c(0, 0.02), fit_mode = fit_mode))
  fits <- fit(groups)
  shared <- fit(groups, "interval_shared_wgd_v1")
  ids <- nk_interval_roster(groups)$physical_interval_id
  fits <- fits[match(ids, fits$physical_interval_id), , drop = FALSE]
  stopifnot(identical(fits$pwgd, c(0, 0.02)), identical(fits$pmis, c(0.001, 0.04)),
            all(fits$wgd_scope == "physical_interval"))
  seeds <- forward_rng_roster(ids, 3L, 1L)
  raw_root <- file.path(test_root, "checkpoints")
  run <- function(selected, grouped = groups, fit_mode = mode) run_checkpointed_forward_runs(
    grouped[[1L]], selected, 4L, 2L, 3L, 1L, 1L, checkpoint_dir = raw_root,
    fit_mode = fit_mode, seed_roster = seeds)
  shared_runs <- run(shared, fit_mode = "interval_shared_wgd_v1")
  independent_runs <- run(fits)
  stopifnot(sum(attr(independent_runs, "checkpoint_receipts")$reused) == 3L)
  for (i in seq_along(independent_runs)) {
    stopifnot(all(vapply(independent_runs[[i]], function(x) identical(x$p_mis, fits$pmis[i]) && identical(x$p_wgd, fits$pwgd[i]), logical(1))))
  }
  expect_error(run(shared), "requested fit mode")
  summary <- function(runs, selected, grouped = groups) summarize_forward_validation(
    runs, grouped[[1L]], names(grouped)[1L], n_null_pairs = 2L, distance_metric = "wasserstein",
    fit_mode = mode, fit_rows = selected, seed_roster = seeds, checkpoint_dir = raw_root)
  checked <- summary(independent_runs, fits)
  resumed <- summary(independent_runs, fits)
  stopifnot(all(resumed$interval_summary$check_method == "posterior_predictive_check_v1"),
            identical(resumed$interval_summary$p_value, resumed$interval_summary$posterior_predictive_p),
            all(is.finite(resumed$interval_summary$p_value)),
            !any(c("calibration_status", "inference_eligible", "diagnostic_p", "q_value", "neutrality_rejected") %in% names(resumed$interval_summary)),
            !"calibration" %in% names(resumed),
            all(vapply(resumed$interval_checks, function(x) isTRUE(attr(x, "checkpoint_reused")), logical(1))))
  for (i in seq_along(independent_runs)) {
    direct <- posterior_predictive_endpoint_check(groups[[1L]][[i]]$end_karyotype,
      lapply(independent_runs[[i]], `[[`, "final_matrix"), n_pairs = 2L, distance_metric = "wasserstein",
      subsample_seeds = forward_rng_lookup(seeds, ids[i], "endpoint_subsample", seq_along(independent_runs[[i]])),
      null_pair_seed = forward_rng_lookup(seeds, ids[i], "endpoint_null_pairs", 0L))
    actual <- checked$interval_checks[[i]]
    stopifnot(identical(actual$observed_distances, direct$observed_distances),
              identical(actual$null_distances, direct$null_distances),
              checked$interval_summary$test_wasserstein[i] == median(direct$observed_distances),
              checked$interval_summary$p_value[i] == mean(direct$null_distances >= median(direct$observed_distances)))
  }
  changed_groups <- groups
  changed_groups$group_a[[2L]]$end_karyotype[2L, 1L] <- 4
  changed_fits <- fit(changed_groups)
  changed_runs <- run(changed_fits, changed_groups)
  changed_summary <- summary(changed_runs, changed_fits, changed_groups)
  stopifnot(all(attr(changed_runs, "checkpoint_receipts")$reused),
            identical(checked$interval_summary$fit_dependency_id[1L], changed_summary$interval_summary$fit_dependency_id[1L]),
            isTRUE(attr(changed_summary$interval_checks[[1L]], "checkpoint_reused")))
  summary_dir <- file.path(test_root, "summaries")
  dir.create(summary_dir)
  grouped_path <- file.path(test_root, "groups.Rds")
  saveRDS(groups, grouped_path)
  result <- list(forward_group_id = "group_a", fit_row = fits, settings = list(fit_mode = mode, check_method = "posterior_predictive_check_v1"),
                 forward_validation = resumed, seed_roster = seeds)
  saveRDS(forward_validation_summary_object(result), file.path(summary_dir, "forward_validation_group_a_summary.Rds"))
  csv <- write_neutral_probability_csv(summary_dir, grouped_path, file.path(test_root, "diagnostics.csv"), fit_mode = mode)
  stopifnot(setequal(csv$pwgd, c(0, 0.02)), all(csv$check_method == "posterior_predictive_check_v1"),
            all(is.finite(csv$p_value)), identical(csv$p_value, csv$posterior_predictive_p))
  expect_error(write_neutral_probability_csv(summary_dir, grouped_path, file.path(test_root, "wrong_mode.csv"), fit_mode = "interval_shared_wgd_v1"), "rejects legacy, mixed")
  unavailable_groups <- groups
  unavailable_groups$group_a[[2L]]$end_karyotype[1L, 1L] <- 5
  unavailable_fits <- fit(unavailable_groups)
  unavailable_runs <- run(unavailable_fits, unavailable_groups)
  unavailable_summary <- summary(unavailable_runs, unavailable_fits, unavailable_groups)
  stopifnot(length(unavailable_runs[[1L]]) == 3L, length(unavailable_runs[[2L]]) == 0L,
            unavailable_summary$interval_summary$forward_status[2L] == "not_simulated_non_estimable_fit",
            is.finite(unavailable_summary$interval_summary$p_value[1L]), is.na(unavailable_summary$interval_summary$p_value[2L]))
  message("PASS: independent own-pair dispatch, original posterior-predictive algebra, checkpoint reuse and explicit mixed estimability.")
}

test_independent_forward_contract()
if ("--independent-only" %in% commandArgs(trailingOnly = TRUE)) {
  message("Focused independent forward contract tests passed; no production or statistical validation was run.")
  quit(status = 0L)
}
groups <- list(group_a = list(fixture("physical_end_a"), fixture("physical_end_b", 3)))
test_root <- tempfile("interval-forward-smoke-")
dir.create(test_root)
original_nll <- get_interval_nll_markov
original_simulator <- simulate_interval_forward
original_runtime <- forward_checkpoint_runtime

# Construct known, different interval optima and one shared nonzero WGD optimum.
# The likelihood kernel itself is covered by the fitter's numerical tests.
get_interval_nll_markov <- function(p_mis, p_wgd, K0, KT, n_div, states = 0:8, eps = 1e-9) {
  target <- if (KT[1L, 1L] == 2) 0.001 else 0.04
  1e5 * ((p_mis - target)^2 + (p_wgd - 0.02)^2)
}
fit_groups <- function(grouped) nk_select_interval_fits(nk_interval_fit_grid(grouped, c(0.001, 0.04), c(0, 0.02)))
fits <- fit_groups(groups)
stopifnot(length(unique(fits$pmis)) == 2L, all(fits$pwgd == 0.02))
ids <- nk_interval_roster(groups)$physical_interval_id
seeds <- forward_rng_roster(ids, 3L, 17L)
stopifnot(!anyDuplicated(seeds$seed), identical(seeds, forward_rng_roster(rev(ids), 3L, 17L)))
larger_seeds <- forward_rng_roster(c(ids, nk_physical_interval_id(fixture("physical_end_c"))), 4L, 17L)
for (i in seq_len(nrow(seeds))) {
  stopifnot(identical(seeds$seed[i], forward_rng_lookup(larger_seeds, seeds$physical_interval_id[i], seeds$stage[i], seeds$rep_idx[i])))
}
bad_seeds <- seeds
bad_seeds$seed[2L] <- bad_seeds$seed[1L]
expect_error(forward_rng_validate(bad_seeds), "collision")
bad_seeds <- seeds
bad_seeds$seed[1L] <- bad_seeds$seed[1L] + 1L
expect_error(forward_rng_validate(bad_seeds), "identity mismatch")
message("PASS: stable interval/stage/replicate seeds, collision checks and roster extension invariance.")

fake_calls <- 0L
fail_seed <- NA_integer_
simulate_interval_forward <- function(interval, p_mis, p_wgd, bottleneck_size,
                                      expansion_factor, record_every, divide_all_each_step, seed) {
  if (!is.na(fail_seed) && seed == fail_seed) stop("Injected synthetic task failure")
  fake_calls <<- fake_calls + 1L
  set.seed(seed)
  mat <- interval$start_karyotype[sample(seq_len(nrow(interval$start_karyotype)), bottleneck_size, replace = TRUE), , drop = FALSE]
  mat[1L, 1L] <- mat[1L, 1L] + sample(0:1, 1L)
  list(interval = interval, p_mis = p_mis, p_wgd = p_wgd, final_matrix = mat)
}
fake_runtime_tag <- "fake_task_test_only"
forward_checkpoint_runtime <- function() list(fixture_runtime = fake_runtime_tag, rng_kind = RNGkind())
run_fixture <- function(path, grouped = groups, selected = fit_groups(grouped), ...) {
  run_checkpointed_forward_runs(grouped[[1L]], selected, 4L, 2L, 3L, 1L, 17L,
                                 checkpoint_dir = path, fit_mode = "interval_shared_wgd_v1",
                                 seed_roster = seeds, ...)
}
without_receipts <- function(x) {
  attr(x, "checkpoint_receipts") <- NULL
  attr(x, "seed_roster") <- NULL
  x
}
baseline_path <- file.path(test_root, "raw")
baseline <- run_fixture(baseline_path)
stopifnot(fake_calls == 6L, all(!attr(baseline, "checkpoint_receipts")$reused))
resumed <- run_fixture(baseline_path, batch_size = 2L)
stopifnot(fake_calls == 6L, all(attr(resumed, "checkpoint_receipts")$reused),
          identical(without_receipts(resumed), without_receipts(baseline)))
reordered_groups <- list(group_a = rev(groups$group_a))
reordered <- run_fixture(baseline_path, reordered_groups, fits[2:1, , drop = FALSE])
stopifnot(fake_calls == 6L, all(attr(reordered, "checkpoint_receipts")$reused),
          identical(without_receipts(reordered), rev(without_receipts(baseline))))
message("PASS: exact interval parameters and raw task reuse survive fit/interval row reordering.")

# A neighbor's changed endpoint alters its selected p_mis and the coupled group
# profile. The untouched interval retains all raw simulation tasks.
changed_groups <- groups
changed_groups$group_a[[2L]]$end_karyotype[1L, 1L] <- 2
changed_fits <- fit_groups(changed_groups)
changed <- run_fixture(baseline_path, changed_groups, changed_fits)
stopifnot(sum(attr(changed, "checkpoint_receipts")$reused) == 3L, fake_calls == 9L)
# Changing observed data without changing selected parameters invalidates the
# profile/check contract, while actual forward simulator inputs remain equal.
profile_groups <- groups
profile_groups$group_a[[2L]]$end_karyotype[2L, 1L] <- 4
profile_fits <- fit_groups(profile_groups)
profile_only <- run_fixture(baseline_path, profile_groups, profile_fits)
stopifnot(all(attr(profile_only, "checkpoint_receipts")$reused), fake_calls == 9L,
          !identical(profile_fits$fit_id, fits$fit_id))
start_groups <- groups
start_groups$group_a[[1L]]$start_karyotype[1L, 1L] <- 4
start_changed <- run_fixture(baseline_path, start_groups)
stopifnot(sum(attr(start_changed, "checkpoint_receipts")$reused) == 3L, fake_calls == 12L)
message("PASS: input/parameter changes invalidate only affected raw simulations; profile-only changes reuse all.")
fake_runtime_tag <- "fake_changed_runtime"
runtime_changed <- run_fixture(baseline_path)
stopifnot(all(!attr(runtime_changed, "checkpoint_receipts")$reused))
fake_runtime_tag <- "fake_task_test_only"
runtime_restored <- run_fixture(baseline_path)
stopifnot(all(attr(runtime_restored, "checkpoint_receipts")$reused))
original_counts <- matrix_to_named_counts
matrix_to_named_counts <- function(M) original_counts(M)
kernel_changed <- run_fixture(baseline_path)
stopifnot(all(!attr(kernel_changed, "checkpoint_receipts")$reused))
matrix_to_named_counts <- original_counts
message("PASS: changed runtime or simulator kernel invalidates all dependent tasks without deleting previous contracts.")

partial_path <- file.path(test_root, "partial")
partial <- run_fixture(partial_path, interval_indices = 2L, replicate_indices = 1L, load_results = FALSE)
stopifnot(!partial$complete, nrow(partial$receipts) == 1L)
finished <- run_fixture(partial_path, load_results = FALSE)
stopifnot(finished$complete, sum(finished$receipts$reused) == 1L)
failure_path <- file.path(test_root, "failure")
fail_seed <- forward_rng_lookup(seeds, ids[1L], "forward_simulation", 2L)
expect_error(run_fixture(failure_path, batch_size = 3L), "Injected synthetic task failure")
stopifnot(nrow(forward_checkpoint_receipts(failure_path)) == 1L)
fail_seed <- NA_integer_
repaired <- run_fixture(failure_path, load_results = FALSE)
stopifnot(repaired$complete, sum(repaired$receipts$reused) == 1L)
task_path <- repaired$receipts$checkpoint_path[1L]
corrupt <- readRDS(task_path)
corrupt$result$final_matrix[1L, 1L] <- 99
saveRDS(corrupt, task_path)
expect_error(run_fixture(failure_path), "payload checksum mismatch")
message("PASS: interrupted tasks resume and corrupt completed payloads fail closed.")

summarize_fixture <- function(runs, grouped = groups, selected = fit_groups(grouped), n_pairs = 2L) {
  summarize_forward_validation(without_receipts(runs), grouped[[1L]], names(grouped)[1L],
                                 n_null_pairs = n_pairs, distance_metric = "wasserstein",
                                 fit_mode = "interval_shared_wgd_v1", fit_rows = selected,
                                 seed_roster = seeds, checkpoint_dir = baseline_path)
}
summary <- summarize_fixture(baseline)
summary_again <- summarize_fixture(baseline)
stopifnot(all(vapply(summary_again$interval_checks, function(x) isTRUE(attr(x, "checkpoint_reused")), logical(1))),
          identical(summary$interval_summary$p_value, summary$interval_summary$posterior_predictive_p),
          all(summary$interval_summary$check_method == "posterior_predictive_check_v1"))
reordered_summary <- summarize_fixture(reordered, reordered_groups, fits[2:1, , drop = FALSE])
by_id <- function(df) { out <- df[order(df$physical_interval_id), setdiff(names(df), "interval_idx"), drop = FALSE]; rownames(out) <- NULL; out }
stopifnot(identical(by_id(summary$interval_summary), by_id(reordered_summary$interval_summary)))
profile_summary <- summarize_fixture(profile_only, profile_groups, profile_fits)
stopifnot(!identical(summary$interval_summary$fit_id, profile_summary$interval_summary$fit_id))
changed_check <- summarize_fixture(baseline, n_pairs = 1L)
stopifnot(all(!vapply(changed_check$interval_checks, function(x) isTRUE(attr(x, "checkpoint_reused")), logical(1))))
extinct <- without_receipts(baseline)
extinct[[1L]][[2L]]$final_matrix <- matrix(numeric(), nrow = 0L, ncol = 2L)
extinct_summary <- summarize_fixture(extinct)
stopifnot(identical(extinct_summary$interval_checks[[1L]]$simulation_rep_indices, c(1L, 3L)),
          setequal(extinct_summary$ploidy$rep_idx[extinct_summary$ploidy$interval_idx == 1L & extinct_summary$ploidy$type == "simulated_matched"], c(1L, 3L)))
message("PASS: posterior-predictive checks cache independently and reorder deterministically.")

summary_dir <- file.path(test_root, "summaries")
dir.create(summary_dir)
saveRDS(groups, file.path(test_root, "groups.Rds"))
result <- list(forward_group_id = "group_a", fit_row = fits,
               settings = list(fit_mode = "interval_shared_wgd_v1"), forward_validation = summary, seed_roster = seeds)
summary_path <- file.path(summary_dir, "forward_validation_group_a_summary.Rds")
saveRDS(forward_validation_summary_object(result), summary_path)
csv <- write_neutral_probability_csv(summary_dir, file.path(test_root, "groups.Rds"), file.path(test_root, "summary.csv"), fit_mode = "interval_shared_wgd_v1")
idx <- match(csv$physical_interval_id, fits$physical_interval_id)
stopifnot(identical(csv$cin_rate, fits$pmis[idx]), identical(csv$fit_id, fits$fit_id[idx]),
          identical(csv$p_value, csv$posterior_predictive_p), all(csv$check_method == "posterior_predictive_check_v1"))
expect_error(write_neutral_probability_csv(summary_dir, file.path(test_root, "groups.Rds"), file.path(test_root, "legacy.csv")), "cannot consume")
expect_error(select_best_group_fits(fits), "cannot consume")
bad_summary <- forward_validation_summary_object(result)
bad_summary$interval_summary$pmis[1L] <- 0.5
saveRDS(bad_summary, summary_path)
expect_error(write_neutral_probability_csv(summary_dir, file.path(test_root, "groups.Rds"), file.path(test_root, "bad.csv"), fit_mode = "interval_shared_wgd_v1"), "numeric or identity mismatch")
bad_summary <- forward_validation_summary_object(result)
bad_summary$interval_summary$p_value[1L] <- 0.01
saveRDS(bad_summary, summary_path)
expect_error(write_neutral_probability_csv(summary_dir, file.path(test_root, "groups.Rds"), file.path(test_root, "bad.csv"), fit_mode = "interval_shared_wgd_v1"), "must preserve")
saveRDS(forward_validation_summary_object(result), summary_path)
saveRDS(forward_validation_summary_object(result), file.path(summary_dir, "forward_validation_duplicate_summary.Rds"))
expect_error(write_neutral_probability_csv(summary_dir, file.path(test_root, "groups.Rds"), file.path(test_root, "bad.csv"), fit_mode = "interval_shared_wgd_v1"), "exactly once")
message("PASS: interval CSV preserves original p-values and exact fit identities, rejecting mixed/stale rows.")

ploidy <- build_ploidy_diagnostics(groups, fits, fit_mode = "interval_shared_wgd_v1")
marginal <- build_marginal_cn_diagnostics(groups, fits, fit_mode = "interval_shared_wgd_v1")
stopifnot(setequal(unique(ploidy$predicted$physical_interval_id), ids),
          setequal(unique(marginal$predicted$physical_interval_id), ids),
          all(ploidy$predicted$pmis == fits$pmis[match(ploidy$predicted$physical_interval_id, fits$physical_interval_id)]),
          all(marginal$predicted$pmis == fits$pmis[match(marginal$predicted$physical_interval_id, fits$physical_interval_id)]))
message("PASS: native model diagnostics retain each interval's selected parameter pair.")

# Preserve all non-estimable intervals as unavailable rows rather than dropping
# them or simulating an arbitrary replacement pair.
get_interval_nll_markov <- function(p_mis, p_wgd, K0, KT, n_div, states = 0:8, eps = 1e-9) Inf
unavailable <- fit_groups(groups)
unavailable_summary <- summarize_fixture(list(list(), list()), selected = unavailable)
stopifnot(nrow(unavailable_summary$interval_summary) == 2L,
          all(unavailable_summary$interval_summary$forward_status == "not_simulated_non_estimable_fit"),
          all(is.na(unavailable_summary$interval_summary$posterior_predictive_p)))
message("PASS: non-estimable fits retain complete interval accounting without invented simulations.")

if ("--contracts-only" %in% commandArgs(trailingOnly = TRUE)) {
  unlink(test_root, recursive = TRUE)
  message("Focused synthetic algebra/dispatch/parity/resume checks passed; simulator and production runs were skipped.")
  quit(status = 0L)
}

# One real simulator smoke is small enough for an interactive test. Separate
# roots exercise worker-count invariance; every population is only four cells.
get_interval_nll_markov <- original_nll
simulate_interval_forward <- original_simulator
forward_checkpoint_runtime <- original_runtime
real_settings <- list(group_intervals = groups$group_a, fit_row = fits, bottleneck_size = 4L,
                      expansion_factor = 2L, n_reps = 2L, record_every = 1L, seed = 17L,
                      fit_mode = "interval_shared_wgd_v1", seed_roster = forward_rng_roster(ids, 2L, 17L))
real_serial <- do.call(run_checkpointed_forward_runs, c(real_settings, list(checkpoint_dir = file.path(test_root, "real_serial"))))
real_parallel <- do.call(run_checkpointed_forward_runs, c(real_settings, list(checkpoint_dir = file.path(test_root, "real_parallel"), n_cores = 2L, batch_size = 2L)))
stopifnot(identical(without_receipts(real_serial), without_receipts(real_parallel)))
real_resumed <- do.call(run_checkpointed_forward_runs, c(real_settings, list(checkpoint_dir = file.path(test_root, "real_parallel"), n_cores = 1L)))
stopifnot(all(attr(real_resumed, "checkpoint_receipts")$reused),
          identical(without_receipts(real_serial), without_receipts(real_resumed)))
message("PASS: real simulator serial/two-worker/resume equality with stable interval streams.")

# Thin entrypoints must publish the same contract and resume raw tasks when only
# diagnostic settings change. These are still tiny synthetic simulations.
saveRDS(fits, file.path(test_root, "fits.Rds"))
runner_path <- file.path(test_root, "runner.Rds")
runner_summary <- file.path(test_root, "runner_summary.Rds")
runner <- run_forward_validation_main(file.path(test_root, "groups.Rds"), file.path(test_root, "fits.Rds"),
                                      "group_a", runner_path, runner_summary, bottleneck_size = 4L,
                                      expansion_factor = 2L, n_reps = 2L, record_every = 1L, n_null_pairs = 1L,
                                      n_cores = 1L, seed = 17L, distance_metric = "wasserstein", fit_mode = "interval_shared_wgd_v1")
rerun <- load_forward_validation_result("group_a", file.path(test_root, "groups.Rds"), file.path(test_root, "fits.Rds"),
                                        runner_path, runner_summary, bottleneck_size = 4L, expansion_factor = 2L,
                                        n_reps = 2L, record_every = 1L, n_null_pairs = 2L, n_cores = 1L,
                                        seed = 17L, distance_metric = "wasserstein", fit_mode = "interval_shared_wgd_v1")
stopifnot(identical(runner$settings$fit_mode, "interval_shared_wgd_v1"),
          all(rerun$checkpoint_receipts$reused),
          identical(rerun$forward_validation$interval_summary$p_value, rerun$forward_validation$interval_summary$posterior_predictive_p))
message("PASS: runner and cache loader use strict interval mode and resume raw simulations across check-setting changes.")

all_groups <- c(groups, list(group_b = list(fixture("physical_end_c", 3))))
all_fits <- fit_groups(all_groups)
all_group_path <- file.path(test_root, "all_groups.Rds")
all_fit_path <- file.path(test_root, "all_fits.Rds")
saveRDS(all_groups, all_group_path)
saveRDS(all_fits, all_fit_path)
run_all <- function() run_forward_validation_all_main(
  all_group_path, all_fit_path, output_dir = file.path(test_root, "all_runs"),
  summary_output_dir = file.path(test_root, "all_summaries"), bottleneck_size = 4L,
  expansion_factor = 2L, n_reps = 2L, record_every = 1L, n_null_pairs = 1L,
  n_cores = 1L, seed = 17L, distance_metric = "wasserstein", fit_mode = "interval_shared_wgd_v1",
  aggregate_output_path = file.path(test_root, "all_diagnostics.csv"))
all_paths <- run_all()
first_all <- utils::read.csv(file.path(test_root, "all_diagnostics.csv"), stringsAsFactors = FALSE)
stopifnot(length(all_paths) == 2L, nrow(first_all) == 3L, identical(first_all$p_value, first_all$posterior_predictive_p))
saveRDS(list(group_b = all_groups$group_b, group_a = rev(all_groups$group_a)), all_group_path)
saveRDS(all_fits[rev(seq_len(nrow(all_fits))), , drop = FALSE], all_fit_path)
all_reordered_paths <- run_all()
stopifnot(all(vapply(all_reordered_paths, function(path) all(readRDS(path)$checkpoint_receipts$reused), logical(1))))
second_all <- utils::read.csv(file.path(test_root, "all_diagnostics.csv"), stringsAsFactors = FALSE)
stopifnot(identical(by_id(first_all), by_id(second_all)))
message("PASS: all-group entrypoint dispatches once per coupled group, publishes the exact roster and reuses tasks after group/interval reordering.")

unlink(test_root, recursive = TRUE)
message("All interval forward checkpoint smoke tests passed.")

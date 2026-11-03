source("R/estimate_pmis.R")
source("R/interval_fit_utils.R")
source("R/workflow_utils.R")
source("R/forward_checkpoint_utils.R")
source("R/interval_calibration.R")

same <- function(x, y) stopifnot(isTRUE(all.equal(x, y, tolerance = 1e-12, check.attributes = FALSE)))
fails <- function(expr, pattern = NULL) {
  error <- tryCatch(force(expr), error = identity)
  stopifnot(inherits(error, "error"))
  if (!is.null(pattern)) stopifnot(grepl(pattern, conditionMessage(error), fixed = TRUE))
}
fixture <- function(endpoint = c(2, 2, 4), nchr = 2L, divisions = 2L) {
  list(connected_tree = "calibration_fixture", start = "physical_start", end = "physical_end",
    start_ids = "sample_start", end_ids = "sample_end", start_alias = NULL, n_collapsed_passages = 0L,
    elapsed_passages = 1L, conditions = c(fixture = 1L),
    start_karyotype = matrix(2, nrow = 3L, ncol = nchr),
    end_karyotype = matrix(rep(endpoint, nchr), ncol = nchr), n_divisions = divisions)
}

# The manual statistic differs from a median of simulated-pair distances.
y <- matrix(c(4, 2), nrow = 1L)
bank <- list(matrix(c(2, 2), nrow = 1L), matrix(c(2, 2), nrow = 1L), matrix(c(8, 2), nrow = 1L))
manual <- nk_calibration_statistic(y, bank)
same(manual$distances, c(2, 2, 4))
same(manual$statistic, 2)
fails(nk_calibration_statistic(y, c(bank, list(matrix(numeric(), 0L, 2L)))))

test_contracts <- function() {
  scratch <- tempfile("interval-calibration-")
  dir.create(scratch)
  original_simulator <- simulate_interval_forward
  original_nll <- get_interval_nll_markov
  on.exit({
    assign("simulate_interval_forward", original_simulator, envir = globalenv())
    assign("get_interval_nll_markov", original_nll, envir = globalenv())
    unlink(scratch, recursive = TRUE)
  }, add = TRUE)
  assign("calibration_test_calls", integer(), envir = globalenv())
  assign("simulate_interval_forward", function(interval, p_mis, p_wgd, seed, ...) {
    stopifnot(is.null(interval$end_karyotype))
    calibration_test_calls <<- c(calibration_test_calls, seed)
    set.seed(seed)
    endpoint <- sample(c(2, 4), 6L, replace = TRUE, prob = c(1 - 8 * p_wgd, 8 * p_wgd))
    list(final_matrix = matrix(rep(endpoint, ncol(interval$start_karyotype)), ncol = ncol(interval$start_karyotype)))
  }, envir = globalenv())
  interval <- fixture()
  args <- list(interval = interval, pmis_grid = c(.001, .04), pwgd_grid = c(0, .05),
    checkpoint_dir = file.path(scratch, "resumed"), n_outer = 5L, reference_bank_size = 3L,
    bank_namespace = "test-bank-a", outer_namespace = "test-outer-a", master_seed = 19L,
    bottleneck_size = 4L, expansion_factor = 2L, record_every = 1L, group_id = "fixture_group")
  run <- function(extra = list(), base = args) do.call(nk_calibrate_physical_interval, utils::modifyList(base, extra, keep.null = TRUE))
  partial_bank <- run(list(max_reference_tasks = 1L, max_outer_tasks = 0L))
  stopifnot(partial_bank$summary$calibration_status == "pending", partial_bank$new_reference_tasks == 1L,
    partial_bank$summary$n_outer_pending == 5L, is.na(partial_bank$summary$p_value),
    sum(partial_bank$reference_receipts$status == "complete") == 1L,
    sum(partial_bank$reference_receipts$status == "pending") == 2L)
  partial_outer <- run(list(max_outer_tasks = 1L))
  stopifnot(partial_outer$tasks_processed == 1L, partial_outer$summary$n_outer_complete == 1L,
    partial_outer$summary$n_outer_pending == 4L, is.na(partial_outer$summary$p_value))
  next_outer <- run(list(max_outer_tasks = 1L))
  stopifnot(next_outer$tasks_processed == 1L, next_outer$summary$n_outer_complete == 2L, sum(next_outer$outer$reused) == 1L)
  resumed <- run()
  uninterrupted <- run(list(checkpoint_dir = file.path(scratch, "uninterrupted")))
  same(resumed$summary, uninterrupted$summary)
  same(resumed$outer[c("rep_idx", "status", "statistic", "pmis", "pwgd", "fit_id", "bank_id")],
       uninterrupted$outer[c("rep_idx", "status", "statistic", "pmis", "pwgd", "fit_id", "bank_id")])
  stopifnot(resumed$summary$calibration_status == "complete", resumed$summary$pilot_only,
    !resumed$summary$inference_eligible, resumed$summary$n_outer_complete == 5L,
    identical(resumed$summary$check_method, "interval_refit_parametric_bootstrap_v1"))
  b <- sum(resumed$outer$statistic >= resumed$observed$statistic)
  same(resumed$summary$p_value, (1 + b) / 6)
  same(c(resumed$summary$mc_tail_lower, resumed$summary$mc_tail_upper), stats::binom.test(b, 5L)$conf.int)
  stopifnot(sum(resumed$selected_pair_frequencies$n) == 5L, all(resumed$refit_diagnostics$n_refitted == 5L))
  before <- length(calibration_test_calls)
  cached <- run(list(max_outer_tasks = 0L, max_reference_tasks = 0L))
  stopifnot(length(calibration_test_calls) == before, all(cached$outer$reused), cached$new_reference_tasks == 0L)
  same(cached$summary, resumed$summary)

  # Independent external check of every pseudo refit and its fixed-bank median.
  for (i in seq_len(nrow(resumed$outer))) {
    row <- resumed$outer[i, ]
    ep <- readRDS(row$endpoint_path)$data$endpoint
    stopifnot(nrow(ep) == nrow(interval$end_karyotype))
    raw <- expand.grid(pmis = args$pmis_grid, pwgd = args$pwgd_grid)
    raw$nll <- mapply(function(p, q) get_interval_nll_markov(p, q, interval$start_karyotype, ep, interval$n_divisions), raw$pmis, raw$pwgd)
    best <- raw[order(raw$nll, raw$pwgd, raw$pmis), ][1L, ]
    same(c(row$pmis, row$pwgd), c(best$pmis, best$pwgd))
    statistic <- readRDS(row$statistic_path)$data
    same(row$statistic, stats::median(statistic$distances))
  }

  # Changing the observed endpoint cannot seed or identify comparison banks.
  prep <- function(value) nk_calibration_prepare(value, args$pmis_grid, args$pwgd_grid, args$n_outer,
    args$reference_bank_size, args$bank_namespace, args$outer_namespace, args$master_seed, 4L, 2L, 1L)
  context <- prep(interval)
  altered <- interval
  altered$end_karyotype[,] <- 3
  context_changed <- prep(altered)
  same(context$bank_contract, context_changed$bank_contract)
  same(context$seed_plan, context_changed$seed_plan)
  stopifnot(!"end_karyotype" %in% names(context$conditioning),
    !any(c("fit_id", "interval_grid_id", "semantic_input_sha256") %in% names(context$bank_contract)))
  old_outer_calls <- sum(calibration_test_calls >= 2 * 2^29)
  different_bank <- run(list(bank_namespace = "test-bank-b"))
  same(different_bank$seed_plan$roster$seed[different_bank$seed_plan$roster$stream %in% c("outer_simulation", "outer_sample")],
       resumed$seed_plan$roster$seed[resumed$seed_plan$roster$stream %in% c("outer_simulation", "outer_sample")])
  stopifnot(sum(calibration_test_calls >= 2 * 2^29) == old_outer_calls,
    different_bank$new_reference_tasks > 0L,
    different_bank$summary$bank_contract_sha256 != resumed$summary$bank_contract_sha256)
  same(different_bank$outer$endpoint_path, resumed$outer$endpoint_path)

  # More than 150,000 declared seeds: preflight without any simulator work.
  full_grid <- expand.grid(pmis = 10^seq(-8, -1, length.out = 19), pwgd = c(0, .001, .01, .05))
  large_seeds <- nk_calibration_seed_plan("physical_fixture", "conditioned_fixture", full_grid, 1000L, 1000L, "bank", "outer", 1L)
  stopifnot(nrow(large_seeds$roster) == 154000L, !anyDuplicated(large_seeds$roster$seed), large_seeds$collision_count == 0L)
  reorder <- nk_calibration_seed_plan("physical_fixture", "conditioned_fixture", full_grid[nrow(full_grid):1L, ], 1000L, 1000L, "bank", "outer", 1L)
  same(large_seeds, reorder)
  fails(nk_calibration_seed_plan("fixture", "fixture", full_grid, 2^29 + 1, 1L, "bank", "outer", 1L), "exceeds")

  # Corrupt payloads fail authentication; they are not silently recomputed.
  path <- resumed$reference_receipts$checkpoint_path[1L]
  saved <- readRDS(path)
  corrupt <- saved
  corrupt$data$endpoint[1L, 1L] <- 8
  saveRDS(corrupt, path, compress = "xz")
  fails(run(), "Stale or corrupt fit checkpoint")
  saveRDS(saved, path, compress = "xz")

  # An actual interrupt after saved bank tasks leaves them reusable on resume.
  assign("calibration_test_interrupt", TRUE, envir = globalenv())
  assign("simulate_interval_forward", function(interval, p_mis, p_wgd, seed, ...) {
    if (seed >= 2 * 2^29 && calibration_test_interrupt) {
      calibration_test_interrupt <<- FALSE
      stop(structure(list(message = "deliberate test interrupt"), class = c("interrupt", "condition")))
    }
    list(final_matrix = matrix(2, 4L, ncol(interval$start_karyotype)))
  }, envir = globalenv())
  interrupt_args <- utils::modifyList(args, list(checkpoint_dir = file.path(scratch, "interrupted"), n_outer = 2L, reference_bank_size = 2L))
  stopped <- tryCatch(run(base = interrupt_args), interrupt = function(e) "interrupted")
  stopifnot(identical(stopped, "interrupted"), length(list.files(file.path(interrupt_args$checkpoint_dir, "reference_endpoints"))) == 2L)
  recovered <- run(base = interrupt_args)
  fresh <- run(list(checkpoint_dir = file.path(scratch, "interrupt-fresh")), base = interrupt_args)
  same(recovered$summary, fresh$summary)
  stopifnot(all(recovered$reference_receipts$reused[recovered$reference_receipts$bank_id == recovered$observed$bank_id]))

  # Zero tails stay positive, and every failure remains in the requested B.
  zero <- run(list(interval = fixture(c(8, 8, 8)), checkpoint_dir = file.path(scratch, "zero-tail"), n_outer = 3L, reference_bank_size = 2L))
  stopifnot(zero$summary$n_tail_ge == 0L, zero$summary$p_value == 1 / 4,
    zero$summary$mc_tail_lower == 0, zero$summary$mc_tail_upper > 0)
  assign("simulate_interval_forward", function(interval, p_mis, p_wgd, seed, ...) {
    list(final_matrix = matrix(2, if (seed >= 2 * 2^29) 0L else 4L, ncol(interval$start_karyotype)))
  }, envir = globalenv())
  extinct <- run(list(checkpoint_dir = file.path(scratch, "extinct"), n_outer = 2L, reference_bank_size = 2L))
  stopifnot(extinct$summary$calibration_status == "extinct", extinct$summary$n_outer_extinct == 2L,
    extinct$summary$n_outer_pending == 0L, all(extinct$outer$failure_stage == "outer_simulation"), is.na(extinct$summary$p_value))
  assign("simulate_interval_forward", function(interval, p_mis, p_wgd, seed, ...) {
    list(final_matrix = matrix(if (seed >= 2 * 2^29) 9 else 2, 4L, ncol(interval$start_karyotype)))
  }, envir = globalenv())
  unsupported <- run(list(checkpoint_dir = file.path(scratch, "unsupported"), n_outer = 2L, reference_bank_size = 2L))
  stopifnot(unsupported$summary$calibration_status == "nonestimable", unsupported$summary$n_outer_nonestimable == 2L,
    all(grepl("outside_0_8", unsupported$outer$detail)), is.na(unsupported$summary$p_value))
  assign("calibration_original_nll", original_nll, envir = globalenv())
  assign("get_interval_nll_markov", function(p_mis, p_wgd, K0, KT, n_div, states = 0:8, eps = 1e-9) {
    if (all(KT == 5)) return(Inf)
    calibration_original_nll(p_mis, p_wgd, K0, KT, n_div, states, eps)
  }, envir = globalenv())
  assign("simulate_interval_forward", function(interval, p_mis, p_wgd, seed, ...) {
    list(final_matrix = matrix(if (seed >= 2 * 2^29) 5 else 2, 4L, ncol(interval$start_karyotype)))
  }, envir = globalenv())
  unavailable <- run(list(checkpoint_dir = file.path(scratch, "unavailable"), n_outer = 2L, reference_bank_size = 2L))
  stopifnot(unavailable$summary$calibration_status == "nonestimable", unavailable$summary$n_outer_nonestimable == 2L,
    all(unavailable$outer$failure_stage == "refit"), all(unavailable$outer$detail == "no_finite_interval_grid"), is.na(unavailable$summary$p_value))
  message("PASS: fixed-bank median, own-pair refits, seed universe, partial/interrupted resume, bank independence, corruption, zero tails and explicit failures.")
}
test_contracts()

# Tiny actual-kernel smoke; reduced resources are not scientific defaults.
smoke_root <- tempfile("interval-calibration-cpp-smoke-")
smoke <- nk_calibrate_physical_interval(fixture(c(2, 2), 22L, 1L), c(0, .01), c(0, .01), smoke_root,
  n_outer = 2L, reference_bank_size = 2L, bank_namespace = "cpp-smoke-bank", outer_namespace = "cpp-smoke-outer",
  bottleneck_size = 4L, expansion_factor = 2L, record_every = 1L)
stopifnot(smoke$summary$calibration_status == "complete", smoke$summary$p_value == 1,
  smoke$summary$pilot_only, !smoke$summary$inference_eligible, smoke$summary$n_endpoint_cells == 2L)
unlink(smoke_root, recursive = TRUE)
message("PASS: tiny actual simulator, native likelihood/refit and transport calibration smoke; no production eligibility inferred.")

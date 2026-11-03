# Conditional, pilot-only parametric-bootstrap checks of one exact native
# physical interval. Source the existing estimate/fit/workflow/checkpoint
# helpers first. The caller owns roster authorization and publication approval.

nk_calibration_method <- function() "interval_refit_parametric_bootstrap_v1"

nk_calibration_limit <- function(x, name) {
  if (identical(x, Inf)) return(x)
  nk_fit_integer(x, name, 0L)
}

nk_calibration_supported <- function(x) {
  is.matrix(x) && is.numeric(x) && nrow(x) > 0L && ncol(x) > 0L &&
    all(is.finite(x)) && all(round(x) %in% 0:8)
}

# Dedicated task allocation, not a replacement RNG. Odd strides permute each
# 2^29 block. The four blocks fit exactly in R's nonnegative signed-int seeds.
# Namespaces can change bank streams without changing outer streams.
nk_calibration_seed_plan <- function(physical_interval_id, conditioning_sha256,
                                     grid, n_outer, reference_bank_size,
                                     bank_namespace, outer_namespace, master_seed) {
  n_outer <- nk_fit_integer(n_outer, "n_outer", 1L)
  reference_bank_size <- nk_fit_integer(reference_bank_size, "reference_bank_size", 1L)
  master_seed <- nk_fit_integer(master_seed, "master_seed", 0L)
  bank_namespace <- nk_fit_ids(bank_namespace, "bank_namespace")
  outer_namespace <- nk_fit_ids(outer_namespace, "outer_namespace")
  grid <- grid[order(grid$pwgd, grid$pmis), c("pmis", "pwgd"), drop = FALSE]
  rownames(grid) <- NULL
  block_size <- 2^29
  n_reference <- nrow(grid) * as.double(reference_bank_size)
  if (!nrow(grid) || n_reference > block_size || n_outer > block_size ||
      2 * (n_reference + n_outer) > .Machine$integer.max) {
    stop("Calibration seed universe exceeds the declared four-block allocation.")
  }
  streams <- c("reference_simulation", "reference_sample", "outer_simulation", "outer_sample")
  allocations <- list()
  rows <- lapply(seq_along(streams), function(i) {
    reference <- i <= 2L
    namespace <- if (reference) bank_namespace else outer_namespace
    identity <- list(schema = "calibration_four_block_affine_v1", physical_interval_id = physical_interval_id,
      conditioning_sha256 = conditioning_sha256, namespace = namespace, master_seed = master_seed,
      stream = streams[i], grid = if (reference) grid else NULL,
      reference_bank_size = if (reference) reference_bank_size else NULL)
    hex <- nk_fit_hash(identity)
    # 28-bit chunks are exactly representable. Products remain below 2^49.
    offset <- strtoi(substr(hex, 1L, 7L), 16L)
    stride <- 2 * (strtoi(substr(hex, 8L, 14L), 16L) %% 2^19) + 1
    n <- if (reference) n_reference else n_outer
    rank <- seq_len(n) - 1
    seed <- as.integer((i - 1) * block_size + (offset + stride * rank) %% block_size)
    allocations[[i]] <<- data.frame(stream = streams[i], namespace = namespace,
      block_start = (i - 1) * block_size, block_size = block_size, offset = offset,
      stride = stride, n_tasks = n, stringsAsFactors = FALSE)
    data.frame(physical_interval_id = physical_interval_id, stream = streams[i], namespace = namespace,
      grid_index = if (reference) rep(seq_len(nrow(grid)), each = reference_bank_size) else NA_integer_,
      rep_idx = if (reference) rep(seq_len(reference_bank_size), nrow(grid)) else seq_len(n_outer),
      pmis = if (reference) rep(grid$pmis, each = reference_bank_size) else NA_real_,
      pwgd = if (reference) rep(grid$pwgd, each = reference_bank_size) else NA_real_,
      master_seed = master_seed, seed = seed, stringsAsFactors = FALSE)
  })
  roster <- do.call(rbind, rows)
  rownames(roster) <- NULL
  if (anyNA(roster$seed) || any(roster$seed < 0L) || anyDuplicated(roster$seed)) {
    stop("Calibration seed preflight found an overflow or collision; no tasks were launched.")
  }
  list(schema = "calibration_four_block_affine_v1", allocations = do.call(rbind, allocations),
       roster = roster, roster_sha256 = nk_fit_hash(roster), collision_count = 0L)
}

nk_calibration_prepare <- function(interval, pmis_grid, pwgd_grid, n_outer, reference_bank_size,
                                   bank_namespace, outer_namespace, master_seed,
                                   bottleneck_size, expansion_factor, record_every) {
  required <- c("nk_interval_fit_grid", "nk_select_interval_fits", "simulate_interval_forward",
                "sample_simulated_endpoint", "endpoint_wasserstein", "forward_checkpoint_runtime")
  missing <- required[!vapply(required, exists, logical(1), mode = "function", inherits = TRUE)]
  if (length(missing)) stop("Source existing estimate/fit/workflow/checkpoint helpers first: ", paste(missing, collapse = ", "))
  physical <- nk_interval_physical_identity(interval)
  input <- nk_interval_semantic_input(interval)
  for (x in list(pmis_grid, pwgd_grid)) {
    if (!is.numeric(x) || !length(x) || any(!is.finite(x)) || any(x < 0) || anyDuplicated(x)) {
      stop("Calibration requires complete, unique, finite nonnegative parameter grids.")
    }
  }
  pmis_grid <- sort(as.double(pmis_grid))
  pwgd_grid <- sort(as.double(pwgd_grid))
  if (8 * max(pmis_grid) + max(pwgd_grid) > 1) stop("Calibration grid violates the existing transition bounds.")
  n_outer <- nk_fit_integer(n_outer, "n_outer", 1L)
  reference_bank_size <- nk_fit_integer(reference_bank_size, "reference_bank_size", 1L)
  settings <- list(bottleneck_size = nk_fit_integer(bottleneck_size, "bottleneck_size", 1L),
    expansion_factor = nk_fit_integer(expansion_factor, "expansion_factor", 1L),
    record_every = nk_fit_integer(record_every, "record_every", 1L), divide_all_each_step = TRUE)
  if (as.double(settings$bottleneck_size) * settings$expansion_factor > .Machine$integer.max) {
    stop("Calibration simulator population size exceeds the existing integer limit.")
  }
  settings$cull_keep <- 1 / settings$expansion_factor
  model <- nk_interval_model_identity(0:8, 1e-9)
  fit_runtime <- nk_interval_runtime_identity()
  simulation_runtime <- forward_checkpoint_runtime()
  kernel <- list(cpp_sha256 = forward_checkpoint_file_hash("ksim2.cpp"),
    functions = forward_checkpoint_function_hash(c("simulate_interval_forward", "sample_initial_population",
      "matrix_to_named_counts", "named_counts_to_matrix", "ensure_ksim_compiled", "sample_simulated_endpoint",
      "nk_calibration_simulate_endpoint", "nk_calibration_supported")))
  statistic <- list(method = nk_calibration_method(), statistic = "median_matched_endpoint_wasserstein_v1",
    distance_metric = "wasserstein", tail = ">=", states = 0:8, eps = 1e-9,
    transport_version = as.character(utils::packageVersion("transport")),
    functions = forward_checkpoint_function_hash(c("endpoint_wasserstein", "nk_calibration_statistic")))
  # No observed endpoint or fit identity can enter this conditioning object.
  conditioning <- list(physical_interval_id = paste0("nki_", nk_fit_hash(physical)),
    start_karyotype = input$start_karyotype, n_divisions = input$n_divisions,
    n_endpoint_cells = nrow(input$end_karyotype), pmis_grid = pmis_grid, pwgd_grid = pwgd_grid,
    settings = settings, kernel_sha256 = nk_fit_hash(kernel), runtime_sha256 = nk_fit_hash(simulation_runtime))
  conditional_sha <- nk_fit_hash(conditioning)
  grid <- expand.grid(pmis = pmis_grid, pwgd = pwgd_grid, KEEP.OUT.ATTRS = FALSE)
  seed_plan <- nk_calibration_seed_plan(conditioning$physical_interval_id, conditional_sha, grid,
    n_outer, reference_bank_size, bank_namespace, outer_namespace, master_seed)
  bank_contract <- list(schema = "neutral_fixed_reference_bank_v1", conditioning_sha256 = conditional_sha,
    grid = grid, reference_bank_size = reference_bank_size, bank_namespace = bank_namespace,
    master_seed = master_seed, seed_schema = seed_plan$schema,
    seed_allocations = seed_plan$allocations[1:2, , drop = FALSE])
  state <- new.env(parent = emptyenv())
  state$banks <- list()
  state$bank_receipts <- list()
  state$reference_new <- 0L
  list(physical = physical, input = input, conditioning = conditioning, conditioning_sha256 = conditional_sha,
    grid = grid, model = model, fit_runtime = fit_runtime, simulation_runtime = simulation_runtime,
    kernel = kernel, statistic = statistic, seed_plan = seed_plan, bank_contract = bank_contract,
    bank_contract_sha256 = nk_fit_hash(bank_contract), n_outer = n_outer,
    reference_bank_size = reference_bank_size, state = state)
}

# Keep only the sufficient matched endpoint. No trajectory or outcome data is
# supplied to reference simulations. Extinction is retained, never resampled.
nk_calibration_simulate_endpoint <- function(context, pmis, pwgd, simulation_seed, sample_seed) {
  tryCatch({
    cnd <- context$conditioning
    simulation <- do.call(simulate_interval_forward, c(list(
      interval = list(start_karyotype = cnd$start_karyotype, n_divisions = cnd$n_divisions),
      p_mis = pmis, p_wgd = pwgd, seed = simulation_seed), cnd$settings))
    final <- simulation$final_matrix
    if (!is.matrix(final) || !is.numeric(final) || ncol(final) != ncol(cnd$start_karyotype)) {
      return(list(status = "failed", detail = "invalid_simulated_endpoint_shape", endpoint = NULL, n_final_cells = NA_integer_))
    }
    if (!nrow(final)) return(list(status = "extinct", detail = "empty_simulated_population", endpoint = final, n_final_cells = 0L))
    if (!nk_calibration_supported(final)) {
      return(list(status = "nonestimable", detail = "simulated_endpoint_outside_0_8_support", endpoint = NULL, n_final_cells = nrow(final)))
    }
    endpoint <- sample_simulated_endpoint(final, cnd$n_endpoint_cells, seed = sample_seed)
    if (nrow(endpoint) != cnd$n_endpoint_cells || !nk_calibration_supported(endpoint)) {
      return(list(status = "failed", detail = "invalid_matched_endpoint", endpoint = NULL, n_final_cells = nrow(final)))
    }
    list(status = "complete", detail = "", endpoint = endpoint, n_final_cells = nrow(final),
         endpoint_sha256 = nk_fit_hash(endpoint))
  }, error = function(e) list(status = "failed", detail = conditionMessage(e), endpoint = NULL, n_final_cells = NA_integer_))
}

nk_calibration_statistic <- function(endpoint, reference_endpoints) {
  if (!length(reference_endpoints) || !nk_calibration_supported(endpoint) ||
      any(!vapply(reference_endpoints, nk_calibration_supported, logical(1))) ||
      any(vapply(reference_endpoints, nrow, integer(1)) != nrow(endpoint))) {
    stop("The calibration median requires every declared matched endpoint; no dropping or clipping is allowed.")
  }
  distances <- vapply(reference_endpoints, function(x) endpoint_wasserstein(endpoint, x, distance_metric = "wasserstein"), numeric(1))
  if (any(!is.finite(distances))) stop("Nonfinite fixed-bank Wasserstein distance.")
  list(statistic = stats::median(distances), distances = distances)
}

nk_calibration_status <- function(status) {
  for (value in c("failed", "nonestimable", "extinct", "pending")) if (value %in% status) return(value)
  "complete"
}

nk_calibration_reference_bank <- function(context, pmis, pwgd, checkpoint_dir, max_reference_tasks) {
  pair <- which(context$grid$pmis == pmis & context$grid$pwgd == pwgd)
  if (length(pair) != 1L) stop("Selected reference-bank pair is not in the declared grid.")
  identity <- list(schema = "neutral_matched_reference_bank_v1", bank_contract_sha256 = context$bank_contract_sha256,
                   grid_index = pair, pmis = pmis, pwgd = pwgd)
  bank_id <- paste0("nkbank_", nk_fit_hash(identity))
  if (!is.null(context$state$banks[[bank_id]])) return(context$state$banks[[bank_id]])
  seeds <- context$seed_plan$roster
  sim <- seeds[seeds$stream == "reference_simulation" & !is.na(seeds$grid_index) & seeds$grid_index == pair, , drop = FALSE]
  sub <- seeds[seeds$stream == "reference_sample" & !is.na(seeds$grid_index) & seeds$grid_index == pair, , drop = FALSE]
  endpoints <- vector("list", context$reference_bank_size)
  receipts <- lapply(seq_len(context$reference_bank_size), function(i) {
    task_identity <- c(identity, list(rep_idx = i, simulation_seed = sim$seed[i], sample_seed = sub$seed[i]))
    sha <- nk_fit_hash(task_identity)
    path <- file.path(checkpoint_dir, "reference_endpoints", paste0(sha, ".Rds"))
    if (!file.exists(path) && context$state$reference_new >= max_reference_tasks) {
      data <- list(status = "pending", detail = "reference_task_budget", n_final_cells = NA_integer_)
      reused <- FALSE
      data_sha <- NA_character_
    } else {
      saved <- nk_fit_checkpoint(task_identity, function() {
        context$state$reference_new <- context$state$reference_new + 1L
        nk_calibration_simulate_endpoint(context, pmis, pwgd, sim$seed[i], sub$seed[i])
      }, checkpoint_dir, "reference_endpoints")
      data <- saved$data
      reused <- identical(saved$cache_status, "reused")
      data_sha <- nk_fit_hash(data)
      endpoints[i] <<- list(data$endpoint)
    }
    data.frame(bank_id = bank_id, rep_idx = i, pmis = pmis, pwgd = pwgd,
      simulation_seed = sim$seed[i], sample_seed = sub$seed[i], status = data$status, detail = data$detail,
      n_final_cells = data$n_final_cells, data_sha256 = data_sha, checkpoint_path = path, reused = reused,
      stringsAsFactors = FALSE)
  })
  receipts <- do.call(rbind, receipts)
  status <- nk_calibration_status(receipts$status)
  data_sha <- nk_fit_hash(receipts[c("rep_idx", "data_sha256", "status")])
  result <- list(bank_id = bank_id, status = status, endpoints = endpoints, receipts = receipts,
                 data_sha256 = data_sha, pmis = pmis, pwgd = pwgd)
  context$state$bank_receipts[[bank_id]] <- receipts
  if (!any(receipts$status == "pending")) context$state$banks[[bank_id]] <- result
  result
}

nk_calibration_fit <- function(interval, context, checkpoint_dir, group_id, reuse_grid = NULL) {
  if (!nk_calibration_supported(interval$start_karyotype) || !nk_calibration_supported(interval$end_karyotype)) {
    return(list(status = "nonestimable", detail = "fit_endpoint_outside_0_8_support", fit = NULL))
  }
  grouped <- stats::setNames(list(list(interval)), group_id)
  grid <- nk_interval_fit_grid(grouped, context$conditioning$pmis_grid, context$conditioning$pwgd_grid,
    checkpoint_dir = file.path(checkpoint_dir, "refits"), states = 0:8, eps = 1e-9,
    fit_mode = "interval_independent_v1", reuse_grid = reuse_grid)
  selected <- nk_select_interval_fits(grid)
  # Preserve the fit contract without duplicating large attributes in task rows.
  fit <- as.data.frame(stats::setNames(lapply(names(selected), function(nm) selected[[nm]]), names(selected)), stringsAsFactors = FALSE)
  list(status = if (is.finite(fit$negll) && is.finite(fit$pmis) && is.finite(fit$pwgd)) "complete" else "nonestimable",
       detail = if (is.finite(fit$negll)) "" else fit$estimability_status,
       fit = fit, grid_sha256 = attr(grid, "grid_sha256"),
       grid_checkpoints = attr(grid, "interval_checkpoints"), profile_checkpoints = attr(selected, "profile_checkpoints"))
}

nk_calibration_evaluate <- function(interval, context, checkpoint_dir, group_id, max_reference_tasks, reuse_grid = NULL) {
  fitted <- nk_calibration_fit(interval, context, checkpoint_dir, group_id, reuse_grid)
  out <- list(status = fitted$status, failure_stage = if (fitted$status == "complete") "" else "refit",
    detail = fitted$detail, fit = fitted$fit, fit_receipts = fitted[c("grid_checkpoints", "profile_checkpoints")],
    statistic = NA_real_, bank_id = NA_character_, bank_data_sha256 = NA_character_, statistic_path = NA_character_)
  if (out$status != "complete") return(out)
  bank <- nk_calibration_reference_bank(context, fitted$fit$pmis, fitted$fit$pwgd, checkpoint_dir, max_reference_tasks)
  out$bank_id <- bank$bank_id
  out$bank_data_sha256 <- bank$data_sha256
  if (bank$status != "complete") {
    out$status <- bank$status
    out$failure_stage <- "reference_bank"
    out$detail <- paste(unique(bank$receipts$detail[bank$receipts$status != "complete"]), collapse = "; ")
    return(out)
  }
  identity <- list(schema = "neutral_fixed_bank_statistic_v1", method = context$statistic,
    endpoint_sha256 = nk_fit_hash(nk_interval_semantic_input(interval)$end_karyotype),
    fit_id = fitted$fit$fit_id, bank_id = bank$bank_id, bank_data_sha256 = bank$data_sha256)
  saved <- nk_fit_checkpoint(identity, function() {
    tryCatch(c(list(status = "complete", detail = ""), nk_calibration_statistic(interval$end_karyotype, bank$endpoints)),
      error = function(e) list(status = "failed", detail = conditionMessage(e), statistic = NA_real_, distances = numeric()))
  }, checkpoint_dir, "statistics")
  out$status <- saved$data$status
  out$detail <- saved$data$detail
  out$failure_stage <- if (out$status == "complete") "" else "statistic"
  out$statistic <- saved$data$statistic
  out$distances <- saved$data$distances
  out$statistic_path <- saved$checkpoint_path
  out
}

nk_calibration_outer_row <- function(data, index, seeds, path = NA_character_, reused = FALSE) {
  out <- data.frame(rep_idx = index, status = data$status, failure_stage = data$failure_stage,
    detail = data$detail, statistic = data$statistic, bank_id = data$bank_id,
    simulation_seed = seeds[1L], sample_seed = seeds[2L], checkpoint_path = path, reused = reused,
    endpoint_path = if (is.null(data$endpoint_path)) NA_character_ else data$endpoint_path,
    statistic_path = if (is.null(data$statistic_path)) NA_character_ else data$statistic_path,
    stringsAsFactors = FALSE)
  numeric_fields <- c("pmis", "pwgd", "negll", "joint_tie_count")
  logical_fields <- c("pmis_at_lower_bound", "pmis_at_upper_bound", "pwgd_at_lower_bound", "pwgd_at_upper_bound")
  for (nm in c(numeric_fields, logical_fields, "fit_id", "fit_dependency_id", "estimability_status")) {
    out[[nm]] <- if (!is.null(data$fit)) data$fit[[nm]] else if (nm %in% numeric_fields) NA_real_ else
      if (nm %in% logical_fields) NA else NA_character_
  }
  out
}

# max_outer_tasks bounds unfinished outer tasks processed in this invocation;
# max_reference_tasks bounds newly simulated reference endpoints. Existing
# complete results are authenticated and reused without consuming these limits.
nk_calibrate_physical_interval <- function(interval, pmis_grid, pwgd_grid, checkpoint_dir,
                                           n_outer, reference_bank_size, bank_namespace, outer_namespace,
                                           master_seed = 1L, max_outer_tasks = Inf, max_reference_tasks = Inf,
                                           bottleneck_size = 2000L, expansion_factor = 32L, record_every = 25L,
                                           group_id = NULL, reuse_grid = NULL, mc_conf_level = 0.95) {
  checkpoint_dir <- nk_fit_ids(checkpoint_dir, "checkpoint_dir")
  max_outer_tasks <- nk_calibration_limit(max_outer_tasks, "max_outer_tasks")
  max_reference_tasks <- nk_calibration_limit(max_reference_tasks, "max_reference_tasks")
  if (length(mc_conf_level) != 1L || !is.numeric(mc_conf_level) || !is.finite(mc_conf_level) ||
      mc_conf_level <= 0 || mc_conf_level >= 1) stop("mc_conf_level must be between zero and one.")
  context <- nk_calibration_prepare(interval, pmis_grid, pwgd_grid, n_outer, reference_bank_size,
    bank_namespace, outer_namespace, master_seed, bottleneck_size, expansion_factor, record_every)
  metadata_group <- if (!is.null(group_id)) group_id else interval$group_id
  if (is.null(group_id)) group_id <- if (!is.null(interval$group_id)) interval$group_id else "__calibration_api_single_interval__"
  group_id <- nk_fit_ids(group_id, "group_id")
  if (!is.null(interval$group_id) && !identical(interval$group_id, group_id)) stop("Calibration group metadata conflicts with the native interval.")
  fit_method <- list(fit_mode = "interval_independent_v1", model = context$model, runtime = context$fit_runtime,
    conditional_fit_input_sha256 = nk_fit_hash(context$input[setdiff(names(context$input), "end_karyotype")]),
    functions = forward_checkpoint_function_hash(c("nk_calibration_fit", "nk_interval_fit_grid", "nk_select_independent_interval_fits",
      "nk_interval_reduce_independent", "nk_calibration_evaluate", "nk_calibration_supported", "nk_calibration_reference_bank")))
  identity <- list(schema = "neutral_interval_calibration_contract_v1", method = context$statistic,
    physical_identity = context$physical, observed_input_sha256 = nk_fit_hash(context$input),
    conditioning = context$conditioning, fit_method = fit_method, kernel = context$kernel,
    simulation_runtime = context$simulation_runtime, bank_contract = context$bank_contract,
    owner_functions = forward_checkpoint_function_hash(c("nk_calibrate_physical_interval", "nk_calibration_status", "nk_calibration_outer_row")),
    seed_plan_sha256 = context$seed_plan$roster_sha256, n_outer = context$n_outer,
    outer_namespace = outer_namespace, master_seed = master_seed, mc_conf_level = mc_conf_level,
    pilot_only = TRUE)
  contract <- nk_fit_checkpoint(identity, function() list(seed_plan = context$seed_plan,
    conditional_model = "observed_start_and_native_divisions", production_eligibility = "requires_separate_pilot_review"),
    checkpoint_dir, "contracts")
  observed <- nk_calibration_evaluate(interval, context, checkpoint_dir, group_id, max_reference_tasks, reuse_grid)
  seeds <- context$seed_plan$roster
  sim_seeds <- seeds$seed[seeds$stream == "outer_simulation"]
  sample_seeds <- seeds$seed[seeds$stream == "outer_sample"]
  work <- 0L
  pending <- function(reason) list(status = "pending", failure_stage = "", detail = reason,
    statistic = NA_real_, bank_id = NA_character_, fit = NULL)
  rows <- lapply(seq_len(context$n_outer), function(i) {
    seed_pair <- c(sim_seeds[i], sample_seeds[i])
    if (is.null(observed$fit) || !is.finite(observed$fit$pmis) || !is.finite(observed$fit$pwgd)) {
      return(nk_calibration_outer_row(pending("observed_fit_unavailable"), i, seed_pair))
    }
    endpoint_identity <- list(schema = "neutral_outer_matched_endpoint_v1", conditioning_sha256 = context$conditioning_sha256,
      generating_pmis = observed$fit$pmis, generating_pwgd = observed$fit$pwgd,
      outer_namespace = outer_namespace, master_seed = master_seed, seed_schema = context$seed_plan$schema,
      rep_idx = i, simulation_seed = seed_pair[1L], sample_seed = seed_pair[2L])
    result_identity <- list(schema = "neutral_outer_refit_statistic_v1", endpoint_identity = endpoint_identity,
      bank_contract_sha256 = context$bank_contract_sha256, fit_method_sha256 = nk_fit_hash(fit_method),
      statistic = context$statistic)
    result_path <- file.path(checkpoint_dir, "outer_results", paste0(nk_fit_hash(result_identity), ".Rds"))
    if (file.exists(result_path)) {
      data <- nk_fit_read_checkpoint(result_path, result_identity)
      if (!is.null(data$fit) && is.finite(data$fit$pmis) && is.finite(data$fit$pwgd)) {
        bank <- nk_calibration_reference_bank(context, data$fit$pmis, data$fit$pwgd, checkpoint_dir, 0L)
        if (data$status == "complete" && (bank$status != "complete" || !identical(bank$data_sha256, data$bank_data_sha256))) {
          stop("Saved outer statistic has missing or inconsistent reference-bank dependencies.")
        }
      }
      return(nk_calibration_outer_row(data, i, seed_pair, result_path, TRUE))
    }
    if (observed$status != "complete" || work >= max_outer_tasks) {
      return(nk_calibration_outer_row(pending(if (observed$status != "complete") "observed_statistic_unavailable" else "outer_task_budget"), i, seed_pair))
    }
    work <<- work + 1L
    endpoint <- nk_fit_checkpoint(endpoint_identity, function() nk_calibration_simulate_endpoint(
      context, observed$fit$pmis, observed$fit$pwgd, seed_pair[1L], seed_pair[2L]), checkpoint_dir, "outer_endpoints")
    if (endpoint$data$status != "complete") {
      data <- list(status = endpoint$data$status, failure_stage = "outer_simulation", detail = endpoint$data$detail,
                   statistic = NA_real_, bank_id = NA_character_, fit = NULL)
    } else {
      pseudo_interval <- interval
      pseudo_interval$end_karyotype <- endpoint$data$endpoint
      data <- nk_calibration_evaluate(pseudo_interval, context, checkpoint_dir, group_id, max_reference_tasks)
    }
    data$endpoint_path <- endpoint$checkpoint_path
    # Pending banks can resume; terminal failures are retained as final tasks.
    if (data$status != "pending") nk_fit_checkpoint(result_identity, function() data, checkpoint_dir, "outer_results")
    nk_calibration_outer_row(data, i, seed_pair, if (data$status == "pending") NA_character_ else result_path)
  })
  outer <- do.call(rbind, rows)
  status <- nk_calibration_status(c(observed$status, outer$status))
  complete <- identical(status, "complete")
  tail_count <- if (complete) sum(outer$statistic >= observed$statistic) else NA_integer_
  mc <- if (complete) stats::binom.test(tail_count, context$n_outer, conf.level = mc_conf_level)$conf.int else c(NA_real_, NA_real_)
  null_quantiles <- if (complete) stats::quantile(outer$statistic, c(.25, .5, .75), names = FALSE) else rep(NA_real_, 3L)
  bank_receipts <- if (length(context$state$bank_receipts)) do.call(rbind, context$state$bank_receipts) else data.frame()
  reference_status <- if (nrow(bank_receipts)) bank_receipts$status else character()
  summary <- data.frame(physical_interval_id = context$conditioning$physical_interval_id,
    group_id = if (is.null(metadata_group)) NA_character_ else metadata_group,
    fit_id = if (is.null(observed$fit)) NA_character_ else observed$fit$fit_id,
    fit_dependency_id = if (is.null(observed$fit)) NA_character_ else observed$fit$fit_dependency_id,
    fit_mode = "interval_independent_v1", fit_scope = "physical_interval", wgd_scope = "physical_interval",
    pmis = if (is.null(observed$fit)) NA_real_ else observed$fit$pmis,
    pwgd = if (is.null(observed$fit)) NA_real_ else observed$fit$pwgd,
    check_method = nk_calibration_method(), calibration_status = status, observed_status = observed$status,
    observed_failure_stage = observed$failure_stage, observed_detail = observed$detail,
    pilot_only = TRUE, inference_eligible = FALSE, distance_metric = "wasserstein",
    test_wasserstein = observed$statistic, null_wasserstein_q25 = null_quantiles[1L],
    null_wasserstein_median = null_quantiles[2L], null_wasserstein_q75 = null_quantiles[3L],
    distance_excess = observed$statistic - null_quantiles[2L],
    p_value = if (complete) (1 + tail_count) / (context$n_outer + 1) else NA_real_,
    n_tail_ge = tail_count, n_outer_requested = context$n_outer, n_outer_complete = sum(outer$status == "complete"),
    n_outer_pending = sum(outer$status == "pending"), n_outer_extinct = sum(outer$status == "extinct"),
    n_outer_nonestimable = sum(outer$status == "nonestimable"), n_outer_failed = sum(outer$status == "failed"),
    reference_bank_size = context$reference_bank_size, n_endpoint_cells = context$conditioning$n_endpoint_cells,
    n_reference_banks_used = length(context$state$bank_receipts),
    n_reference_complete = sum(reference_status == "complete"), n_reference_pending = sum(reference_status == "pending"),
    n_reference_extinct = sum(reference_status == "extinct"), n_reference_nonestimable = sum(reference_status == "nonestimable"),
    n_reference_failed = sum(reference_status == "failed"),
    n_divisions = context$conditioning$n_divisions, p_value_resolution = 1 / (context$n_outer + 1),
    mc_tail_lower = mc[1L], mc_tail_upper = mc[2L], mc_conf_level = mc_conf_level,
    mc_interval_method = "Clopper_Pearson_conditional_tail_probability",
    zero_tail_policy = "finite_MC_plus_one_never_zero_v1", tail_comparison = ">=",
    bank_namespace = bank_namespace, outer_namespace = outer_namespace, master_seed = master_seed,
    conditioning_sha256 = context$conditioning_sha256, bank_contract_sha256 = context$bank_contract_sha256,
    calibration_id = paste0("nkcal_", contract$sha256),
    method_limitation = "Approximate plug-in parametric bootstrap conditional on observed starts and fixed independent banks; finite-grid, boundary and simulator/Markov approximation remain; no joint-interval inference.",
    stringsAsFactors = FALSE)
  for (nm in c("passaging_ontology_run_id", "passaging_ontology_manifest_sha256")) {
    if (!is.null(interval[[nm]])) summary[[nm]] <- nk_fit_ids(interval[[nm]], nm)
  }
  refitted <- outer[is.finite(outer$pmis) & is.finite(outer$pwgd), , drop = FALSE]
  frequencies <- if (nrow(refitted)) stats::aggregate(list(n = refitted$rep_idx),
    by = list(pmis = refitted$pmis, pwgd = refitted$pwgd), FUN = length) else data.frame(pmis = numeric(), pwgd = numeric(), n = integer())
  frequencies$fraction_of_refitted <- frequencies$n / nrow(refitted)
  flag_names <- c("pmis_at_lower_bound", "pmis_at_upper_bound", "pwgd_at_lower_bound", "pwgd_at_upper_bound")
  diagnostics <- data.frame(diagnostic = c("joint_tie", flag_names), n_refitted = nrow(refitted),
    n = c(sum(refitted$joint_tie_count > 1L), vapply(flag_names, function(nm) sum(refitted[[nm]]), numeric(1))), stringsAsFactors = FALSE)
  diagnostics$fraction_of_refitted <- if (nrow(refitted)) diagnostics$n / nrow(refitted) else NA_real_
  result <- list(schema_version = "neutral_interval_calibration_pilot_v1", summary = summary, observed = observed,
    outer = outer, selected_pair_frequencies = frequencies, refit_diagnostics = diagnostics,
    reference_receipts = bank_receipts, seed_plan = context$seed_plan, contract_path = contract$checkpoint_path,
    tasks_processed = work, new_reference_tasks = context$state$reference_new)
  snapshot_identity <- list(schema = "neutral_calibration_snapshot_v1", contract_sha256 = contract$sha256,
    result_sha256 = nk_fit_hash(result))
  saved <- nk_fit_checkpoint(snapshot_identity, function() result, checkpoint_dir, "summaries")
  result$summary_path <- saved$checkpoint_path
  result
}

# Resumable storage for the existing interval/replicate simulator. Reporting
# settings deliberately do not invalidate scientifically identical simulations.
# The caller supplies authorized intervals; caching does not infer or authorize
# biological identities, histories, or ontology memberships.

forward_checkpoint_hash <- function(x) {
  if (!requireNamespace("digest", quietly = TRUE)) {
    stop("Forward checkpoints require the 'digest' package.")
  }
  digest::digest(x, algo = "sha256", serialize = TRUE, serializeVersion = 2L)
}

forward_checkpoint_file_hash <- function(path) {
  if (!file.exists(path)) stop("Missing checkpoint provenance file: ", path)
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}

forward_checkpoint_runtime <- function() {
  if (!requireNamespace("Rcpp", quietly = TRUE)) stop("Rcpp is required.")
  runtime_files <- c(
    file.path(R.home("etc"), "Makeconf"),
    file.path(R.home("lib"), paste0("libR", .Platform$dynlib.ext)),
    system.file("libs", paste0("Rcpp", .Platform$dynlib.ext), package = "Rcpp")
  )
  user_makevars <- Sys.getenv("R_MAKEVARS_USER", unset = path.expand("~/.R/Makevars"))
  runtime_files <- c(runtime_files, user_makevars[file.exists(user_makevars)])
  cxx <- system2(file.path(R.home("bin"), "R"), c("CMD", "config", "CXX"), stdout = TRUE)
  compiler <- Sys.which(strsplit(cxx[1L], "[[:space:]]+")[[1L]][1L])
  if (!nzchar(compiler)) stop("Could not resolve the configured simulator compiler.")
  list(
    R_version = R.version.string,
    platform = R.version$platform,
    Rcpp_version = as.character(utils::packageVersion("Rcpp")),
    compiler_sha256 = forward_checkpoint_file_hash(unname(compiler)),
    rng_kind = RNGkind(),
    runtime_files = stats::setNames(vapply(runtime_files, forward_checkpoint_file_hash, character(1)), runtime_files),
    build_environment = Sys.getenv(c("CC", "CXX", "CXX11", "CXX14", "CXX17", "CXXFLAGS", "PKG_CXXFLAGS", "PKG_CPPFLAGS", "PKG_LIBS", "R_MAKEVARS_USER", "R_MAKEVARS_SITE")),
    compiler_config = system2(file.path(R.home("bin"), "R"), c("CMD", "config", "--all"), stdout = TRUE, stderr = TRUE)
  )
}

forward_checkpoint_integer <- function(x, name, minimum = 1L) {
  if (length(x) != 1L || !is.numeric(x) || !is.finite(x) ||
      x < minimum || x > .Machine$integer.max || x != floor(x)) {
    stop(name, " must be an integer >= ", minimum, ".")
  }
  as.integer(x)
}

forward_checkpoint_is_interval <- function(contract) {
  identical(contract$identity$schema_version, "neutral_forward_checkpoint_v2")
}

# Seeds depend only on immutable interval identity, stage, replicate and the
# declared master seed. Collisions fail closed; no order-dependent probing.
forward_rng_seed <- function(physical_interval_id, stage, rep_idx, master_seed) {
  key <- list(version = "physical_interval_stage_rep_v1", master_seed = master_seed,
              physical_interval_id = physical_interval_id, stage = stage, rep_idx = rep_idx)
  hex <- forward_checkpoint_hash(key)
  value <- strtoi(substr(hex, 1L, 7L), 16L) * 16 + strtoi(substr(hex, 8L, 8L), 16L)
  as.integer(value %% (.Machine$integer.max - 1) + 1)
}

forward_rng_roster <- function(physical_interval_ids, n_reps, master_seed) {
  n_reps <- forward_checkpoint_integer(n_reps, "n_reps")
  master_seed <- forward_checkpoint_integer(master_seed, "master_seed", 0L)
  ids <- sort(as.character(physical_interval_ids))
  if (!length(ids) || anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) {
    stop("The RNG roster requires unique nonempty physical interval identities.")
  }
  out <- do.call(rbind, lapply(ids, function(id) {
    data.frame(physical_interval_id = id,
               stage = c(rep("forward_simulation", n_reps), rep("endpoint_subsample", n_reps), "endpoint_null_pairs"),
               rep_idx = c(seq_len(n_reps), seq_len(n_reps), 0L),
               master_seed = master_seed, stringsAsFactors = FALSE)
  }))
  out$seed <- vapply(seq_len(nrow(out)), function(i) {
    forward_rng_seed(out$physical_interval_id[i], out$stage[i], out$rep_idx[i], master_seed)
  }, integer(1))
  forward_rng_validate(out)
  rownames(out) <- NULL
  out
}

forward_rng_validate <- function(roster) {
  required <- c("physical_interval_id", "stage", "rep_idx", "master_seed", "seed")
  if (!is.data.frame(roster) || !nrow(roster) || !all(required %in% names(roster)) ||
      anyNA(roster[, required, drop = FALSE]) ||
      any(!nzchar(roster$physical_interval_id)) || any(!nzchar(roster$stage)) ||
      anyDuplicated(roster[, c("physical_interval_id", "stage", "rep_idx"), drop = FALSE])) {
    stop("Invalid or duplicate stable RNG roster keys.")
  }
  if (anyDuplicated(roster$seed)) {
    stop("Stable RNG seed collision; choose and record a different master seed before running.")
  }
  expected <- vapply(seq_len(nrow(roster)), function(i) {
    forward_rng_seed(as.character(roster$physical_interval_id[i]), as.character(roster$stage[i]),
                     forward_checkpoint_integer(roster$rep_idx[i], "rep_idx", 0L),
                     forward_checkpoint_integer(roster$master_seed[i], "master_seed", 0L))
  }, integer(1))
  if (!identical(as.integer(roster$seed), expected)) stop("Stable RNG roster seed identity mismatch.")
  invisible(roster)
}

forward_rng_lookup <- function(roster, physical_interval_id, stage, rep_indices) {
  rows <- roster[roster$physical_interval_id == physical_interval_id & roster$stage == stage, , drop = FALSE]
  idx <- match(rep_indices, rows$rep_idx)
  if (anyNA(idx) || anyDuplicated(rows$rep_idx)) stop("Missing or ambiguous stable RNG roster tasks.")
  rows$seed[idx]
}

forward_checkpoint_function_hash <- function(function_names) {
  stats::setNames(vapply(function_names, function(name) {
    fun <- get(name, mode = "function", inherits = TRUE)
    forward_checkpoint_hash(list(formals = formals(fun), body = body(fun)))
  }, character(1)), function_names)
}

forward_checkpoint_interval_contract <- function(group_intervals, fit_rows, settings, seed_roster = NULL) {
  if (!exists("nk_validate_interval_fits", mode = "function")) source("R/interval_fit_utils.R")
  if (!is.data.frame(fit_rows) || length(unique(fit_rows$group_id)) != 1L) {
    stop("Interval forward checkpoints require exactly one coupled group of interval fits.")
  }
  grouped <- stats::setNames(list(group_intervals), unique(as.character(fit_rows$group_id)))
  fit_rows <- nk_validate_interval_fits(fit_rows, grouped)
  ids <- unname(vapply(group_intervals, nk_physical_interval_id, character(1)))
  fit_rows <- fit_rows[match(ids, fit_rows$physical_interval_id), , drop = FALSE]
  if (is.null(seed_roster)) seed_roster <- forward_rng_roster(ids, settings$n_reps, settings$seed)
  forward_rng_validate(seed_roster)
  if (any(seed_roster$master_seed != settings$seed)) stop("RNG roster master seed differs from the requested seed.")
  for (id in ids) {
    forward_rng_lookup(seed_roster, id, "forward_simulation", seq_len(settings$n_reps))
    forward_rng_lookup(seed_roster, id, "endpoint_subsample", seq_len(settings$n_reps))
    forward_rng_lookup(seed_roster, id, "endpoint_null_pairs", 0L)
  }
  settings$seed_formula <- "physical_interval_stage_rep_v1"
  kernel <- list(
    cpp_sha256 = forward_checkpoint_file_hash("ksim2.cpp"),
    functions_sha256 = forward_checkpoint_function_hash(c("simulate_interval_forward", "sample_initial_population", "matrix_to_named_counts", "named_counts_to_matrix", "ensure_ksim_compiled")),
    task_schema = "neutral_raw_forward_task_v1"
  )
  runtime <- forward_checkpoint_runtime()
  # Store only simulator inputs in reusable payloads. Endpoint observations,
  # group profiles and labels belong to the invocation/check contracts below.
  simulation_intervals <- lapply(seq_along(ids), function(i) {
    list(physical_interval_id = ids[i],
         start_karyotype = nk_interval_semantic_input(group_intervals[[i]])$start_karyotype,
         n_divisions = as.integer(group_intervals[[i]]$n_divisions))
  })
  simulation_settings <- settings[c("bottleneck_size", "expansion_factor", "record_every", "divide_all_each_step")]
  raw_identities <- lapply(seq_along(ids), function(i) {
    list(schema_version = "neutral_raw_forward_input_v1", interval = simulation_intervals[[i]],
         pmis = unname(fit_rows$pmis[i]), pwgd = unname(fit_rows$pwgd[i]),
         settings = simulation_settings, kernel = kernel, runtime = runtime)
  })
  identity <- list(schema_version = "neutral_forward_checkpoint_v2", fit_mode = "interval_shared_wgd_v1",
                   group_intervals = group_intervals, fit_row = fit_rows, settings = settings,
                   physical_interval_ids = ids, simulation_intervals = simulation_intervals,
                   raw_identities = raw_identities, seed_roster = seed_roster,
                   runtime = runtime)
  list(identity = identity, identity_sha256 = forward_checkpoint_hash(identity))
}

forward_checkpoint_contract <- function(group_intervals, fit_row, bottleneck_size,
                                        expansion_factor, n_reps, record_every, seed,
                                        fit_mode = c("pooled_group_v1", "interval_shared_wgd_v1"),
                                        seed_roster = NULL) {
  fit_mode <- match.arg(fit_mode)
  if (!length(group_intervals)) stop("No intervals supplied for forward checkpoints.")
  if (!is.data.frame(fit_row) || (fit_mode == "pooled_group_v1" && nrow(fit_row) != 1L) ||
      !nrow(fit_row) ||
      !all(c("pmis", "pwgd") %in% names(fit_row)) ||
      any(!is.finite(c(fit_row$pmis, fit_row$pwgd)))) {
    stop("Forward checkpoints require complete fitted-parameter rows at the declared fit scope.")
  }
  if (fit_mode == "pooled_group_v1" && "physical_interval_id" %in% names(fit_row)) {
    stop("Interval fit rows cannot be read as legacy pooled group fits.")
  }
  for (interval in group_intervals) {
    for (field in c("start_karyotype", "end_karyotype")) {
      value <- interval[[field]]
      if (!is.matrix(value) || !is.numeric(value) || !nrow(value) ||
          !ncol(value) || any(!is.finite(value))) {
        stop("Each checkpoint interval requires a finite nonempty ", field, " matrix.")
      }
    }
    if (ncol(interval$start_karyotype) != ncol(interval$end_karyotype)) {
      stop("Checkpoint interval endpoint chromosome counts differ.")
    }
    forward_checkpoint_integer(interval$n_divisions, "n_divisions", 0L)
  }
  settings <- list(
    bottleneck_size = forward_checkpoint_integer(bottleneck_size, "bottleneck_size"),
    expansion_factor = forward_checkpoint_integer(expansion_factor, "expansion_factor"),
    n_reps = forward_checkpoint_integer(n_reps, "n_reps"),
    record_every = forward_checkpoint_integer(record_every, "record_every"),
    seed = forward_checkpoint_integer(seed, "seed", 0L),
    divide_all_each_step = TRUE,
    seed_formula = "seed + rep_idx + 1000 * interval_idx"
  )
  if (fit_mode == "interval_shared_wgd_v1") {
    return(forward_checkpoint_interval_contract(group_intervals, fit_row, settings, seed_roster))
  }
  if (as.double(seed) + n_reps + 1000 * length(group_intervals) > .Machine$integer.max) {
    stop("The existing deterministic task seed formula would overflow.")
  }
  source_files <- c("R/estimate_pmis.R", "ksim2.cpp", "R/forward_checkpoint_utils.R")
  identity <- list(
    schema_version = "neutral_forward_checkpoint_v1",
    group_intervals = group_intervals,
    fit_row = fit_row,
    settings = settings,
    simulation_source_sha256 = stats::setNames(vapply(source_files, forward_checkpoint_file_hash, character(1)), source_files),
    runtime = forward_checkpoint_runtime()
  )
  list(identity = identity, identity_sha256 = forward_checkpoint_hash(identity))
}

# Complete files are never overwritten. A failed serialization leaves no
# completed task, and a failed neighboring task does not discard saved work.
forward_checkpoint_atomic_save <- function(object, path, overwrite = FALSE) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(paste0(basename(path), ".partial-"), tmpdir = dirname(path))
  on.exit(unlink(temporary), add = TRUE)
  saveRDS(object, temporary, compress = "xz", version = 3L)
  if (file.exists(path) && !isTRUE(overwrite)) stop("Refusing to overwrite an existing checkpoint: ", path)
  if (!file.rename(temporary, path)) stop("Could not publish checkpoint atomically: ", path)
  invisible(path)
}

forward_checkpoint_read_contract <- function(checkpoint_dir, contract_sha256 = NULL) {
  path <- if (!is.null(contract_sha256)) file.path(checkpoint_dir, "contracts", paste0(contract_sha256, ".Rds")) else {
    if (file.exists(file.path(checkpoint_dir, "latest_contract.Rds"))) {
      file.path(checkpoint_dir, "latest_contract.Rds")
    } else file.path(checkpoint_dir, "contract.Rds")
  }
  contract <- tryCatch(readRDS(path), error = function(e) {
    stop("Unreadable forward checkpoint contract: ", path, "; ", conditionMessage(e))
  })
  if (!contract$identity$schema_version %in% c("neutral_forward_checkpoint_v1", "neutral_forward_checkpoint_v2") ||
      !identical(contract$identity_sha256, forward_checkpoint_hash(contract$identity))) {
    stop("Invalid forward checkpoint contract identity: ", path)
  }
  contract
}

forward_checkpoint_task_path <- function(checkpoint_dir, interval_idx, rep_idx, contract = NULL) {
  if (is.null(contract) && file.exists(file.path(checkpoint_dir, "latest_contract.Rds"))) {
    contract <- forward_checkpoint_read_contract(checkpoint_dir)
  }
  if (!is.null(contract) && forward_checkpoint_is_interval(contract)) {
    task <- forward_checkpoint_task_identity(contract, interval_idx, rep_idx)
    return(file.path(checkpoint_dir, "rawsim", paste0(task$rawsim_sha256, ".Rds")))
  }
  file.path(checkpoint_dir, sprintf("interval_%04d_rep_%05d.Rds", interval_idx, rep_idx))
}

forward_checkpoint_task_identity <- function(contract, interval_idx, rep_idx) {
  if (forward_checkpoint_is_interval(contract)) {
    id <- contract$identity$physical_interval_ids[interval_idx]
    identity <- list(physical_interval_id = id, stage = "forward_simulation", rep_idx = as.integer(rep_idx),
                     seed = forward_rng_lookup(contract$identity$seed_roster, id, "forward_simulation", rep_idx),
                     input_sha256 = forward_checkpoint_hash(contract$identity$raw_identities[[interval_idx]]))
    identity$rawsim_sha256 <- forward_checkpoint_hash(identity)
    return(identity)
  }
  list(
    contract_sha256 = contract$identity_sha256,
    interval_idx = as.integer(interval_idx),
    rep_idx = as.integer(rep_idx),
    seed = as.integer(contract$identity$settings$seed + rep_idx + 1000L * interval_idx)
  )
}

forward_checkpoint_read_task <- function(path, contract, interval_idx, rep_idx) {
  task <- tryCatch(readRDS(path), error = function(e) {
    stop("Unreadable forward checkpoint; do not silently reuse or replace it: ", path,
         "; ", conditionMessage(e))
  })
  expected <- forward_checkpoint_task_identity(contract, interval_idx, rep_idx)
  if (!identical(task$task_identity, expected) ||
      !identical(task$result_sha256, forward_checkpoint_hash(task$result))) {
    stop("Forward checkpoint task identity or payload checksum mismatch: ", path)
  }
  expected_interval <- if (forward_checkpoint_is_interval(contract)) contract$identity$simulation_intervals[[interval_idx]] else contract$identity$group_intervals[[interval_idx]]
  expected_fit <- if (forward_checkpoint_is_interval(contract)) contract$identity$fit_row[interval_idx, , drop = FALSE] else contract$identity$fit_row
  if (!identical(task$result$interval, expected_interval) ||
      !identical(task$result$p_mis, expected_fit$pmis) ||
      !identical(task$result$p_wgd, expected_fit$pwgd)) {
    stop("Forward checkpoint result does not match its scientific inputs: ", path)
  }
  task
}

forward_checkpoint_task_receipt <- function(path, task, reused, contract = NULL, interval_idx = NULL) {
  if (!is.null(contract) && forward_checkpoint_is_interval(contract)) {
    return(data.frame(physical_interval_id = task$task_identity$physical_interval_id,
                      interval_idx = as.integer(interval_idx), rep_idx = task$task_identity$rep_idx,
                      seed = task$task_identity$seed, rawsim_sha256 = task$task_identity$rawsim_sha256,
                      contract_sha256 = contract$identity_sha256,
                      result_sha256 = task$result_sha256, checkpoint_sha256 = forward_checkpoint_file_hash(path),
                      checkpoint_path = normalizePath(path, mustWork = TRUE), size_bytes = unname(file.info(path)$size),
                      simulation_elapsed_seconds = task$simulation_elapsed_seconds, reused = reused,
                      stringsAsFactors = FALSE))
  }
  data.frame(
    interval_idx = task$task_identity$interval_idx,
    rep_idx = task$task_identity$rep_idx,
    seed = task$task_identity$seed,
    contract_sha256 = task$task_identity$contract_sha256,
    result_sha256 = task$result_sha256,
    checkpoint_sha256 = forward_checkpoint_file_hash(path),
    checkpoint_path = normalizePath(path, mustWork = TRUE),
    size_bytes = unname(file.info(path)$size),
    simulation_elapsed_seconds = task$simulation_elapsed_seconds,
    reused = reused,
    stringsAsFactors = FALSE
  )
}

forward_checkpoint_tasks <- function(contract) {
  # Interval-major order keeps each native per-interval replicate roster intact.
  data.frame(
    interval_idx = rep(seq_along(contract$identity$group_intervals), each = contract$identity$settings$n_reps),
    rep_idx = rep(seq_len(contract$identity$settings$n_reps), times = length(contract$identity$group_intervals))
  )
}

forward_checkpoint_receipts <- function(checkpoint_dir, contract = forward_checkpoint_read_contract(checkpoint_dir)) {
  tasks <- forward_checkpoint_tasks(contract)
  receipts <- lapply(seq_len(nrow(tasks)), function(i) {
    path <- forward_checkpoint_task_path(checkpoint_dir, tasks$interval_idx[i], tasks$rep_idx[i], contract)
    if (!file.exists(path)) return(NULL)
    task <- forward_checkpoint_read_task(path, contract, tasks$interval_idx[i], tasks$rep_idx[i])
    forward_checkpoint_task_receipt(path, task, reused = TRUE, contract = contract, interval_idx = tasks$interval_idx[i])
  })
  if (!any(lengths(receipts))) return(data.frame())
  do.call(rbind, receipts)
}

forward_checkpoint_indices <- function(indices, universe, name) {
  if (is.null(indices)) return(universe)
  if (!is.numeric(indices) || !length(indices) || any(!is.finite(indices)) ||
      any(indices != as.integer(indices)) || any(!indices %in% universe) || anyDuplicated(indices)) {
    stop("Invalid checkpoint ", name, ".")
  }
  as.integer(indices)
}

# Subsetting retains the native indices in receipts. Callers of the existing
# summary owner must preserve its original interval-specific seed offset.
load_forward_checkpoint_runs <- function(checkpoint_dir, interval_indices = NULL,
                                        replicate_indices = NULL,
                                        contract = forward_checkpoint_read_contract(checkpoint_dir)) {
  interval_indices <- forward_checkpoint_indices(interval_indices, seq_along(contract$identity$group_intervals), "interval_indices")
  replicate_indices <- forward_checkpoint_indices(replicate_indices, seq_len(contract$identity$settings$n_reps), "replicate_indices")
  lapply(interval_indices, function(interval_idx) {
    lapply(replicate_indices, function(rep_idx) {
      path <- forward_checkpoint_task_path(checkpoint_dir, interval_idx, rep_idx, contract)
      result <- forward_checkpoint_read_task(path, contract, interval_idx, rep_idx)$result
      if (forward_checkpoint_is_interval(contract)) result$interval <- contract$identity$group_intervals[[interval_idx]]
      result
    })
  })
}

# Each worker saves its result before returning a small receipt. No completed
# simulation bank must remain in memory when load_results is FALSE.
run_checkpointed_forward_runs <- function(group_intervals, fit_row, bottleneck_size,
                                         expansion_factor, n_reps, record_every, seed,
                                         checkpoint_dir, n_cores = 1L,
                                         batch_size = n_cores, load_results = TRUE,
                                         interval_indices = NULL, replicate_indices = NULL,
                                         fit_mode = c("pooled_group_v1", "interval_shared_wgd_v1"),
                                         seed_roster = NULL) {
  fit_mode <- match.arg(fit_mode)
  n_cores <- forward_checkpoint_integer(n_cores, "n_cores")
  batch_size <- forward_checkpoint_integer(batch_size, "batch_size")
  if (length(checkpoint_dir) != 1L || is.na(checkpoint_dir) || !nzchar(checkpoint_dir)) {
    stop("A nonempty checkpoint_dir is required.")
  }
  contract <- forward_checkpoint_contract(group_intervals, fit_row, bottleneck_size,
                                          expansion_factor, n_reps, record_every, seed,
                                          fit_mode = fit_mode, seed_roster = seed_roster)
  dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
  checkpoint_dir <- normalizePath(checkpoint_dir, mustWork = TRUE)
  contract_path <- file.path(checkpoint_dir, "contract.Rds")
  if (forward_checkpoint_is_interval(contract)) {
    if (file.exists(contract_path)) stop("Legacy pooled and interval checkpoint contracts cannot share a directory.")
    snapshot_path <- file.path(checkpoint_dir, "contracts", paste0(contract$identity_sha256, ".Rds"))
    if (file.exists(snapshot_path)) {
      if (!identical(forward_checkpoint_read_contract(checkpoint_dir, contract$identity_sha256), contract)) {
        stop("Interval checkpoint contract snapshot mismatch.")
      }
    } else forward_checkpoint_atomic_save(contract, snapshot_path)
    # This convenience pointer is replaceable; content-addressed snapshots and
    # task files remain immutable. Every running invocation holds its contract.
    forward_checkpoint_atomic_save(contract, file.path(checkpoint_dir, "latest_contract.Rds"), overwrite = TRUE)
  } else if (file.exists(file.path(checkpoint_dir, "latest_contract.Rds"))) {
    stop("Legacy pooled and interval checkpoint contracts cannot share a directory.")
  } else if (file.exists(contract_path)) {
    existing <- forward_checkpoint_read_contract(checkpoint_dir)
    if (!identical(existing, contract)) {
      stop("Checkpoint scientific-input/model/settings/seed/runtime identity mismatch. Use a new checkpoint directory.")
    }
  } else {
    if (length(list.files(checkpoint_dir, pattern = "^interval_.*[.]Rds$"))) {
      stop("Checkpoint tasks exist without a contract; refusing unverified reuse.")
    }
    forward_checkpoint_atomic_save(contract, contract_path)
  }

  tasks <- forward_checkpoint_tasks(contract)
  interval_indices <- forward_checkpoint_indices(interval_indices, seq_along(group_intervals), "interval_indices")
  replicate_indices <- forward_checkpoint_indices(replicate_indices, seq_len(n_reps), "replicate_indices")
  requested <- tasks$interval_idx %in% interval_indices & tasks$rep_idx %in% replicate_indices
  expected_paths <- vapply(seq_len(nrow(tasks)), function(i) {
    forward_checkpoint_task_path(checkpoint_dir, tasks$interval_idx[i], tasks$rep_idx[i], contract)
  }, character(1))
  extra <- if (forward_checkpoint_is_interval(contract)) character() else setdiff(list.files(checkpoint_dir, pattern = "^interval_.*[.]Rds$", full.names = TRUE), expected_paths)
  if (length(extra)) stop("Unexpected checkpoint task files: ", paste(basename(extra), collapse = ", "))
  receipts <- vector("list", nrow(tasks))
  completed <- file.exists(expected_paths)
  for (i in which(completed)) {
    task <- forward_checkpoint_read_task(expected_paths[i], contract, tasks$interval_idx[i], tasks$rep_idx[i])
    receipts[[i]] <- forward_checkpoint_task_receipt(expected_paths[i], task, reused = TRUE, contract = contract, interval_idx = tasks$interval_idx[i])
  }
  message("Forward checkpoints: reused ", sum(completed), "/", nrow(tasks), " validated tasks.")
  pending <- which(requested & !completed)

  worker <- function(i) {
    interval_idx <- tasks$interval_idx[i]
    rep_idx <- tasks$rep_idx[i]
    task_identity <- forward_checkpoint_task_identity(contract, interval_idx, rep_idx)
    settings <- contract$identity$settings
    interval <- if (forward_checkpoint_is_interval(contract)) contract$identity$simulation_intervals[[interval_idx]] else contract$identity$group_intervals[[interval_idx]]
    fit <- if (forward_checkpoint_is_interval(contract)) contract$identity$fit_row[interval_idx, , drop = FALSE] else contract$identity$fit_row
    started <- proc.time()[["elapsed"]]
    result <- simulate_interval_forward(
      interval = interval,
      p_mis = fit$pmis,
      p_wgd = fit$pwgd,
      bottleneck_size = settings$bottleneck_size,
      expansion_factor = settings$expansion_factor,
      record_every = settings$record_every,
      divide_all_each_step = settings$divide_all_each_step,
      seed = task_identity$seed
    )
    task <- list(task_identity = task_identity, result_sha256 = forward_checkpoint_hash(result),
                 simulation_elapsed_seconds = unname(proc.time()[["elapsed"]] - started), result = result)
    forward_checkpoint_atomic_save(task, expected_paths[i])
    forward_checkpoint_task_receipt(expected_paths[i], task, reused = FALSE, contract = contract, interval_idx = interval_idx)
  }

  if (length(pending)) {
    cl <- NULL
    if (n_cores > 1L) {
      cl <- parallel::makeCluster(min(n_cores, batch_size, length(pending)))
      on.exit(parallel::stopCluster(cl), add = TRUE)
      parallel::clusterCall(cl, setwd, getwd())
      parallel::clusterCall(cl, function(paths) .libPaths(paths), .libPaths())
      parallel::clusterCall(cl, function(kinds) do.call(RNGkind, as.list(kinds)), RNGkind())
      parallel::clusterEvalQ(cl, {
        source("R/estimate_pmis.R")
        source("R/forward_checkpoint_utils.R")
        NULL
      })
      parallel::clusterExport(cl, c("contract", "tasks", "expected_paths"), envir = environment())
      runtimes <- parallel::clusterCall(cl, forward_checkpoint_runtime)
      if (any(!vapply(runtimes, identical, logical(1), contract$identity$runtime))) {
        stop("Checkpoint worker runtime differs from the recorded simulation runtime.")
      }
    }
    batches <- split(pending, ceiling(seq_along(pending) / batch_size))
    for (batch in batches) {
      new_receipts <- if (is.null(cl)) lapply(batch, worker) else parallel::parLapply(cl, batch, worker)
      receipts[batch] <- new_receipts
      completed[batch] <- TRUE
      message("Forward checkpoints: ", sum(completed), "/", nrow(tasks), " tasks complete.")
    }
  }
  receipt_table <- do.call(rbind, receipts)
  rownames(receipt_table) <- NULL
  if (!load_results) {
    return(list(checkpoint_dir = checkpoint_dir, contract = contract, receipts = receipt_table,
                requested_tasks = tasks[requested, , drop = FALSE], complete = all(completed)))
  }
  out <- load_forward_checkpoint_runs(checkpoint_dir, interval_indices, replicate_indices, contract = contract)
  attr(out, "checkpoint_receipts") <- receipt_table
  if (forward_checkpoint_is_interval(contract)) attr(out, "seed_roster") <- contract$identity$seed_roster
  out
}

# The current endpoint check is a plug-in diagnostic, not a calibrated test.
# Cache it separately from raw simulations and the still-pending refit method.
forward_checkpoint_endpoint_check <- function(identity, checkpoint_dir, compute) {
  key <- forward_checkpoint_hash(identity)
  path <- file.path(checkpoint_dir, "endpoint_checks", paste0(key, ".Rds"))
  if (file.exists(path)) {
    saved <- tryCatch(readRDS(path), error = function(e) stop("Unreadable endpoint-check checkpoint: ", path))
    if (!identical(saved$identity, identity) || !identical(saved$result_sha256, forward_checkpoint_hash(saved$result))) {
      stop("Endpoint-check checkpoint identity or payload checksum mismatch: ", path)
    }
    out <- saved$result
    attr(out, "checkpoint_reused") <- TRUE
    return(out)
  }
  out <- compute()
  forward_checkpoint_atomic_save(list(identity = identity, result_sha256 = forward_checkpoint_hash(out), result = out), path)
  attr(out, "checkpoint_reused") <- FALSE
  out
}

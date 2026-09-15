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

forward_checkpoint_contract <- function(group_intervals, fit_row, bottleneck_size,
                                        expansion_factor, n_reps, record_every, seed) {
  if (!length(group_intervals)) stop("No intervals supplied for forward checkpoints.")
  if (!is.data.frame(fit_row) || nrow(fit_row) != 1L ||
      !all(c("pmis", "pwgd") %in% names(fit_row)) ||
      any(!is.finite(c(fit_row$pmis, fit_row$pwgd)))) {
    stop("Forward checkpoints require one complete fitted-parameter row.")
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
forward_checkpoint_atomic_save <- function(object, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(paste0(basename(path), ".partial-"), tmpdir = dirname(path))
  on.exit(unlink(temporary), add = TRUE)
  saveRDS(object, temporary, compress = "xz", version = 3L)
  if (file.exists(path)) stop("Refusing to overwrite an existing checkpoint: ", path)
  if (!file.rename(temporary, path)) stop("Could not publish checkpoint atomically: ", path)
  invisible(path)
}

forward_checkpoint_read_contract <- function(checkpoint_dir) {
  path <- file.path(checkpoint_dir, "contract.Rds")
  contract <- tryCatch(readRDS(path), error = function(e) {
    stop("Unreadable forward checkpoint contract: ", path, "; ", conditionMessage(e))
  })
  if (!identical(contract$identity$schema_version, "neutral_forward_checkpoint_v1") ||
      !identical(contract$identity_sha256, forward_checkpoint_hash(contract$identity))) {
    stop("Invalid forward checkpoint contract identity: ", path)
  }
  contract
}

forward_checkpoint_task_path <- function(checkpoint_dir, interval_idx, rep_idx) {
  file.path(checkpoint_dir, sprintf("interval_%04d_rep_%05d.Rds", interval_idx, rep_idx))
}

forward_checkpoint_task_identity <- function(contract, interval_idx, rep_idx) {
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
  if (!identical(task$result$interval, contract$identity$group_intervals[[interval_idx]]) ||
      !identical(task$result$p_mis, contract$identity$fit_row$pmis) ||
      !identical(task$result$p_wgd, contract$identity$fit_row$pwgd)) {
    stop("Forward checkpoint result does not match its scientific inputs: ", path)
  }
  task
}

forward_checkpoint_task_receipt <- function(path, task, reused) {
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

forward_checkpoint_receipts <- function(checkpoint_dir) {
  contract <- forward_checkpoint_read_contract(checkpoint_dir)
  tasks <- forward_checkpoint_tasks(contract)
  receipts <- lapply(seq_len(nrow(tasks)), function(i) {
    path <- forward_checkpoint_task_path(checkpoint_dir, tasks$interval_idx[i], tasks$rep_idx[i])
    if (!file.exists(path)) return(NULL)
    task <- forward_checkpoint_read_task(path, contract, tasks$interval_idx[i], tasks$rep_idx[i])
    forward_checkpoint_task_receipt(path, task, reused = TRUE)
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
                                        replicate_indices = NULL) {
  contract <- forward_checkpoint_read_contract(checkpoint_dir)
  interval_indices <- forward_checkpoint_indices(interval_indices, seq_along(contract$identity$group_intervals), "interval_indices")
  replicate_indices <- forward_checkpoint_indices(replicate_indices, seq_len(contract$identity$settings$n_reps), "replicate_indices")
  lapply(interval_indices, function(interval_idx) {
    lapply(replicate_indices, function(rep_idx) {
      path <- forward_checkpoint_task_path(checkpoint_dir, interval_idx, rep_idx)
      forward_checkpoint_read_task(path, contract, interval_idx, rep_idx)$result
    })
  })
}

# Each worker saves its result before returning a small receipt. No completed
# simulation bank must remain in memory when load_results is FALSE.
run_checkpointed_forward_runs <- function(group_intervals, fit_row, bottleneck_size,
                                         expansion_factor, n_reps, record_every, seed,
                                         checkpoint_dir, n_cores = 1L,
                                         batch_size = n_cores, load_results = TRUE,
                                         interval_indices = NULL, replicate_indices = NULL) {
  n_cores <- forward_checkpoint_integer(n_cores, "n_cores")
  batch_size <- forward_checkpoint_integer(batch_size, "batch_size")
  if (length(checkpoint_dir) != 1L || is.na(checkpoint_dir) || !nzchar(checkpoint_dir)) {
    stop("A nonempty checkpoint_dir is required.")
  }
  contract <- forward_checkpoint_contract(group_intervals, fit_row, bottleneck_size,
                                          expansion_factor, n_reps, record_every, seed)
  dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
  checkpoint_dir <- normalizePath(checkpoint_dir, mustWork = TRUE)
  contract_path <- file.path(checkpoint_dir, "contract.Rds")
  if (file.exists(contract_path)) {
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
    forward_checkpoint_task_path(checkpoint_dir, tasks$interval_idx[i], tasks$rep_idx[i])
  }, character(1))
  extra <- setdiff(list.files(checkpoint_dir, pattern = "^interval_.*[.]Rds$", full.names = TRUE), expected_paths)
  if (length(extra)) stop("Unexpected checkpoint task files: ", paste(basename(extra), collapse = ", "))
  receipts <- vector("list", nrow(tasks))
  completed <- file.exists(expected_paths)
  for (i in which(completed)) {
    task <- forward_checkpoint_read_task(expected_paths[i], contract, tasks$interval_idx[i], tasks$rep_idx[i])
    receipts[[i]] <- forward_checkpoint_task_receipt(expected_paths[i], task, reused = TRUE)
  }
  message("Forward checkpoints: reused ", sum(completed), "/", nrow(tasks), " validated tasks.")
  pending <- which(requested & !completed)

  worker <- function(i) {
    interval_idx <- tasks$interval_idx[i]
    rep_idx <- tasks$rep_idx[i]
    task_identity <- forward_checkpoint_task_identity(contract, interval_idx, rep_idx)
    settings <- contract$identity$settings
    started <- proc.time()[["elapsed"]]
    result <- simulate_interval_forward(
      interval = contract$identity$group_intervals[[interval_idx]],
      p_mis = contract$identity$fit_row$pmis,
      p_wgd = contract$identity$fit_row$pwgd,
      bottleneck_size = settings$bottleneck_size,
      expansion_factor = settings$expansion_factor,
      record_every = settings$record_every,
      divide_all_each_step = settings$divide_all_each_step,
      seed = task_identity$seed
    )
    task <- list(task_identity = task_identity, result_sha256 = forward_checkpoint_hash(result),
                 simulation_elapsed_seconds = unname(proc.time()[["elapsed"]] - started), result = result)
    forward_checkpoint_atomic_save(task, expected_paths[i])
    forward_checkpoint_task_receipt(expected_paths[i], task, reused = FALSE)
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
  out <- load_forward_checkpoint_runs(checkpoint_dir, interval_indices, replicate_indices)
  attr(out, "checkpoint_receipts") <- receipt_table
  out
}

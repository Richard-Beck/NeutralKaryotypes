# Versioned interval fitting around get_interval_nll_markov(). The caller owns
# authorization of the native physical roster; these helpers never infer
# ontology membership, paths or sample identity from a display label.

nk_fit_hash <- function(x) {
  if (!requireNamespace("digest", quietly = TRUE)) stop("Interval fitting requires 'digest'.")
  digest::digest(x, algo = "sha256", serialize = TRUE, serializeVersion = 2L)
}

nk_fit_integer <- function(x, name, minimum = 0L) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < minimum ||
      x > .Machine$integer.max || x != floor(x)) stop(name, " must be one integer >= ", minimum, ".")
  as.integer(x)
}

nk_fit_ids <- function(x, name, multiple = FALSE) {
  if (!is.character(x) || !length(x) || anyNA(x) || any(!nzchar(x)) ||
      anyDuplicated(x) || (!multiple && length(x) != 1L)) stop("Missing or ambiguous exact native ", name, ".")
  if (multiple) sort(unname(x), method = "radix") else unname(x)
}

nk_interval_physical_identity <- function(interval) {
  required <- c("connected_tree", "start", "end", "start_ids", "end_ids", "start_alias", "n_collapsed_passages")
  if (!is.list(interval) || !all(required %in% names(interval))) {
    stop("Interval identity requires exact native context, endpoints, sample IDs and proxy metadata.")
  }
  alias <- interval$start_alias
  if (is.null(alias)) alias <- character()
  if (length(alias)) alias <- nk_fit_ids(alias, "start_alias")
  collapsed <- nk_fit_integer(interval$n_collapsed_passages, "n_collapsed_passages")
  if ((length(alias) == 0L) != (collapsed == 0L)) stop("Conflicting native proxy/collapsed-path metadata.")
  role <- if (is.null(interval$comparison_type)) "native_successive" else
    nk_fit_ids(interval$comparison_type, "comparison_type")
  if (!role %in% c("native_successive", "explicit_full_window")) stop("Unknown native comparison_type: ", role)
  if (identical(interval$start, interval$end)) stop("Physical interval endpoints must be distinct.")
  path <- interval$model_path_passage_ids
  if (is.null(path)) {
    if (identical(role, "explicit_full_window")) stop("A full window requires its exact native model path.")
    path <- character()
  } else {
    # Path order is meaningful; do not sort it as an endpoint sample set.
    nk_fit_ids(path, "model_path_passage_ids", multiple = TRUE)
    path <- unname(path)
    if (!identical(path[1L], interval$start) || !identical(tail(path, 1L), interval$end) ||
        (length(alias) && (length(path) < 3L || !identical(path[2L], alias)))) {
      stop("Explicit native model path conflicts with its endpoints or proxy.")
    }
  }
  list(
    schema = "native_physical_interval_v1",
    connected_tree = nk_fit_ids(interval$connected_tree, "connected_tree"),
    start = nk_fit_ids(interval$start, "start"),
    end = nk_fit_ids(interval$end, "end"),
    start_ids = nk_fit_ids(interval$start_ids, "start_ids", multiple = TRUE),
    end_ids = nk_fit_ids(interval$end_ids, "end_ids", multiple = TRUE),
    start_alias = alias,
    n_collapsed_passages = collapsed,
    comparison_type = role,
    path_contract = if (length(path)) "explicit_native_path" else "native_endpoint_proxy_contract",
    model_path_passage_ids = path
  )
}

nk_physical_interval_id <- function(interval) {
  paste0("nki_", nk_fit_hash(nk_interval_physical_identity(interval)))
}

nk_interval_semantic_input <- function(interval) {
  matrices <- lapply(c("start_karyotype", "end_karyotype"), function(field) {
    x <- interval[[field]]
    if (!is.matrix(x) || !is.numeric(x) || !nrow(x) || !ncol(x) || any(!is.finite(x))) {
      stop("Each interval requires a finite nonempty ", field, " matrix.")
    }
    matrix(as.double(x), nrow = nrow(x), ncol = ncol(x))
  })
  if (ncol(matrices[[1L]]) != ncol(matrices[[2L]])) stop("Endpoint chromosome counts differ.")
  exposure <- interval$conditions
  if (!is.null(exposure)) {
    exposure_names <- names(exposure)
    if (!is.numeric(exposure) || any(!is.finite(exposure)) || any(exposure < 0) ||
        is.null(exposure_names) || anyNA(exposure_names) || any(!nzchar(exposure_names)) ||
        anyDuplicated(exposure_names)) stop("Invalid native path exposure arithmetic.")
    exposure <- stats::setNames(as.numeric(exposure), exposure_names)
    exposure <- exposure[order(names(exposure), method = "radix")]
  }
  list(start_karyotype = matrices[[1L]], end_karyotype = matrices[[2L]],
       n_divisions = nk_fit_integer(interval$n_divisions, "n_divisions"),
       elapsed_passages = if (is.null(interval$elapsed_passages)) NULL else
         nk_fit_integer(interval$elapsed_passages, "elapsed_passages"),
       conditions = exposure)
}

nk_interval_roster <- function(grouped_intervals) {
  ids <- names(grouped_intervals)
  if (!is.list(grouped_intervals) || !length(grouped_intervals) || is.null(ids) ||
      anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids) || any(lengths(grouped_intervals) == 0L)) {
    stop("A nonempty, uniquely named native interval group roster is required.")
  }
  identities <- list()
  rows <- lapply(seq_along(grouped_intervals), function(g) {
    do.call(rbind, lapply(seq_along(grouped_intervals[[g]]), function(i) {
      interval <- grouped_intervals[[g]][[i]]
      if (!is.null(interval$group_id) && !identical(interval$group_id, ids[g])) stop("Native interval group_id conflicts with its roster.")
      physical <- nk_interval_physical_identity(interval)
      input <- nk_interval_semantic_input(interval)
      id <- paste0("nki_", nk_fit_hash(physical))
      identities[[id]] <<- physical
      row <- data.frame(physical_interval_id = id, group_id = ids[g], interval_idx = i,
                        connected_tree = physical$connected_tree,
                        native_start = physical$start, native_end = physical$end,
                        start_alias = if (length(physical$start_alias)) physical$start_alias else NA_character_,
                        comparison_type = physical$comparison_type, path_contract = physical$path_contract,
                        n_start_cells = nrow(input$start_karyotype), n_end_cells = nrow(input$end_karyotype),
                        n_divisions = input$n_divisions, stringsAsFactors = FALSE)
      row$start_ids <- I(list(physical$start_ids))
      row$end_ids <- I(list(physical$end_ids))
      row$model_path_passage_ids <- I(list(physical$model_path_passage_ids))
      row
    }))
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  if (anyDuplicated(out$physical_interval_id)) stop("Duplicate physical comparisons in native interval roster; resolve memberships or conflicting path/division metadata upstream.")
  # A missing ordinary path cannot be used to silently disambiguate a second
  # comparison with the same endpoints and proxy but different path/role data.
  endpoint_keys <- vapply(out$physical_interval_id, function(id) {
    x <- identities[[id]]
    nk_fit_hash(x[c("connected_tree", "start", "end", "start_ids", "end_ids", "start_alias")])
  }, character(1))
  if (anyDuplicated(endpoint_keys)) stop("Ambiguous repeated physical endpoints/proxy with differing paths or comparison roles; an authorized disambiguation is required.")
  attr(out, "physical_identities") <- identities
  out
}

nk_interval_model_identity <- function(states, eps) {
  function_names <- c("get_interval_nll_markov", "hist_mat_from_matrix", "predict_marginal_cn_probs",
                      "make_cn_transition", "mat_pow", "nll_multinomial")
  definitions <- lapply(function_names, function(name) {
    fn <- get(name, mode = "function", inherits = TRUE)
    list(formals = formals(fn), body = body(fn))
  })
  list(schema = "native_interval_markov_v1", states = states, eps = eps,
       functions = stats::setNames(definitions, function_names))
}

nk_interval_runtime_identity <- function() {
  lib <- file.path(R.home("lib"), paste0(c("libR", "libRblas", "libRlapack"), .Platform$dynlib.ext))
  numerical_libraries <- c(BLAS = unname(extSoftVersion()["BLAS"]), LAPACK = La_library())
  lib <- unique(c(lib, unname(numerical_libraries)))
  lib <- lib[!is.na(lib) & file.exists(lib)]
  list(R_version = R.version.string, platform = R.version$platform,
       numerical_libraries = numerical_libraries,
       numerical_environment = Sys.getenv(c("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS")),
       runtime_sha256 = stats::setNames(vapply(lib, function(path) {
         digest::digest(file = path, algo = "sha256", serialize = FALSE)
       }, character(1)), basename(lib)))
}

nk_fit_atomic_save <- function(object, path, overwrite = FALSE) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(paste0(basename(path), ".partial-"), tmpdir = dirname(path))
  on.exit(unlink(tmp), add = TRUE)
  saveRDS(object, tmp, compress = "xz", version = 3L)
  if (!overwrite && file.exists(path)) stop("Refusing to overwrite existing fit checkpoint: ", path)
  if (!file.rename(tmp, path)) stop("Could not atomically save fit object: ", path)
  invisible(path)
}

nk_fit_read_checkpoint <- function(path, identity) {
  object <- tryCatch(readRDS(path), error = function(e) stop("Unreadable fit checkpoint: ", path, "; ", conditionMessage(e)))
  if (!is.list(object) || !identical(object$identity, identity) ||
      !identical(object$identity_sha256, nk_fit_hash(identity)) ||
      !identical(object$data_sha256, nk_fit_hash(object$data))) stop("Stale or corrupt fit checkpoint: ", path)
  object$data
}

nk_fit_checkpoint <- function(identity, compute, checkpoint_dir = NULL, family = "interval_grids") {
  sha <- nk_fit_hash(identity)
  path <- if (is.null(checkpoint_dir)) NA_character_ else file.path(checkpoint_dir, family, paste0(sha, ".Rds"))
  reused <- !is.na(path) && file.exists(path)
  data <- if (reused) nk_fit_read_checkpoint(path, identity) else compute()
  if (!reused && !is.na(path)) nk_fit_atomic_save(list(identity = identity, identity_sha256 = sha,
                                                    data = data, data_sha256 = nk_fit_hash(data)), path)
  list(data = data, sha256 = sha, checkpoint_path = path, cache_status = if (reused) "reused" else "computed")
}

nk_interval_grid_content_hash <- function(grid) {
  out <- grid[order(grid$physical_interval_id, grid$pwgd, grid$pmis), , drop = FALSE]
  # as.list.data.frame retains custom attributes, including the checksum itself.
  nk_fit_hash(stats::setNames(lapply(names(out), function(name) out[[name]]), names(out)))
}

# A complete surface is independent of the selection policy. Authenticate a
# prior surface before reuse; historical shared-WGD selected fits are never
# relabeled as independently fitted pairs.
nk_interval_reusable_surface <- function(grid, physical_interval_id, input_sha256, settings) {
  if (is.null(grid)) return(NULL)
  old_settings <- attr(grid, "settings")
  source_settings <- attr(grid, "surface_settings")
  source_schema <- attr(grid, "surface_schema")
  if (is.null(source_settings)) {
    source_settings <- old_settings
    source_schema <- "neutral_interval_grid_checkpoint_v1"
  }
  if (is.null(old_settings) || is.null(source_settings) || is.null(source_schema) ||
      !identical(old_settings$model_sha256, nk_fit_hash(attr(grid, "model_identity"))) ||
      !identical(old_settings$runtime_sha256, nk_fit_hash(attr(grid, "runtime_identity")))) {
    stop("Reusable interval grid lacks authenticated model/runtime metadata.")
  }
  part <- grid[grid$physical_interval_id == physical_interval_id, , drop = FALSE]
  if (!nrow(part)) return(NULL)
  identity <- list(schema = source_schema, physical_interval_id = physical_interval_id,
                   semantic_input_sha256 = part$semantic_input_sha256[1L], settings = source_settings)
  if (!identical(part$interval_grid_id[1L], paste0("nkig_", nk_fit_hash(identity)))) stop("Reusable interval surface identity checksum mismatch.")
  scientific <- c("pmis_grid", "pwgd_grid", "model_sha256", "runtime_sha256")
  if (!identical(old_settings[scientific], settings[scientific]) ||
      !identical(part$semantic_input_sha256[1L], input_sha256)) return(NULL)
  part <- part[order(part$pwgd, part$pmis), , drop = FALSE]
  data <- data.frame(pmis = part$pmis, pwgd = part$pwgd, negll = part$negll)
  list(data = data, source_interval_grid_id = part$interval_grid_id[1L])
}

# Every interval surface remains reusable when another interval changes. The
# historical shared mode alone couples the subsequent group profile reduction.
nk_interval_fit_grid <- function(grouped_intervals,
                                 pmis_grid = 10^(seq(-8, -1, length.out = 19)),
                                 pwgd_grid = c(0, 0.001, 0.01, 0.05),
                                 checkpoint_dir = NULL, states = 0:8, eps = 1e-9,
                                 fit_mode = "interval_shared_wgd_v1", reuse_grid = NULL) {
  fit_mode <- match.arg(fit_mode, c("interval_shared_wgd_v1", "interval_independent_v1"))
  if (!is.null(reuse_grid)) {
    if (is.null(attr(reuse_grid, "grid_sha256"))) stop("Reusable interval grid requires its content checksum.")
    nk_interval_grid_validate(reuse_grid)
  }
  roster <- nk_interval_roster(grouped_intervals)
  for (x in list(pmis_grid, pwgd_grid)) {
    if (!is.numeric(x) || !length(x) || any(!is.finite(x)) || any(x < 0) || anyDuplicated(x)) stop("Parameter grids require unique finite nonnegative values.")
  }
  if (!is.numeric(states) || !length(states) || any(!is.finite(states)) ||
      any(states != as.integer(states)) || !identical(as.integer(states), seq.int(min(states), max(states))) || min(states) < 0L) {
    stop("States must be nonnegative consecutive integers in ascending order.")
  }
  if (!is.numeric(eps) || length(eps) != 1L || !is.finite(eps) || eps <= 0) stop("eps must be finite and positive.")
  states <- as.integer(states)
  pmis_grid <- sort(as.double(pmis_grid))
  pwgd_grid <- sort(as.double(pwgd_grid))
  if (max(states) * max(pmis_grid) + max(pwgd_grid) > 1) stop("Grid violates the existing transition probability bounds.")
  params <- expand.grid(pmis = pmis_grid, pwgd = pwgd_grid, KEEP.OUT.ATTRS = FALSE)
  model <- nk_interval_model_identity(states, eps)
  runtime <- nk_interval_runtime_identity()
  settings <- list(fit_mode = fit_mode, pmis_grid = pmis_grid, pwgd_grid = pwgd_grid,
                   model_sha256 = nk_fit_hash(model), runtime_sha256 = nk_fit_hash(runtime))
  independent <- identical(fit_mode, "interval_independent_v1")
  surface_settings <- if (independent) settings[setdiff(names(settings), "fit_mode")] else settings
  surface_schema <- if (independent) "neutral_interval_grid_checkpoint_v2" else "neutral_interval_grid_checkpoint_v1"
  cache <- list()
  rows <- lapply(seq_len(nrow(roster)), function(i) {
    row <- roster[i, , drop = FALSE]
    interval <- grouped_intervals[[row$group_id]][[row$interval_idx]]
    input <- nk_interval_semantic_input(interval)
    identity <- list(schema = surface_schema, physical_interval_id = row$physical_interval_id,
                     semantic_input_sha256 = nk_fit_hash(input), settings = surface_settings)
    reused_surface <- nk_interval_reusable_surface(reuse_grid, row$physical_interval_id, identity$semantic_input_sha256, settings)
    saved <- nk_fit_checkpoint(identity, function() {
      if (!is.null(reused_surface)) return(reused_surface$data)
      out <- params
      out$negll <- vapply(seq_len(nrow(params)), function(j) {
        get_interval_nll_markov(p_mis = params$pmis[j], p_wgd = params$pwgd[j],
          K0 = input$start_karyotype, KT = input$end_karyotype, n_div = input$n_divisions, states = states, eps = eps)
      }, numeric(1))
      out
    }, checkpoint_dir)
    cache[[i]] <<- data.frame(physical_interval_id = row$physical_interval_id, interval_grid_id = paste0("nkig_", saved$sha256),
                              semantic_input_sha256 = identity$semantic_input_sha256,
                              checkpoint_path = saved$checkpoint_path,
                              cache_status = if (!is.null(reused_surface)) "reused" else saved$cache_status,
                              reuse_kind = if (saved$cache_status == "reused") "checkpoint" else if (!is.null(reused_surface)) "authenticated_surface" else "computed",
                              source_interval_grid_id = if (!is.null(reused_surface)) reused_surface$source_interval_grid_id else NA_character_)
    out <- saved$data
    out$physical_interval_id <- row$physical_interval_id
    out$group_id <- row$group_id
    out$interval_idx <- row$interval_idx
    out$interval_grid_id <- paste0("nkig_", saved$sha256)
    out$semantic_input_sha256 <- identity$semantic_input_sha256
    out$n_start_cells <- row$n_start_cells
    out$n_end_cells <- row$n_end_cells
    out$n_divisions <- row$n_divisions
    out$fit_mode <- fit_mode
    out$fit_schema <- "neutral_interval_grid_v1"
    out
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  attr(out, "roster") <- roster
  attr(out, "settings") <- settings
  attr(out, "surface_settings") <- surface_settings
  attr(out, "surface_schema") <- surface_schema
  attr(out, "model_identity") <- model
  attr(out, "runtime_identity") <- runtime
  attr(out, "interval_checkpoints") <- do.call(rbind, cache)
  attr(out, "checkpoint_dir") <- checkpoint_dir
  attr(out, "grid_sha256") <- nk_interval_grid_content_hash(out)
  out
}

nk_interval_grid_validate <- function(grid) {
  required <- c("physical_interval_id", "group_id", "interval_idx", "interval_grid_id", "semantic_input_sha256", "pmis", "pwgd", "negll",
                "n_start_cells", "n_end_cells", "n_divisions", "fit_mode", "fit_schema")
  if (!is.data.frame(grid) || !nrow(grid) || !all(required %in% names(grid)) ||
      anyNA(grid[c("physical_interval_id", "group_id", "interval_grid_id", "pmis", "pwgd")]) ||
      length(unique(grid$fit_mode)) != 1L || !all(grid$fit_mode %in% c("interval_shared_wgd_v1", "interval_independent_v1")) ||
      !all(grid$fit_schema == "neutral_interval_grid_v1")) {
    stop("Expected a versioned interval likelihood grid; pooled fits cannot be used as interval fits.")
  }
  if (!is.numeric(grid$negll) || any(!is.finite(grid$pmis)) || any(!is.finite(grid$pwgd)) ||
      any(grid$pmis < 0 | grid$pwgd < 0)) stop("Invalid interval grid numeric values.")
  if (!is.null(attr(grid, "grid_sha256")) && !identical(attr(grid, "grid_sha256"), nk_interval_grid_content_hash(grid))) stop("Interval grid content checksum mismatch.")
  params <- expand.grid(pmis = sort(unique(grid$pmis)), pwgd = sort(unique(grid$pwgd)), KEEP.OUT.ATTRS = FALSE)
  for (id in unique(grid$physical_interval_id)) {
    part <- grid[grid$physical_interval_id == id, , drop = FALSE]
    part <- part[order(part$pwgd, part$pmis), , drop = FALSE]
    if (nrow(part) != nrow(params) || !identical(as.double(part$pmis), params$pmis) ||
        !identical(as.double(part$pwgd), params$pwgd) ||
        any(vapply(part[c("group_id", "interval_idx", "interval_grid_id", "semantic_input_sha256", "n_start_cells", "n_end_cells", "n_divisions")], function(x) length(unique(x)) != 1L, logical(1)))) {
      stop("Incomplete, duplicate or conflicting interval parameter grid: ", id)
    }
  }
  invisible(TRUE)
}

# Exact equality defines ties; selection is lowest shared q, then lowest p.
# These grid profiles have no confidence-interval interpretation.
nk_interval_reduce_group <- function(data) {
  ids <- sort(unique(data$physical_interval_id), method = "radix")
  p <- sort(unique(data$pmis))
  q <- sort(unique(data$pwgd))
  surfaces <- lapply(ids, function(id) {
    d <- data[data$physical_interval_id == id, , drop = FALSE]
    d <- d[order(d$pwgd, d$pmis), , drop = FALSE]
    v <- d$negll
    v[!is.finite(v)] <- Inf
    matrix(v, nrow = length(p), ncol = length(q))
  })
  conditional <- vapply(surfaces, function(m) apply(m, 2L, min), numeric(length(q)))
  conditional <- matrix(conditional, nrow = length(q), ncol = length(ids))
  objective <- apply(conditional, 1L, sum)
  best_q <- if (any(is.finite(objective))) which(objective == min(objective)) else integer()
  chosen_q <- if (length(best_q)) best_q[1L] else NA_integer_
  group_profile <- data.frame(pwgd = q, negll = objective,
                              is_minimum = seq_along(q) %in% best_q,
                              is_selected = seq_along(q) == chosen_q)
  group_profile$is_selected[is.na(group_profile$is_selected)] <- FALSE
  interval_profiles <- do.call(rbind, lapply(seq_along(ids), function(i) {
    data.frame(physical_interval_id = ids[i], pwgd = q, negll = conditional[, i])
  }))
  pmis_profiles <- do.call(rbind, lapply(seq_along(ids), function(i) {
    others <- if (length(ids) == 1L) rep(0, length(q)) else apply(conditional[, -i, drop = FALSE], 1L, sum)
    prof <- apply(sweep(surfaces[[i]], 2L, others, "+"), 1L, min)
    data.frame(physical_interval_id = ids[i], pmis = p, negll = prof,
               delta_negll = if (any(is.finite(prof))) prof - min(prof) else rep(NA_real_, length(p)))
  }))
  selected <- do.call(rbind, lapply(seq_along(ids), function(i) {
    values <- if (is.na(chosen_q)) rep(Inf, length(p)) else surfaces[[i]][, chosen_q]
    best_p <- if (any(is.finite(values))) which(values == min(values)) else integer()
    chosen_p <- if (length(best_p)) best_p[1L] else NA_integer_
    data.frame(physical_interval_id = ids[i],
      pmis = if (length(best_p)) p[chosen_p] else NA_real_,
      pwgd = if (length(best_q)) q[chosen_q] else NA_real_,
      negll = if (length(best_p)) values[chosen_p] else Inf,
      group_negll = if (length(best_q)) objective[chosen_q] else Inf,
      pmis_tie_count = length(best_p), pwgd_tie_count = length(best_q),
      pmis_tied = length(best_p) > 1L, pwgd_tied = length(best_q) > 1L,
      pmis_flat = length(best_p) > 1L && length(best_p) == length(p),
      pwgd_flat = length(best_q) > 1L && length(best_q) == length(q),
      pmis_grid_singleton = length(p) == 1L, pwgd_grid_singleton = length(q) == 1L,
      pmis_at_lower_bound = !is.na(chosen_p) && chosen_p == 1L,
      pmis_at_upper_bound = !is.na(chosen_p) && chosen_p == length(p),
      pwgd_at_lower_bound = !is.na(chosen_q) && chosen_q == 1L,
      pwgd_at_upper_bound = !is.na(chosen_q) && chosen_q == length(q),
      n_finite_grid_points = sum(is.finite(surfaces[[i]])),
      estimability_status = if (!length(best_q)) "no_finite_shared_wgd_profile" else if (length(best_q) > 1L || length(best_p) > 1L) "tied_grid_optimum" else "unique_grid_optimum",
      uncertainty_status = "not_calibrated", stringsAsFactors = FALSE)
  }))
  list(selected = selected, pwgd_profiles = group_profile,
       conditional_interval_profiles = interval_profiles, pmis_profiles = pmis_profiles)
}

nk_interval_reduce_independent <- function(data) {
  ids <- unique(data$physical_interval_id)
  if (length(ids) != 1L) stop("Independent reduction requires exactly one physical interval.")
  data <- data[order(data$pwgd, data$pmis), , drop = FALSE]
  p <- sort(unique(data$pmis))
  q <- sort(unique(data$pwgd))
  values <- data$negll
  values[!is.finite(values)] <- Inf
  surface <- matrix(values, nrow = length(p), ncol = length(q))
  optimum <- min(values)
  # The column-major order is ascending q, then ascending p, including ties.
  best <- if (is.finite(optimum)) which(surface == optimum, arr.ind = TRUE) else matrix(integer(), ncol = 2L)
  p_index <- if (nrow(best)) best[1L, 1L] else NA_integer_
  q_index <- if (nrow(best)) best[1L, 2L] else NA_integer_
  p_profile <- apply(surface, 1L, min)
  q_profile <- apply(surface, 2L, min)
  conditional <- if (nrow(best)) surface[, q_index] else rep(Inf, length(p))
  p_ties <- if (nrow(best)) sum(conditional == optimum) else 0L
  q_ties <- if (nrow(best)) sum(q_profile == optimum) else 0L
  profile <- function(parameter, objective, index) {
    out <- data.frame(physical_interval_id = ids, value = parameter, negll = objective,
                      delta_negll = if (is.finite(optimum)) objective - optimum else NA_real_,
                      is_minimum = is.finite(objective) & objective == optimum,
                      is_selected = seq_along(parameter) == index, stringsAsFactors = FALSE)
    out$is_selected[is.na(out$is_selected)] <- FALSE
    out
  }
  pmis_profile <- profile(p, p_profile, p_index)
  names(pmis_profile)[2L] <- "pmis"
  pwgd_profile <- profile(q, q_profile, q_index)
  names(pwgd_profile)[2L] <- "pwgd"
  selected <- data.frame(physical_interval_id = ids,
    pmis = if (nrow(best)) p[p_index] else NA_real_,
    pwgd = if (nrow(best)) q[q_index] else NA_real_, negll = optimum,
    group_negll = NA_real_, joint_tie_count = nrow(best),
    pmis_tie_count = p_ties, pwgd_tie_count = q_ties,
    pmis_tied = p_ties > 1L, pwgd_tied = q_ties > 1L,
    pmis_flat = p_ties > 1L && p_ties == length(p),
    pwgd_flat = q_ties > 1L && q_ties == length(q),
    pmis_profile_tie_count = sum(pmis_profile$is_minimum),
    pmis_profile_flat = length(p) > 1L && all(pmis_profile$is_minimum),
    pmis_grid_singleton = length(p) == 1L, pwgd_grid_singleton = length(q) == 1L,
    pmis_at_lower_bound = !is.na(p_index) && p_index == 1L,
    pmis_at_upper_bound = !is.na(p_index) && p_index == length(p),
    pwgd_at_lower_bound = !is.na(q_index) && q_index == 1L,
    pwgd_at_upper_bound = !is.na(q_index) && q_index == length(q),
    n_finite_grid_points = sum(is.finite(values)),
    estimability_status = if (!nrow(best)) "no_finite_interval_grid" else if (nrow(best) > 1L) "tied_grid_optimum" else "unique_grid_optimum",
    uncertainty_status = "not_calibrated", stringsAsFactors = FALSE)
  list(selected = selected, pmis_profiles = pmis_profile, pwgd_profiles = pwgd_profile)
}

nk_interval_independent_fit_id <- function(fits, i) {
  paste0("nkif_", nk_fit_hash(list(fit_mode = fits$fit_mode[i],
    physical_interval_id = fits$physical_interval_id[i], interval_grid_id = fits$interval_grid_id[i],
    interval_profile_id = fits$interval_profile_id[i], fit_dependency_id = fits$fit_dependency_id[i],
    tie_policy = fits$tie_policy[i], pmis = fits$pmis[i], pwgd = fits$pwgd[i], negll = fits$negll[i])))
}

nk_select_independent_interval_fits <- function(grid, checkpoint_dir) {
  profiles <- list()
  caches <- list()
  result <- lapply(sort(unique(grid$physical_interval_id), method = "radix"), function(id) {
    data <- grid[grid$physical_interval_id == id, , drop = FALSE]
    data <- data[order(data$pwgd, data$pmis), , drop = FALSE]
    core_names <- c("physical_interval_id", "interval_grid_id", "pmis", "pwgd", "negll")
    core <- as.data.frame(stats::setNames(lapply(core_names, function(name) data[[name]]), core_names), stringsAsFactors = FALSE)
    identity <- list(schema = "neutral_independent_interval_profile_checkpoint_v1",
                     surface_sha256 = nk_fit_hash(core), fit_mode = "interval_independent_v1",
                     tie_policy = "exact_minimum_lowest_pwgd_then_pmis_v1",
                     reduction_sha256 = nk_fit_hash(list(formals(nk_interval_reduce_independent), body(nk_interval_reduce_independent))))
    saved <- nk_fit_checkpoint(identity, function() nk_interval_reduce_independent(core), checkpoint_dir, "interval_profiles")
    profile_id <- paste0("nkip_", saved$sha256)
    profiles[[id]] <<- saved$data[c("pmis_profiles", "pwgd_profiles")]
    caches[[id]] <<- data.frame(physical_interval_id = id, interval_profile_id = profile_id,
                                checkpoint_path = saved$checkpoint_path, cache_status = saved$cache_status)
    out <- saved$data$selected
    for (name in c("group_id", "interval_idx", "interval_grid_id", "semantic_input_sha256", "n_start_cells", "n_end_cells", "n_divisions")) out[[name]] <- data[[name]][1L]
    out$interval_profile_id <- profile_id
    out$group_profile_id <- NA_character_
    out$fit_schema <- "neutral_interval_fit_v1"
    out$fit_mode <- "interval_independent_v1"
    out$fit_scope <- "physical_interval"
    out$wgd_scope <- "physical_interval"
    out$tie_policy <- identity$tie_policy
    out$fit_dependency_id <- paste0("nkdep_", nk_fit_hash(list(fit_mode = out$fit_mode,
      physical_interval_id = id, interval_grid_id = out$interval_grid_id,
      interval_profile_id = profile_id, tie_policy = out$tie_policy)))
    out$fit_id <- nk_interval_independent_fit_id(out, 1L)
    out
  })
  out <- do.call(rbind, result)
  rownames(out) <- NULL
  attr(out, "profiles") <- profiles
  attr(out, "profile_checkpoints") <- do.call(rbind, caches)
  for (name in c("interval_checkpoints", "roster", "settings", "model_identity", "runtime_identity")) attr(out, name) <- attr(grid, name)
  nk_validate_interval_fits(out)
}

nk_select_interval_fits <- function(grid, checkpoint_dir = attr(grid, "checkpoint_dir")) {
  if (is.data.frame(grid) && "fit_schema" %in% names(grid) && all(grid$fit_schema == "neutral_interval_fit_v1")) {
    return(nk_validate_interval_fits(grid))
  }
  nk_interval_grid_validate(grid)
  if (identical(unique(grid$fit_mode), "interval_independent_v1")) return(nk_select_independent_interval_fits(grid, checkpoint_dir))
  profiles <- list()
  caches <- list()
  result <- lapply(sort(unique(grid$group_id), method = "radix"), function(group) {
    data <- grid[grid$group_id == group, , drop = FALSE]
    data <- data[order(data$physical_interval_id, data$pwgd, data$pmis), , drop = FALSE]
    core_names <- c("physical_interval_id", "interval_grid_id", "pmis", "pwgd", "negll")
    core <- as.data.frame(stats::setNames(lapply(core_names, function(name) data[[name]]), core_names),
                          stringsAsFactors = FALSE)
    identity <- list(schema = "neutral_group_wgd_profile_checkpoint_v1", surfaces_sha256 = nk_fit_hash(core),
                     tie_policy = "exact_minimum_lowest_pwgd_then_pmis_v1",
                     reduction_sha256 = nk_fit_hash(list(formals(nk_interval_reduce_group), body(nk_interval_reduce_group))))
    saved <- nk_fit_checkpoint(identity, function() nk_interval_reduce_group(core), checkpoint_dir, "group_profiles")
    profiles[[group]] <<- saved$data[c("pwgd_profiles", "conditional_interval_profiles", "pmis_profiles")]
    caches[[group]] <<- data.frame(group_id = group, group_profile_id = paste0("nkgp_", saved$sha256),
                                   checkpoint_path = saved$checkpoint_path, cache_status = saved$cache_status)
    out <- saved$data$selected
    meta <- data[match(out$physical_interval_id, data$physical_interval_id), , drop = FALSE]
    out$group_id <- group
    for (name in c("interval_idx", "interval_grid_id", "semantic_input_sha256", "n_start_cells", "n_end_cells", "n_divisions")) out[[name]] <- meta[[name]]
    out$group_profile_id <- paste0("nkgp_", saved$sha256)
    out$fit_schema <- "neutral_interval_fit_v1"
    out$fit_mode <- "interval_shared_wgd_v1"
    out$fit_scope <- "physical_interval"
    out$wgd_scope <- "group_shared_profile"
    out$tie_policy <- identity$tie_policy
    out$fit_id <- vapply(seq_len(nrow(out)), function(i) paste0("nkif_", nk_fit_hash(list(
      physical_interval_id = out$physical_interval_id[i], interval_grid_id = out$interval_grid_id[i],
      group_profile_id = out$group_profile_id[i], tie_policy = out$tie_policy[i],
      pmis = out$pmis[i], pwgd = out$pwgd[i], negll = out$negll[i]))), character(1))
    out
  })
  out <- do.call(rbind, result)
  rownames(out) <- NULL
  attr(out, "profiles") <- profiles
  attr(out, "group_checkpoints") <- do.call(rbind, caches)
  attr(out, "interval_checkpoints") <- attr(grid, "interval_checkpoints")
  attr(out, "roster") <- attr(grid, "roster")
  attr(out, "settings") <- attr(grid, "settings")
  attr(out, "model_identity") <- attr(grid, "model_identity")
  attr(out, "runtime_identity") <- attr(grid, "runtime_identity")
  nk_validate_interval_fits(out)
}

nk_validate_interval_fits <- function(fits, grouped_intervals = NULL) {
  required <- c("physical_interval_id", "group_id", "interval_idx", "pmis", "pwgd", "negll", "fit_id",
                "interval_grid_id", "semantic_input_sha256", "group_profile_id", "fit_schema", "fit_mode", "fit_scope", "wgd_scope", "tie_policy", "estimability_status")
  if (!is.data.frame(fits) || !nrow(fits) || !all(required %in% names(fits)) ||
      anyNA(fits[c("physical_interval_id", "group_id", "fit_id", "interval_grid_id", "semantic_input_sha256", "fit_mode")]) ||
      anyDuplicated(fits$physical_interval_id) || anyDuplicated(fits$fit_id) ||
      !all(fits$fit_schema == "neutral_interval_fit_v1") || length(unique(fits$fit_mode)) != 1L ||
      !all(fits$fit_mode %in% c("interval_shared_wgd_v1", "interval_independent_v1")) ||
      !all(fits$fit_scope == "physical_interval") ||
      !all(fits$tie_policy == "exact_minimum_lowest_pwgd_then_pmis_v1")) stop("Invalid or duplicate selected interval fit contract.")
  independent <- identical(unique(fits$fit_mode), "interval_independent_v1")
  if (independent) {
    if (!all(c("interval_profile_id", "fit_dependency_id", "joint_tie_count") %in% names(fits)) ||
        anyNA(fits[c("interval_profile_id", "fit_dependency_id", "joint_tie_count")]) ||
        !all(fits$wgd_scope == "physical_interval") || !all(is.na(fits$group_profile_id)) ||
        anyDuplicated(fits$interval_profile_id) || anyDuplicated(fits$fit_dependency_id)) stop("Invalid independent interval fit dependency contract.")
    expected_dependency <- vapply(seq_len(nrow(fits)), function(i) paste0("nkdep_", nk_fit_hash(list(
      fit_mode = fits$fit_mode[i], physical_interval_id = fits$physical_interval_id[i],
      interval_grid_id = fits$interval_grid_id[i], interval_profile_id = fits$interval_profile_id[i],
      tie_policy = fits$tie_policy[i]))), character(1))
    if (!identical(fits$fit_dependency_id, expected_dependency)) stop("Independent interval fit dependency checksum mismatch.")
  } else {
    if (anyNA(fits$group_profile_id) || !all(fits$wgd_scope == "group_shared_profile")) stop("Invalid shared WGD fit dependency contract.")
    for (group in unique(fits$group_id)) {
      x <- fits[fits$group_id == group, , drop = FALSE]
      if (length(unique(x$pwgd)) != 1L || length(unique(x$group_profile_id)) != 1L) stop("Conflicting shared WGD fit within group: ", group)
    }
  }
  finite <- is.finite(fits$pmis) & is.finite(fits$pwgd) & is.finite(fits$negll)
  unavailable_status <- if (independent) "no_finite_interval_grid" else "no_finite_shared_wgd_profile"
  if (!all(fits$estimability_status %in% c("unique_grid_optimum", "tied_grid_optimum", unavailable_status))) stop("Unknown interval fit estimability status.")
  unavailable <- fits$estimability_status == unavailable_status
  if (any(!finite & !unavailable) || any(finite & unavailable) || any(fits$pmis[finite] < 0 | fits$pwgd[finite] < 0)) stop("Interval fit estimability status conflicts with selected parameters.")
  expected_ids <- if (independent) vapply(seq_len(nrow(fits)), function(i) nk_interval_independent_fit_id(fits, i), character(1)) else
    vapply(seq_len(nrow(fits)), function(i) paste0("nkif_", nk_fit_hash(list(
    physical_interval_id = fits$physical_interval_id[i], interval_grid_id = fits$interval_grid_id[i],
    group_profile_id = fits$group_profile_id[i], tie_policy = fits$tie_policy[i],
    pmis = fits$pmis[i], pwgd = fits$pwgd[i], negll = fits$negll[i]))), character(1))
  if (!identical(fits$fit_id, expected_ids)) stop("Selected interval fit identity checksum mismatch.")
  if (!is.null(grouped_intervals)) {
    roster <- nk_interval_roster(grouped_intervals)
    if (!setequal(roster$physical_interval_id, fits$physical_interval_id)) stop("Selected fit and native physical interval rosters differ.")
    m <- match(fits$physical_interval_id, roster$physical_interval_id)
    if (any(fits$group_id != roster$group_id[m])) stop("Selected fit group disagrees with native interval roster.")
    fits$interval_idx <- roster$interval_idx[m]
    expected_inputs <- vapply(seq_len(nrow(fits)), function(i) {
      nk_fit_hash(nk_interval_semantic_input(grouped_intervals[[fits$group_id[i]]][[fits$interval_idx[i]]]))
    }, character(1))
    if (!identical(fits$semantic_input_sha256, expected_inputs)) stop("Selected interval fits are stale for the native endpoint data or divisions.")
  }
  fits
}

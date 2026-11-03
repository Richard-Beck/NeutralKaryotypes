# Build the ancestor path from the component root to a passage.
# Inputs: `passage_id` scalar character, `parent_map` named character vector of parent ids.
# Output: character vector ordered from root to `passage_id`.
find_passage_ancestors <- function(passage_id, parent_map) {
  out <- character()
  cur <- passage_id
  while (!is.null(cur) && !is.na(cur)) {
    out <- c(cur, out)
    if (!(cur %in% names(parent_map))) {
      break
    }
    cur <- unname(parent_map[cur])
  }
  out
}

# Find the shared prefix across several root-to-node paths.
# Inputs: `paths` list of character vectors.
# Output: character vector giving the shared prefix.
common_path_prefix <- function(paths) {
  if (!length(paths)) {
    return(character())
  }

  min_len <- min(vapply(paths, length, integer(1)))
  prefix <- character()
  for (i in seq_len(min_len)) {
    vals <- vapply(paths, `[`, character(1), i)
    if (length(unique(vals)) != 1L) {
      break
    }
    prefix <- c(prefix, vals[1])
  }
  prefix
}

# Drop the upstream part of a path so it starts at a requested node.
# Inputs: `path` character vector, `node_id` scalar character that must appear in `path`.
# Output: suffix of `path` beginning at `node_id`.
trim_path_to_node <- function(path, node_id) {
  node_idx <- match(node_id, path)
  if (is.na(node_idx)) {
    stop("Node ", node_id, " was not present in the supplied path.")
  }
  path[node_idx:length(path)]
}

# Build a connected passage tree trimmed to the MRCA of sampled passages.
# Inputs: component passage ids, precomputed path map, collapsed passaging table, and sample map.
# Output: list describing one connected tree plus sampled ids and trimmed paths.
build_component_tree <- function(component_passage_ids, path_map, passaging_tbl, sample_map) {
  component_paths <- path_map[component_passage_ids]
  common_prefix <- common_path_prefix(component_paths)
  if (!length(common_prefix)) {
    stop("Could not determine a shared ancestry path for component.")
  }

  mrca_id <- tail(common_prefix, 1)
  trimmed_paths <- lapply(component_paths, trim_path_to_node, node_id = mrca_id)
  subtree_nodes <- sort(unique(unlist(trimmed_paths)))
  subtree_edges <- passaging_tbl[
    passaging_tbl$passage_id %in% subtree_nodes,
    c("passage_id", "passage_from"),
    drop = FALSE
  ]
  subtree_edges$passage_from[!(subtree_edges$passage_from %in% subtree_nodes)] <- NA_character_
  subtree_samples <- sample_map[sample_map$passage_id %in% subtree_nodes, , drop = FALSE]

  list(
    component_root_id = common_prefix[1],
    mrca_id = mrca_id,
    mrca_depth = length(common_prefix),
    sampled_passage_ids = sort(component_passage_ids),
    sampled_ids = sort(unique(subtree_samples$samples[subtree_samples$passage_id %in% component_passage_ids])),
    tree_passage_ids = subtree_nodes,
    tree_edges = subtree_edges[order(subtree_edges$passage_id), , drop = FALSE],
    tree_samples = subtree_samples[order(subtree_samples$passage_id, subtree_samples$samples), , drop = FALSE],
    trimmed_paths = trimmed_paths
  )
}

# Convert sampled ids into connected trees over the collapsed passage graph.
# Inputs: `db_col` from `collapse_db()`, `sample_ids` character vector of karyotyped ids.
# Output: named list with `sample_map`, `connected_trees`, and `component_summary`.
build_connected_trees <- function(db_col, sample_ids) {
  sample_ids <- unique(as.character(sample_ids))
  sample_ids <- intersect(sample_ids, db_col$samples$samples)

  sample_map <- db_col$samples[match(sample_ids, db_col$samples$samples), , drop = FALSE]
  sample_map <- unique(sample_map[, c("samples", "passage_id"), drop = FALSE])
  sample_map <- sample_map[order(sample_map$passage_id, sample_map$samples), , drop = FALSE]

  passaging_tbl <- db_col$passaging
  parent_of <- setNames(passaging_tbl$passage_from, passaging_tbl$passage_id)
  unmapped_passage_ids <- setdiff(unique(sample_map$passage_id), passaging_tbl$passage_id)
  if (length(unmapped_passage_ids) > 0) {
    warning(
      length(unmapped_passage_ids),
      " sampled passage_id values are absent from db_col$passaging and will be skipped."
    )
  }

  sampled_passage_ids <- setdiff(unique(sample_map$passage_id), unmapped_passage_ids)
  passage_paths <- setNames(
    lapply(sampled_passage_ids, find_passage_ancestors, parent_map = parent_of),
    sampled_passage_ids
  )
  root_id <- vapply(passage_paths, function(path) path[1], character(1))
  connected_passage_sets <- split(sampled_passage_ids, root_id)

  connected_trees <- lapply(
    connected_passage_sets,
    build_component_tree,
    path_map = passage_paths,
    passaging_tbl = passaging_tbl,
    sample_map = sample_map
  )

  component_summary <- do.call(rbind, lapply(seq_along(connected_trees), function(i) {
    tree <- connected_trees[[i]]
    data.frame(
      component = names(connected_trees)[i],
      component_root_id = tree$component_root_id,
      mrca_id = tree$mrca_id,
      trimmed_pre_mrca = tree$component_root_id != tree$mrca_id,
      mrca_depth = tree$mrca_depth,
      n_sampled_passages = length(tree$sampled_passage_ids),
      n_sampled_ids = length(tree$sampled_ids),
      n_tree_passages = length(tree$tree_passage_ids),
      stringsAsFactors = FALSE
    )
  }))

  list(
    sample_map = sample_map,
    connected_trees = connected_trees,
    component_summary = component_summary
  )
}

# Resolve one representative media id per passage using the most common sample annotation.
# Inputs: `passage_ids` character vector, `samples_tbl` with `passage_id` and `media_id`.
# Output: data.frame with one row per passage and columns `passage_id` and `media_id`.
resolve_passage_media <- function(passage_ids, samples_tbl) {
  media_ids <- lapply(passage_ids, function(id) {
    samples_tbl$media_id[samples_tbl$passage_id %in% id]
  })
  resolved_media <- vapply(media_ids, function(ids) {
    ids <- ids[!is.na(ids) & nzchar(ids)]
    if (!length(ids)) {
      return(NA_character_)
    }
    agg <- table(ids)
    names(agg)[which.max(agg)]
  }, character(1))

  data.frame(
    passage_id = passage_ids,
    media_id = unname(resolved_media),
    stringsAsFactors = FALSE
  )
}

# Label media rows with coarse experimental conditions using hand-written rules.
# Inputs: `media_tbl` from `build_media_col()`.
# Output: `media_tbl` with an added `condition` column.
annotate_media_conditions <- function(media_tbl) {
  filters <- list(
    phosphate = !is.na(media_tbl$EnergySource2) &
      media_tbl$EnergySource2 == "Phosphates" &
      media_tbl$EnergySource2_pct < 100 &
      !is.na(media_tbl$EnergySource2_pct),
    oxygen = media_tbl$oxygen_pct < 20.5 &
      !is.na(media_tbl$oxygen_pct),
    glutamine = !is.na(media_tbl$EnergySource) &
      media_tbl$EnergySource == "L_glutamine" &
      media_tbl$EnergySource_nM < 2000000 &
      !is.na(media_tbl$EnergySource_nM),
    glucose = !is.na(media_tbl$EnergySource) &
      media_tbl$EnergySource == "Glucose" &
      media_tbl$EnergySource_nM < 11101500 &
      !is.na(media_tbl$EnergySource_nM)
  )
  filters <- do.call(cbind, filters)

  media_tbl$condition <- apply(filters, 1, function(row_flags) {
    if (!any(row_flags)) {
      return("control")
    }
    deprivations <- names(row_flags)[row_flags]
    paste(deprivations, collapse = "_")
  })

  media_tbl
}

# Attach per-passage condition annotations to each connected tree.
# Inputs: `connected_trees`, `db_samples` from `db_col$samples`, and raw `media_tbl`.
# Output: list with updated `connected_trees`, `passage_media`, and annotated `media_col`.
assign_conditions_to_trees <- function(connected_trees, db_samples, media_tbl) {
  all_ids <- unique(unlist(lapply(connected_trees, `[[`, "tree_passage_ids")))
  passage_media <- resolve_passage_media(all_ids, db_samples)
  media_col <- build_media_col(unique(passage_media$media_id), media_tbl)
  media_col <- annotate_media_conditions(media_col)
  passage_media$condition <- media_col$condition[match(passage_media$media_id, media_col$media_id)]

  connected_trees <- lapply(connected_trees, function(tree) {
    tree$tree_conditions <- passage_media[
      match(tree$tree_passage_ids, passage_media$passage_id),
      ,
      drop = FALSE
    ]
    tree$tree_conditions <- tree$tree_conditions[order(tree$tree_conditions$passage_id), , drop = FALSE]
    tree
  })

  list(
    connected_trees = connected_trees,
    passage_media = passage_media,
    media_col = media_col
  )
}

# Summarize intervals that already contain zero-copy chromosome states.
# Inputs: `simulation_intervals` list from `build_simulation_intervals()`.
# Output: data.frame or `NULL` if no zero-copy states are found.
summarize_zero_copy_states <- function(simulation_intervals) {
  out <- do.call(rbind, lapply(simulation_intervals, function(interval) {
    start_zero <- which(colSums(interval$start_karyotype == 0, na.rm = TRUE) > 0)
    end_zero <- which(colSums(interval$end_karyotype == 0, na.rm = TRUE) > 0)
    all_zero_chr <- sort(unique(c(start_zero, end_zero)))
    if (!length(all_zero_chr)) {
      return(NULL)
    }

    zero_rows <- unique(c(
      rownames(interval$start_karyotype)[rowSums(interval$start_karyotype == 0, na.rm = TRUE) > 0],
      rownames(interval$end_karyotype)[rowSums(interval$end_karyotype == 0, na.rm = TRUE) > 0]
    ))

    data.frame(
      connected_tree = interval$connected_tree,
      start = interval$start,
      end = interval$end,
      condition = interval$condition,
      zero_copy_chromosomes = paste0("chr", all_zero_chr, collapse = ", "),
      zero_copy_samples = paste(zero_rows, collapse = ", "),
      stringsAsFactors = FALSE
    )
  }))

  if (is.null(out) || !nrow(out)) {
    return(NULL)
  }
  out
}

# Pick the best-fitting parameter row within each group.
# Inputs: `group_fit_df` produced by `estimate_group_pmis_grid()`.
# Output: data.frame with one best row per group.
select_best_group_fits <- function(group_fit_df) {
  assert_legacy_pooled_fits(group_fit_df)
  best <- do.call(rbind, lapply(split(group_fit_df, group_fit_df$group_id), function(df) {
    df[which.min(df$negll), c("group_id", "pmis", "pwgd", "negll", "n_members")]
  }))
  rownames(best) <- NULL
  best
}

# Pick the best-fitting `p_mis` separately within each `(group, p_wgd)` slice.
# Inputs: `group_fit_df` produced by `estimate_group_pmis_grid()`.
# Output: data.frame with one best row per `(group_id, pwgd)` combination.
select_best_group_fits_by_pwgd <- function(group_fit_df) {
  assert_legacy_pooled_fits(group_fit_df)
  best <- do.call(rbind, lapply(
    split(group_fit_df, list(group_fit_df$group_id, group_fit_df$pwgd), drop = TRUE),
    function(df) df[which.min(df$negll), c("group_id", "pmis", "pwgd", "negll", "n_members")]
  ))
  rownames(best) <- NULL
  best
}

# The named legacy readers stay available to frozen clients, but never reduce
# interval rows into a pooled optimum or silently accept a mixed generation.
assert_legacy_pooled_fits <- function(fits) {
  if (!is.data.frame(fits) || "physical_interval_id" %in% names(fits) ||
      ("fit_mode" %in% names(fits) && anyNA(fits$fit_mode)) ||
      ("fit_mode" %in% names(fits) && any(fits$fit_mode != "pooled_group_v1"))) {
    stop("Legacy pooled fit readers cannot consume interval or mixed fit rows.")
  }
  invisible(fits)
}

select_forward_fits <- function(fit_table, grouped_intervals,
                                fit_mode = c("pooled_group_v1", "interval_shared_wgd_v1")) {
  fit_mode <- match.arg(fit_mode)
  if (fit_mode == "pooled_group_v1") return(select_best_group_fits(fit_table))
  if (!exists("nk_validate_interval_fits", mode = "function")) source("R/interval_fit_utils.R")
  selected <- nk_select_interval_fits(fit_table)
  nk_validate_interval_fits(selected, grouped_intervals)
}

# Keep opt-in interval outputs separate from frozen pooled artifacts. Explicit
# caller paths are still authoritative; these values apply only when omitted.
forward_validation_default_paths <- function(fit_mode = c("pooled_group_v1", "interval_shared_wgd_v1")) {
  fit_mode <- match.arg(fit_mode)
  interval_mode <- fit_mode == "interval_shared_wgd_v1"
  output_dir <- if (interval_mode) "results/interval_forward_validation" else "results"
  list(
    group_fit_path = if (interval_mode) "results/interval_fit_grid.Rds" else "results/group_fit_df.Rds",
    output_dir = output_dir,
    output_path = file.path(output_dir, "forward_validation.Rds"),
    summary_output_dir = if (interval_mode) "result_summaries/interval_forward_validation" else "result_summaries",
    aggregate_output_path = if (interval_mode) "data/neutral_interval_diagnostics.csv" else "data/neutral_probability.csv"
  )
}

# Build ploidy diagnostics comparing observed endpoints to Markov predictions.
# Inputs: grouped intervals and the best-per-`p_wgd` fit table.
# Output: list with `predicted` and `observed` data.frames ready for plotting.
build_ploidy_diagnostics <- function(grouped_intervals, best_by_pwgd,
                                    fit_mode = c("pooled_group_v1", "interval_shared_wgd_v1")) {
  fit_mode <- match.arg(fit_mode)
  if (fit_mode == "interval_shared_wgd_v1") {
    if (!exists("nk_validate_interval_fits", mode = "function")) source("R/interval_fit_utils.R")
    fits <- nk_validate_interval_fits(best_by_pwgd, grouped_intervals)
    if (any(!is.finite(fits$pmis) | !is.finite(fits$pwgd))) stop("Prediction diagnostics require estimable interval fits.")
    rows <- lapply(seq_len(nrow(fits)), function(i) {
      fit <- fits[i, , drop = FALSE]
      interval <- grouped_intervals[[fit$group_id]][[fit$interval_idx]]
      # Reuse the existing prediction kernel with this interval's already
      # selected pair; no nuisance refitting or regrouping is performed here.
      predicted <- predict_group_ploidy_distribution(list(interval), fit$pmis, fit$pwgd)
      observed <- stats::aggregate(list(value = rep(1, nrow(interval$end_karyotype))),
                                   list(ploidy = rowMeans(round(interval$end_karyotype))), sum)
      observed$value <- observed$value / sum(observed$value)
      predicted <- data.frame(ploidy = predicted$ploidy, value = predicted$prob)
      bind <- function(df, type) {
        df$type <- type
        for (field in c("physical_interval_id", "group_id", "interval_idx", "fit_id", "pmis", "pwgd")) df[[field]] <- fit[[field]]
        df
      }
      list(predicted = bind(predicted, "predicted"), observed = bind(observed, "observed"))
    })
    return(list(predicted = do.call(rbind, lapply(rows, `[[`, "predicted")),
                observed = do.call(rbind, lapply(rows, `[[`, "observed"))))
  }
  assert_legacy_pooled_fits(best_by_pwgd)
  ploidy_pred_diag <- do.call(rbind, lapply(seq_len(nrow(best_by_pwgd)), function(i) {
    group_id <- best_by_pwgd$group_id[i]
    pred_df <- predict_group_ploidy_distribution(
      grouped_intervals[[group_id]],
      p_mis = best_by_pwgd$pmis[i],
      p_wgd = best_by_pwgd$pwgd[i]
    )

    pred_df <- aggregate(
      list(value = pred_df$prob),
      by = list(
        group_id = rep(group_id, nrow(pred_df)),
        pwgd = rep(best_by_pwgd$pwgd[i], nrow(pred_df)),
        ploidy = pred_df$ploidy
      ),
      sum
    )
    pred_df$type <- "predicted"
    pred_df
  }))

  ploidy_obs_diag <- do.call(rbind, lapply(names(grouped_intervals), function(group_id) {
    obs_df <- observed_group_ploidy_distribution(grouped_intervals[[group_id]])
    obs_df <- aggregate(
      list(value = rep(1, nrow(obs_df))),
      by = list(group_id = rep(group_id, nrow(obs_df)), ploidy = obs_df$ploidy),
      sum
    )
    obs_df$value <- obs_df$value / sum(obs_df$value)
    obs_df
  }))

  list(
    predicted = ploidy_pred_diag,
    observed = ploidy_obs_diag
  )
}

# Build marginal chromosome diagnostics comparing fitted predictions to observed endpoints.
# Inputs: grouped intervals, best-per-`p_wgd` fit table, and chromosome indices to display.
# Output: list with `predicted` and `observed` data.frames for plotting.
build_marginal_cn_diagnostics <- function(grouped_intervals,
                                          best_by_pwgd,
                                          chr_idx = 1:2,
                                          states = 0:8,
                                          fit_mode = c("pooled_group_v1", "interval_shared_wgd_v1")) {
  fit_mode <- match.arg(fit_mode)
  if (fit_mode == "interval_shared_wgd_v1") {
    if (!exists("nk_validate_interval_fits", mode = "function")) source("R/interval_fit_utils.R")
    best_by_pwgd <- nk_validate_interval_fits(best_by_pwgd, grouped_intervals)
    if (any(!is.finite(best_by_pwgd$pmis) | !is.finite(best_by_pwgd$pwgd))) stop("Prediction diagnostics require estimable interval fits.")
  } else assert_legacy_pooled_fits(best_by_pwgd)
  chr_idx <- as.integer(chr_idx)
  chr_labels <- paste0("chr", chr_idx)

  predicted <- do.call(rbind, lapply(seq_len(nrow(best_by_pwgd)), function(i) {
    group_id <- best_by_pwgd$group_id[i]
    group_intervals <- grouped_intervals[[group_id]]

    indices <- if (fit_mode == "interval_shared_wgd_v1") best_by_pwgd$interval_idx[i] else seq_along(group_intervals)
    do.call(rbind, lapply(indices, function(interval_id) {
      interval <- group_intervals[[interval_id]]
      P_pred <- predict_marginal_cn_probs(
        K0 = interval$start_karyotype,
        n_div = interval$n_divisions,
        pmis = best_by_pwgd$pmis[i],
        pwgd = best_by_pwgd$pwgd[i],
        states = states
      )

      do.call(rbind, lapply(chr_idx, function(chr) {
        row <- data.frame(
          group_id = group_id,
          interval_id = interval_id,
          pwgd = best_by_pwgd$pwgd[i],
          chromosome = paste0("chr", chr),
          copy_number = states,
          value = P_pred[, chr],
          type = "predicted",
          stringsAsFactors = FALSE
        )
        if (fit_mode == "interval_shared_wgd_v1") {
          row$physical_interval_id <- best_by_pwgd$physical_interval_id[i]
          row$fit_id <- best_by_pwgd$fit_id[i]
          row$pmis <- best_by_pwgd$pmis[i]
        }
        row
      }))
    }))
  }))

  observed <- do.call(rbind, lapply(names(grouped_intervals), function(group_id) {
    group_intervals <- grouped_intervals[[group_id]]

    do.call(rbind, lapply(seq_along(group_intervals), function(interval_id) {
      interval <- group_intervals[[interval_id]]
      KT <- round(as.matrix(interval$end_karyotype))

      do.call(rbind, lapply(chr_idx, function(chr) {
        counts <- table(factor(KT[, chr], levels = states))
        probs <- as.numeric(counts) / sum(counts)
        row <- data.frame(
          group_id = group_id,
          interval_id = interval_id,
          chromosome = paste0("chr", chr),
          copy_number = states,
          value = probs,
          type = "observed",
          stringsAsFactors = FALSE
        )
        if (fit_mode == "interval_shared_wgd_v1") row$physical_interval_id <- nk_physical_interval_id(interval)
        row
      }))
    }))
  }))

  list(
    predicted = predicted[predicted$chromosome %in% chr_labels, , drop = FALSE],
    observed = observed[observed$chromosome %in% chr_labels, , drop = FALSE]
  )
}

# Draw a matched-size sample from a simulated endpoint matrix.
# Inputs: `simulated_matrix`, desired `target_n`, and optional RNG `seed`.
# Output: matrix with `target_n` sampled rows.
sample_simulated_endpoint <- function(simulated_matrix, target_n, seed = NULL) {
  simulated_matrix <- as.matrix(simulated_matrix)
  if (!nrow(simulated_matrix) || target_n <= 0L) {
    return(simulated_matrix[0, , drop = FALSE])
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }

  idx <- sample(
    seq_len(nrow(simulated_matrix)),
    size = as.integer(target_n),
    replace = nrow(simulated_matrix) < target_n
  )
  simulated_matrix[idx, , drop = FALSE]
}

chromosome_lengths_bp <- c(
  chr1 = 248956422,
  chr2 = 242193529,
  chr3 = 198295559,
  chr4 = 190214555,
  chr5 = 181538259,
  chr6 = 170805979,
  chr7 = 159345973,
  chr8 = 145138636,
  chr9 = 138394717,
  chr10 = 133797422,
  chr11 = 135086622,
  chr12 = 133275309,
  chr13 = 114364328,
  chr14 = 107043718,
  chr15 = 101991189,
  chr16 = 90338345,
  chr17 = 83257441,
  chr18 = 80373285,
  chr19 = 58617616,
  chr20 = 64444167,
  chr21 = 46709983,
  chr22 = 50818468
)

chromosome_distance_weights <- function(sample_matrix) {
  chr_names <- colnames(sample_matrix)
  if (is.null(chr_names) || !all(chr_names %in% names(chromosome_lengths_bp))) {
    chr_names <- paste0("chr", seq_len(ncol(sample_matrix)))
  }
  weights <- chromosome_lengths_bp[chr_names]
  if (any(is.na(weights))) {
    stop("Missing chromosome lengths for: ", paste(chr_names[is.na(weights)], collapse = ", "))
  }
  weights / sum(weights)
}

# Compute Wasserstein distance between two endpoint karyotype samples.
# Inputs: `sample_a` and `sample_b` cell-by-chromosome matrices.
# Output: scalar Wasserstein distance.
endpoint_wasserstein <- function(sample_a,
                                 sample_b,
                                 distance_metric = c("chrom_weighted_wasserstein", "wasserstein")) {
  if (!requireNamespace("transport", quietly = TRUE)) {
    stop("endpoint_wasserstein() requires the 'transport' package.")
  }

  distance_metric <- match.arg(distance_metric)
  sample_a <- round(as.matrix(sample_a))
  sample_b <- round(as.matrix(sample_b))
  if (!nrow(sample_a) || !nrow(sample_b)) {
    return(NA_real_)
  }
  if (ncol(sample_a) != ncol(sample_b)) {
    stop("Endpoint samples must have the same number of columns.")
  }
  coord_scale <- rep(1, ncol(sample_a))
  if (identical(distance_metric, "chrom_weighted_wasserstein")) {
    coord_scale <- sqrt(chromosome_distance_weights(sample_a))
  }

  make_wpp <- function(M) {
    samstr <- table(apply(M, 1, paste, collapse = "."))
    coords <- matrix(
      as.numeric(unlist(strsplit(names(samstr), split = "[.]"))),
      ncol = ncol(M),
      byrow = TRUE
    )
    coords <- sweep(coords, 2, coord_scale, `*`)
    mass <- as.numeric(samstr)
    mass <- mass / sum(mass)
    transport::wpp(coords, mass = mass)
  }

  transport::wasserstein(make_wpp(sample_a), make_wpp(sample_b))
}

# Posterior-predictive endpoint check using matched-size simulated subsamples.
# Inputs: observed endpoint matrix, a list of simulated endpoint matrices, and RNG controls.
# Output: list with observed distances, simulated null distances, summary statistics, and sampled endpoints.
posterior_predictive_endpoint_check <- function(observed_matrix,
                                                simulated_matrices,
                                                n_pairs = 100,
                                                seed = 1,
                                                distance_metric = c("chrom_weighted_wasserstein", "wasserstein"),
                                                subsample_seeds = NULL, null_pair_seed = NULL) {
  distance_metric <- match.arg(distance_metric)
  observed_matrix <- round(as.matrix(observed_matrix))
  nonempty <- vapply(simulated_matrices, function(M) nrow(as.matrix(M)) > 0L, logical(1))
  if (!is.null(subsample_seeds)) {
    if (length(subsample_seeds) != length(simulated_matrices) || anyNA(subsample_seeds) || is.null(null_pair_seed)) {
      stop("Stable endpoint-check seeds must account for every simulation replicate and null-pair selection.")
    }
    subsample_seeds <- subsample_seeds[nonempty]
  }
  simulated_matrices <- simulated_matrices[nonempty]
  n_observed <- nrow(observed_matrix)

  if (!n_observed || length(simulated_matrices) < 2L) {
    out <- list(
      observed_distances = numeric(),
      null_distances = numeric(),
      summary = list(
        n_observed = n_observed,
        n_simulations = length(simulated_matrices),
        sampled_simulated_n = n_observed,
        test_statistic = NA_real_,
        null_median = NA_real_,
        null_q25 = NA_real_,
        null_q75 = NA_real_,
        posterior_predictive_p = NA_real_,
        distance_metric = distance_metric
      ),
      ploidy = numeric()
    )
    if (!is.null(subsample_seeds)) out$simulation_rep_indices <- which(nonempty)
    return(out)
  }

  sampled_sims <- lapply(seq_along(simulated_matrices), function(i) {
    sample_simulated_endpoint(
      simulated_matrix = simulated_matrices[[i]],
      target_n = n_observed,
      seed = if (is.null(subsample_seeds)) seed + i else subsample_seeds[i]
    )
  })

  observed_distances <- vapply(sampled_sims, function(sim_sample) {
    endpoint_wasserstein(observed_matrix, sim_sample, distance_metric = distance_metric)
  }, numeric(1))
  test_statistic <- median(observed_distances, na.rm = TRUE)

  pair_idx <- utils::combn(seq_along(sampled_sims), 2L)
  if (ncol(pair_idx) > n_pairs) {
    set.seed(if (is.null(null_pair_seed)) seed + 10000L else null_pair_seed)
    pair_idx <- pair_idx[, sample(seq_len(ncol(pair_idx)), size = n_pairs, replace = FALSE), drop = FALSE]
  }
  null_distances <- vapply(seq_len(ncol(pair_idx)), function(k) {
    endpoint_wasserstein(
      sampled_sims[[pair_idx[1, k]]],
      sampled_sims[[pair_idx[2, k]]],
      distance_metric = distance_metric
    )
  }, numeric(1))

  out <- list(
    observed_distances = observed_distances,
    null_distances = null_distances,
    summary = list(
      n_observed = n_observed,
      n_simulations = length(sampled_sims),
      sampled_simulated_n = n_observed,
      test_statistic = test_statistic,
      null_median = median(null_distances, na.rm = TRUE),
      null_q25 = unname(stats::quantile(null_distances, probs = 0.25, na.rm = TRUE)),
      null_q75 = unname(stats::quantile(null_distances, probs = 0.75, na.rm = TRUE)),
      posterior_predictive_p = mean(null_distances >= test_statistic, na.rm = TRUE),
      distance_metric = distance_metric
    ),
    sampled_sims = sampled_sims
  )
  if (!is.null(subsample_seeds)) out$simulation_rep_indices <- which(nonempty)
  out
}

# Summarize forward simulations against observed endpoints for one fitted group.
# Inputs: nested `forward_runs`, the corresponding `group_intervals`, group id, and posterior-predictive settings.
# Output: list with replicate summary, interval summary, and ploidy data for plotting.
summarize_forward_validation <- function(forward_runs,
                                         group_intervals,
                                         group_id,
                                         n_null_pairs = 100,
                                         seed_base = 10000,
                                         distance_metric = c("chrom_weighted_wasserstein", "wasserstein"),
                                         fit_mode = c("pooled_group_v1", "interval_shared_wgd_v1"),
                                         fit_rows = NULL, seed_roster = NULL, checkpoint_dir = NULL) {
  fit_mode <- match.arg(fit_mode)
  distance_metric <- match.arg(distance_metric)
  interval_mode <- fit_mode == "interval_shared_wgd_v1"
  if (length(forward_runs) != length(group_intervals)) stop("Forward run and native interval rosters differ.")
  if (interval_mode) {
    if (!exists("nk_validate_interval_fits", mode = "function")) source("R/interval_fit_utils.R")
    if (!exists("forward_rng_validate", mode = "function")) source("R/forward_checkpoint_utils.R")
    fit_rows <- nk_validate_interval_fits(fit_rows, stats::setNames(list(group_intervals), group_id))
    interval_ids <- unname(vapply(group_intervals, nk_physical_interval_id, character(1)))
    fit_rows <- fit_rows[match(interval_ids, fit_rows$physical_interval_id), , drop = FALSE]
    fit_rows <- data.frame(fit_rows, check.names = FALSE, stringsAsFactors = FALSE)
    rownames(fit_rows) <- NULL
    forward_rng_validate(seed_roster)
    n_null_pairs <- forward_checkpoint_integer(n_null_pairs, "n_null_pairs")
    check_kernel <- list(version = "plugin_endpoint_check_v2",
                         functions_sha256 = forward_checkpoint_function_hash(c("posterior_predictive_endpoint_check", "sample_simulated_endpoint", "endpoint_wasserstein", "chromosome_distance_weights")),
                         chromosome_lengths_bp = chromosome_lengths_bp,
                         R_version = R.version.string, rng_kind = RNGkind(),
                         transport_version = as.character(utils::packageVersion("transport")))
  }
  interval_checks <- lapply(seq_along(forward_runs), function(interval_idx) {
    sim_matrices <- lapply(forward_runs[[interval_idx]], `[[`, "final_matrix")
    if (!interval_mode) {
      return(posterior_predictive_endpoint_check(
        observed_matrix = group_intervals[[interval_idx]]$end_karyotype,
        simulated_matrices = sim_matrices, n_pairs = n_null_pairs,
        seed = seed_base + interval_idx * 100, distance_metric = distance_metric))
    }
    id <- interval_ids[interval_idx]
    fit <- fit_rows[interval_idx, , drop = FALSE]
    for (sim in forward_runs[[interval_idx]]) {
      if (!identical(sim$p_mis, fit$pmis) || !identical(sim$p_wgd, fit$pwgd) ||
          !identical(nk_physical_interval_id(sim$interval), id) ||
          !identical(nk_interval_semantic_input(sim$interval), nk_interval_semantic_input(group_intervals[[interval_idx]]))) {
        stop("Forward simulations do not match the exact interval fit and inputs.")
      }
    }
    subsample_seeds <- forward_rng_lookup(seed_roster, id, "endpoint_subsample", seq_along(sim_matrices))
    pair_seed <- forward_rng_lookup(seed_roster, id, "endpoint_null_pairs", 0L)
    identity <- list(schema_version = "neutral_plugin_endpoint_checkpoint_v1", physical_interval_id = id,
                     fit_id = fit$fit_id, pmis = fit$pmis, pwgd = fit$pwgd,
                     observed_matrix = group_intervals[[interval_idx]]$end_karyotype,
                     simulation_sha256 = vapply(sim_matrices, forward_checkpoint_hash, character(1)),
                     n_null_pairs = n_null_pairs, distance_metric = distance_metric,
                     subsample_seeds = subsample_seeds, null_pair_seed = pair_seed, kernel = check_kernel)
    compute <- function() posterior_predictive_endpoint_check(
      observed_matrix = group_intervals[[interval_idx]]$end_karyotype,
      simulated_matrices = sim_matrices, n_pairs = n_null_pairs, distance_metric = distance_metric,
      subsample_seeds = subsample_seeds, null_pair_seed = pair_seed)
    if (is.null(checkpoint_dir)) compute() else forward_checkpoint_endpoint_check(identity, checkpoint_dir, compute)
  })

  forward_summary <- do.call(rbind, lapply(seq_along(interval_checks), function(interval_idx) {
    check <- interval_checks[[interval_idx]]
    sim_sizes <- vapply(forward_runs[[interval_idx]], function(sim) nrow(sim$final_matrix), numeric(1))
    sim_zero_frac <- vapply(forward_runs[[interval_idx]], function(sim) {
      mean(sim$final_matrix == 0, na.rm = TRUE)
    }, numeric(1))
    sim_mean_ploidy <- vapply(forward_runs[[interval_idx]], function(sim) {
      if (nrow(sim$final_matrix)) mean(rowMeans(sim$final_matrix), na.rm = TRUE) else NA_real_
    }, numeric(1))

    data.frame(
      group_id = group_id,
      interval_idx = interval_idx,
      n_observed_cells = check$summary$n_observed,
      n_simulations = check$summary$n_simulations,
      n_compared_cells = check$summary$sampled_simulated_n,
      mean_final_cells = mean(sim_sizes, na.rm = TRUE),
      mean_ploidy = mean(sim_mean_ploidy, na.rm = TRUE),
      frac_zero_chr = mean(sim_zero_frac, na.rm = TRUE),
      test_wasserstein = check$summary$test_statistic,
      null_wasserstein_q25 = check$summary$null_q25,
      null_wasserstein_median = check$summary$null_median,
      null_wasserstein_q75 = check$summary$null_q75,
      posterior_predictive_p = check$summary$posterior_predictive_p,
      distance_metric = check$summary$distance_metric,
      stringsAsFactors = FALSE
    )
  }))

  interval_summary <- forward_summary[, c(
    "group_id", "interval_idx", "n_observed_cells", "n_simulations", "n_compared_cells",
    "test_wasserstein", "null_wasserstein_q25", "null_wasserstein_median",
    "null_wasserstein_q75", "posterior_predictive_p", "distance_metric"
  )]

  simulated_ploidy <- do.call(rbind, lapply(seq_along(interval_checks), function(interval_idx) {
    check <- interval_checks[[interval_idx]]
    if (!length(check$sampled_sims)) {
      return(NULL)
    }
    do.call(rbind, lapply(seq_along(check$sampled_sims), function(rep_idx) {
      sampled_sim <- check$sampled_sims[[rep_idx]]
      data.frame(
        interval_idx = interval_idx,
        rep_idx = if (is.null(check$simulation_rep_indices)) rep_idx else check$simulation_rep_indices[rep_idx],
        ploidy = rowMeans(sampled_sim),
        type = "simulated_matched",
        stringsAsFactors = FALSE
      )
    }))
  }))

  observed_ploidy <- do.call(rbind, lapply(seq_along(group_intervals), function(interval_idx) {
    data.frame(
      interval_idx = interval_idx,
      rep_idx = NA_integer_,
      ploidy = rowMeans(group_intervals[[interval_idx]]$end_karyotype),
      type = "observed",
      stringsAsFactors = FALSE
    )
  }))

  distance_plot_df <- do.call(rbind, lapply(seq_along(interval_checks), function(interval_idx) {
    check <- interval_checks[[interval_idx]]
    rbind(
      data.frame(
        interval_idx = rep(interval_idx, length(check$observed_distances)),
        distance = check$observed_distances,
        comparison = rep("observed_vs_simulated", length(check$observed_distances)),
        distance_metric = rep(check$summary$distance_metric, length(check$observed_distances)),
        stringsAsFactors = FALSE
      ),
      data.frame(
        interval_idx = rep(interval_idx, length(check$null_distances)),
        distance = check$null_distances,
        comparison = rep("simulated_vs_simulated", length(check$null_distances)),
        distance_metric = rep(check$summary$distance_metric, length(check$null_distances)),
        stringsAsFactors = FALSE
      )
    )
  }))

  out <- list(
    replicate_summary = forward_summary,
    interval_summary = interval_summary,
    ploidy = rbind(simulated_ploidy, observed_ploidy),
    distance_plot = distance_plot_df,
    interval_checks = interval_checks
  )
  if (interval_mode) {
    identity_columns <- c("physical_interval_id", "fit_id", "fit_mode", "fit_scope", "wgd_scope", "pmis", "pwgd", "negll", "interval_grid_id", "group_profile_id", "semantic_input_sha256", "fit_schema", "estimability_status")
    identity_columns <- intersect(identity_columns, names(fit_rows))
    ontology_fields <- c("passaging_ontology_run_id", "passaging_ontology_manifest_sha256")
    ontology_values <- lapply(ontology_fields, function(field) {
      present <- vapply(group_intervals, function(interval) !is.null(interval[[field]]), logical(1))
      if (!any(present)) return(NULL)
      if (!all(present)) stop("Mixed missing ontology provenance in the interval forward roster: ", field)
      values <- vapply(group_intervals, function(interval) {
        value <- interval[[field]]
        if (!is.character(value) || length(value) != 1L || is.na(value) || !nzchar(value)) stop("Invalid ontology provenance: ", field)
        value
      }, character(1))
      if (length(unique(values)) != 1L) stop("Mixed ontology generations in the interval forward roster: ", field)
      values
    })
    names(ontology_values) <- ontology_fields
    attach_fit <- function(df) {
      for (field in identity_columns) df[[field]] <- fit_rows[[field]][df$interval_idx]
      for (field in ontology_fields) if (!is.null(ontology_values[[field]])) df[[field]] <- ontology_values[[field]][df$interval_idx]
      df
    }
    for (field in c("replicate_summary", "interval_summary", "ploidy", "distance_plot")) {
      out[[field]] <- attach_fit(out[[field]])
    }
    for (field in c("replicate_summary", "interval_summary")) {
      df <- out[[field]]
      df$diagnostic_p <- df$posterior_predictive_p
      df$posterior_predictive_p <- NULL
      df$p_value <- NA_real_
      df$q_value <- NA_real_
      df$check_method <- "plugin_endpoint_check_v2"
      df$calibration_status <- "pending_refit_calibration"
      df$inference_eligible <- FALSE
      df$neutrality_rejected <- NA
      df$forward_status <- ifelse(!is.finite(df$pmis) | !is.finite(df$pwgd), "not_simulated_non_estimable_fit",
                                  ifelse(df$n_simulations < 2L, "insufficient_surviving_simulations", "plugin_diagnostic_complete"))
      out[[field]] <- df
    }
    # A future shared-WGD calibration is coupled to the complete group inputs;
    # this pending contract is never interpreted as completed null calibration.
    ord <- order(fit_rows$physical_interval_id)
    calibration_identity <- list(schema_version = "shared_wgd_calibration_pending_v1",
                                 status = "pending_method_and_feasibility_decision",
                                 fit_rows = fit_rows[ord, setdiff(names(fit_rows), "interval_idx"), drop = FALSE],
                                 interval_inputs = unname(lapply(group_intervals[ord], nk_interval_semantic_input)))
    rownames(calibration_identity$fit_rows) <- NULL
    out$calibration <- list(status = calibration_identity$status, inference_eligible = FALSE,
                            identity_sha256 = forward_checkpoint_hash(calibration_identity), identity = calibration_identity)
    out$seed_roster <- seed_roster
  }
  out
}

# Interval-mode publication requires exact fit/check/native roster equality.
# It exports the plug-in diagnostic separately from unavailable inferential
# p/q values, so an old BH consumer cannot promote it automatically.
write_interval_neutral_probability_csv <- function(summary_dir, grouped_intervals_path, output_path, summary_pattern) {
  if (!exists("nk_validate_interval_fits", mode = "function")) source("R/interval_fit_utils.R")
  if (!file.exists(grouped_intervals_path)) stop("Interval summary publication requires the native grouped intervals.")
  grouped <- readRDS(grouped_intervals_path)
  roster <- nk_interval_roster(grouped)
  paths <- list.files(summary_dir, pattern = summary_pattern, full.names = TRUE)
  if (!length(paths)) stop("No interval summaries found for the required native roster.")
  tables <- lapply(paths, function(path) {
    res <- readRDS(path)
    if (!identical(res$settings$fit_mode, "interval_shared_wgd_v1") ||
        !res$forward_group_id %in% names(grouped)) {
      stop("Interval summary writer rejects legacy, mixed or unknown-group summaries: ", path)
    }
    fits <- nk_validate_interval_fits(res$fit_row, grouped[res$forward_group_id])
    df <- res$interval_summary
    parity <- c("fit_id", "group_id", "pmis", "pwgd", "negll", "fit_mode", "fit_scope", "wgd_scope", "interval_grid_id", "group_profile_id", "semantic_input_sha256", "fit_schema", "estimability_status")
    required <- c("physical_interval_id", parity,
                  "diagnostic_p", "p_value", "q_value", "check_method", "calibration_status", "inference_eligible", "neutrality_rejected")
    if (!is.data.frame(df) || !all(required %in% names(df)) || anyDuplicated(df$physical_interval_id) ||
        !setequal(df$physical_interval_id, fits$physical_interval_id)) {
      stop("Interval fit/check roster or schema mismatch: ", path)
    }
    idx <- match(df$physical_interval_id, fits$physical_interval_id)
    for (field in parity) {
      if (!identical(df[[field]], fits[[field]][idx])) stop("Interval summary/fit numeric or identity mismatch in ", field, ": ", path)
    }
    if (anyNA(df$inference_eligible) || any(df$inference_eligible) || any(!is.na(df$p_value)) || any(!is.na(df$q_value)) ||
        any(!is.na(df$neutrality_rejected)) || any(df$calibration_status != "pending_refit_calibration") ||
        any(df$check_method != "plugin_endpoint_check_v2")) {
      stop("Plug-in endpoint summaries cannot be promoted as calibrated inference: ", path)
    }
    df$result_file <- basename(path)
    df
  })
  out <- do.call(rbind, tables)
  if (anyDuplicated(out$physical_interval_id) || !setequal(out$physical_interval_id, roster$physical_interval_id)) {
    stop("Published interval summaries must equal the complete native physical roster exactly once.")
  }
  for (field in intersect(c("passaging_ontology_run_id", "passaging_ontology_manifest_sha256"), names(out))) {
    if (anyNA(out[[field]]) || length(unique(out[[field]])) != 1L) stop("Mixed ontology provenance in interval summary publication: ", field)
  }
  out <- out[order(out$physical_interval_id), , drop = FALSE]
  idx <- match(out$physical_interval_id, roster$physical_interval_id)
  if (any(out$group_id != roster$group_id[idx])) stop("Summary group identities disagree with native physical intervals.")
  out$interval_idx <- roster$interval_idx[idx]
  intervals <- lapply(seq_len(nrow(out)), function(i) grouped[[out$group_id[i]]][[out$interval_idx[i]]])
  out$condition <- out$group_id
  out$replicate_id <- vapply(intervals, `[[`, character(1), "end")
  out$ancestor <- vapply(intervals, `[[`, character(1), "start")
  out$delta_pass <- vapply(intervals, function(interval) {
    if (is.null(interval$elapsed_passages)) NA_real_ else as.numeric(interval$elapsed_passages)
  }, numeric(1))
  out$cin_rate <- out$pmis
  out$p_wgd <- out$pwgd
  out$fit_negll <- out$negll
  out$summary_schema <- "neutral_interval_forward_summary_v1"
  rownames(out) <- NULL
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(out, output_path, row.names = FALSE)
  message("Wrote ", nrow(out), " interval plug-in diagnostic rows; calibrated inference remains pending: ", output_path)
  out
}

# Aggregate cached forward-validation summaries into a flat neutrality results table.
# Inputs: summary directory/path pattern, grouped intervals, and output CSV path.
# Output: data.frame written to `output_path`, compatible with the old neutral_probability.csv columns.
write_neutral_probability_csv <- function(summary_dir = forward_validation_default_paths(fit_mode)$summary_output_dir,
                                          grouped_intervals_path = "core_data/grouped_intervals.Rds",
                                          output_path = forward_validation_default_paths(fit_mode)$aggregate_output_path,
                                          summary_pattern = "^forward_validation_.*_summary[.]Rds$",
                                          fit_mode = c("pooled_group_v1", "interval_shared_wgd_v1")) {
  fit_mode <- match.arg(fit_mode)
  if (fit_mode == "interval_shared_wgd_v1") {
    return(write_interval_neutral_probability_csv(summary_dir, grouped_intervals_path, output_path, summary_pattern))
  }
  out_cols <- c(
    "condition", "replicate_id", "ancestor", "p_value", "cin_rate", "delta_pass",
    "p_wgd", "fit_negll", "interval_idx", "test_wasserstein",
    "null_wasserstein_median", "distance_metric", "result_file"
  )

  forward_result_files <- list.files(summary_dir, pattern = summary_pattern, full.names = TRUE)
  grouped_intervals <- if (file.exists(grouped_intervals_path)) {
    readRDS(grouped_intervals_path)
  } else {
    list()
  }

  if (!length(forward_result_files)) {
    out <- data.frame(matrix(ncol = length(out_cols), nrow = 0))
    names(out) <- out_cols
    dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(out, output_path, row.names = FALSE)
    message("No cached forward-validation summaries found; wrote empty ", output_path)
    return(out)
  }

  forward_validation_summary <- do.call(rbind, lapply(forward_result_files, function(path) {
    res <- readRDS(path)
    if ((!is.null(res$settings$fit_mode) && !identical(res$settings$fit_mode, "pooled_group_v1")) ||
        "physical_interval_id" %in% names(res$interval_summary)) {
      stop("Legacy pooled summary writer cannot consume interval or mixed summaries: ", path)
    }
    df <- res$interval_summary
    if (is.null(df) || !nrow(df)) {
      return(NULL)
    }
    if (is.null(df$distance_metric)) {
      df$distance_metric <- if (!is.null(res$settings$distance_metric)) {
        res$settings$distance_metric
      } else {
        "wasserstein"
      }
    }

    df$forward_group_id <- res$forward_group_id
    df$result_file <- basename(path)
    if (!is.null(res$fit_row) && nrow(res$fit_row)) {
      df$cin_rate <- res$fit_row$pmis[1]
      df$p_wgd <- res$fit_row$pwgd[1]
      df$fit_negll <- res$fit_row$negll[1]
    } else {
      df$cin_rate <- NA_real_
      df$p_wgd <- NA_real_
      df$fit_negll <- NA_real_
    }

    if (res$forward_group_id %in% names(grouped_intervals)) {
      interval_meta <- grouped_intervals[[res$forward_group_id]]
      df$condition <- res$forward_group_id
      df$replicate_id <- vapply(df$interval_idx, function(i) interval_meta[[i]]$end, character(1))
      df$ancestor <- vapply(df$interval_idx, function(i) interval_meta[[i]]$start, character(1))
      df$delta_pass <- vapply(df$interval_idx, function(i) interval_meta[[i]]$elapsed_passages, numeric(1))
    } else {
      df$condition <- res$forward_group_id
      df$replicate_id <- NA_character_
      df$ancestor <- NA_character_
      df$delta_pass <- NA_real_
    }

    df
  }))

  if (is.null(forward_validation_summary) || !nrow(forward_validation_summary)) {
    out <- data.frame(matrix(ncol = length(out_cols), nrow = 0))
    names(out) <- out_cols
  } else {
    forward_validation_summary$p_value <- forward_validation_summary$posterior_predictive_p
    forward_validation_summary <- forward_validation_summary[
      order(forward_validation_summary$forward_group_id, forward_validation_summary$interval_idx),
      ,
      drop = FALSE
    ]
    out <- forward_validation_summary[, out_cols]
  }

  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(out, output_path, row.names = FALSE)
  message("Wrote ", nrow(out), " aggregate neutrality result rows to ", output_path)
  out
}

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
    !is.na(media_tbl$EnergySource2) & media_tbl$EnergySource2_pct < 100,
    media_tbl$oxygen_pct < 20.5,
    media_tbl$EnergySource == "L_glutamine" &
      media_tbl$EnergySource_nM < 2000000 &
      !is.na(media_tbl$EnergySource_nM),
    media_tbl$EnergySource == "Glucose" &
      media_tbl$EnergySource_nM < 1110150 &
      !is.na(media_tbl$EnergySource_nM)
  )
  filters <- do.call(cbind, filters)

  media_tbl$condition <- apply(filters, 1, function(row_flags) {
    if (sum(row_flags) == 0L) {
      return("control")
    }
    deprivations <- c("phosphate", "oxygen", "glutamine", "glucose")[row_flags]
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
  best <- do.call(rbind, lapply(
    split(group_fit_df, list(group_fit_df$group_id, group_fit_df$pwgd), drop = TRUE),
    function(df) df[which.min(df$negll), c("group_id", "pmis", "pwgd", "negll", "n_members")]
  ))
  rownames(best) <- NULL
  best
}

# Build ploidy diagnostics comparing observed endpoints to Markov predictions.
# Inputs: grouped intervals and the best-per-`p_wgd` fit table.
# Output: list with `predicted` and `observed` data.frames ready for plotting.
build_ploidy_diagnostics <- function(grouped_intervals, best_by_pwgd) {
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
                                          states = 0:8) {
  chr_idx <- as.integer(chr_idx)
  chr_labels <- paste0("chr", chr_idx)

  predicted <- do.call(rbind, lapply(seq_len(nrow(best_by_pwgd)), function(i) {
    group_id <- best_by_pwgd$group_id[i]
    group_intervals <- grouped_intervals[[group_id]]

    do.call(rbind, lapply(seq_along(group_intervals), function(interval_id) {
      interval <- group_intervals[[interval_id]]
      P_pred <- predict_marginal_cn_probs(
        K0 = interval$start_karyotype,
        n_div = interval$n_divisions,
        pmis = best_by_pwgd$pmis[i],
        pwgd = best_by_pwgd$pwgd[i],
        states = states
      )

      do.call(rbind, lapply(chr_idx, function(chr) {
        data.frame(
          group_id = group_id,
          interval_id = interval_id,
          pwgd = best_by_pwgd$pwgd[i],
          chromosome = paste0("chr", chr),
          copy_number = states,
          value = P_pred[, chr],
          type = "predicted",
          stringsAsFactors = FALSE
        )
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
        data.frame(
          group_id = group_id,
          interval_id = interval_id,
          chromosome = paste0("chr", chr),
          copy_number = states,
          value = probs,
          type = "observed",
          stringsAsFactors = FALSE
        )
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

# Compute Wasserstein distance between two endpoint karyotype samples.
# Inputs: `sample_a` and `sample_b` cell-by-chromosome matrices.
# Output: scalar Wasserstein distance.
endpoint_wasserstein <- function(sample_a, sample_b) {
  if (!requireNamespace("transport", quietly = TRUE)) {
    stop("endpoint_wasserstein() requires the 'transport' package.")
  }

  sample_a <- round(as.matrix(sample_a))
  sample_b <- round(as.matrix(sample_b))
  if (!nrow(sample_a) || !nrow(sample_b)) {
    return(NA_real_)
  }
  if (ncol(sample_a) != ncol(sample_b)) {
    stop("Endpoint samples must have the same number of columns.")
  }

  make_wpp <- function(M) {
    samstr <- table(apply(M, 1, paste, collapse = "."))
    coords <- matrix(
      as.numeric(unlist(strsplit(names(samstr), split = "[.]"))),
      ncol = ncol(M),
      byrow = TRUE
    )
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
                                                seed = 1) {
  observed_matrix <- round(as.matrix(observed_matrix))
  simulated_matrices <- Filter(function(M) nrow(as.matrix(M)) > 0, simulated_matrices)
  n_observed <- nrow(observed_matrix)

  if (!n_observed || length(simulated_matrices) < 2L) {
    return(list(
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
        posterior_predictive_p = NA_real_
      ),
      ploidy = numeric()
    ))
  }

  sampled_sims <- lapply(seq_along(simulated_matrices), function(i) {
    sample_simulated_endpoint(
      simulated_matrix = simulated_matrices[[i]],
      target_n = n_observed,
      seed = seed + i
    )
  })

  observed_distances <- vapply(sampled_sims, function(sim_sample) {
    endpoint_wasserstein(observed_matrix, sim_sample)
  }, numeric(1))
  test_statistic <- median(observed_distances, na.rm = TRUE)

  pair_idx <- utils::combn(seq_along(sampled_sims), 2L)
  if (ncol(pair_idx) > n_pairs) {
    set.seed(seed + 10000L)
    pair_idx <- pair_idx[, sample(seq_len(ncol(pair_idx)), size = n_pairs, replace = FALSE), drop = FALSE]
  }
  null_distances <- vapply(seq_len(ncol(pair_idx)), function(k) {
    endpoint_wasserstein(
      sampled_sims[[pair_idx[1, k]]],
      sampled_sims[[pair_idx[2, k]]]
    )
  }, numeric(1))

  list(
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
      posterior_predictive_p = mean(null_distances >= test_statistic, na.rm = TRUE)
    ),
    sampled_sims = sampled_sims
  )
}

# Summarize forward simulations against observed endpoints for one fitted group.
# Inputs: nested `forward_runs`, the corresponding `group_intervals`, group id, and posterior-predictive settings.
# Output: list with replicate summary, interval summary, and ploidy data for plotting.
summarize_forward_validation <- function(forward_runs,
                                         group_intervals,
                                         group_id,
                                         n_null_pairs = 100,
                                         seed_base = 10000) {
  interval_checks <- lapply(seq_along(forward_runs), function(interval_idx) {
    sim_matrices <- lapply(forward_runs[[interval_idx]], `[[`, "final_matrix")
    posterior_predictive_endpoint_check(
      observed_matrix = group_intervals[[interval_idx]]$end_karyotype,
      simulated_matrices = sim_matrices,
      n_pairs = n_null_pairs,
      seed = seed_base + interval_idx * 100
    )
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
      stringsAsFactors = FALSE
    )
  }))

  interval_summary <- forward_summary[, c(
    "group_id", "interval_idx", "n_observed_cells", "n_simulations", "n_compared_cells",
    "test_wasserstein", "null_wasserstein_q25", "null_wasserstein_median",
    "null_wasserstein_q75", "posterior_predictive_p"
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
        rep_idx = rep_idx,
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
        interval_idx = interval_idx,
        distance = check$observed_distances,
        comparison = "observed_vs_simulated",
        stringsAsFactors = FALSE
      ),
      data.frame(
        interval_idx = interval_idx,
        distance = check$null_distances,
        comparison = "simulated_vs_simulated",
        stringsAsFactors = FALSE
      )
    )
  }))

  list(
    replicate_summary = forward_summary,
    interval_summary = interval_summary,
    ploidy = rbind(simulated_ploidy, observed_ploidy),
    distance_plot = distance_plot_df,
    interval_checks = interval_checks
  )
}

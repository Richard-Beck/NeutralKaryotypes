#!/usr/bin/env Rscript

parse_args <- function(args) {
  opts <- list(
    csv = "data/neutral_probability.csv",
    out_dir = "figures",
    output_name = "08_replicate_grouped_neglog10p.png",
    annotated_csv = "neutral_probability_plot_annotated.csv",
    karyotypes_path = "core_data/karyotypes.Rds",
    db_col_path = "core_data/db_col.Rds",
    media_path = "core_data/media_raw.Rds",
    karyotyped_samples_path = "core_data/karyotyped_samples.txt",
    heatmap_subdir = "cell_line_heatmaps",
    heatmap_enrichment_csv = "cluster_condition_enrichment.csv",
    heatmap_cluster_selection_csv = "cluster_selection_summary.csv",
    heatmap_enrichment_alpha = 0.05,
    heatmap_width = 2200,
    heatmap_height = 1400,
    heatmap_distance = "chrom_weighted",
    heatmap_cluster_mode = "auto",
    heatmap_cluster_score = "silhouette",
    heatmap_k_min = 2,
    heatmap_k_max = 10,
    width = 1980,
    height = 1152,
    res = 180
  )

  for (arg in args) {
    if (!startsWith(arg, "--") || !grepl("=", arg, fixed = TRUE)) {
      stop("Arguments must use --name=value syntax. Invalid argument: ", arg)
    }
    parts <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1]]
    key <- parts[1]
    value <- paste(parts[-1], collapse = "=")
    if (!key %in% names(opts)) {
      stop("Unknown argument: --", key)
    }
    opts[[key]] <- value
  }

  opts$width <- as.integer(opts$width)
  opts$height <- as.integer(opts$height)
  opts$heatmap_width <- as.integer(opts$heatmap_width)
  opts$heatmap_height <- as.integer(opts$heatmap_height)
  opts$heatmap_enrichment_alpha <- as.numeric(opts$heatmap_enrichment_alpha)
  opts$heatmap_k_min <- as.integer(opts$heatmap_k_min)
  opts$heatmap_k_max <- as.integer(opts$heatmap_k_max)
  opts$res <- as.integer(opts$res)
  opts
}

required_columns <- c(
  "condition",
  "replicate_id",
  "ancestor",
  "p_value",
  "cin_rate",
  "delta_pass"
)

normalize_cell_line_label <- function(value) {
  text <- as.character(value)

  if (grepl("^MDA-231", text) || grepl("^MDA-MB-231", text)) {
    return("MDA-231")
  }
  if (grepl("^SUM-?159", text)) {
    ploidy_match <- regexec("SUM-?159.*_([0-9]N)(_|$)", text)
    ploidy_hit <- regmatches(text, ploidy_match)[[1]]
    if (length(ploidy_hit) >= 2) {
      if (grepl("(^|_)NLS(_|$)", text)) {
        return(paste0("SUM-159_NLS_", ploidy_hit[2]))
      }
      return(paste0("SUM-159_", ploidy_hit[2]))
    }

    ploidy_match <- regexpr("SUM-?159_[0-9]N", text)
    if (ploidy_match[1] > 0) {
      return(gsub("^SUM159", "SUM-159", regmatches(text, ploidy_match)))
    }
    return("SUM-159")
  }
  if (grepl("^SNU-668", text)) {
    return("SNU-668")
  }
  if (grepl("^HGC-27", text)) {
    return("HGC-27")
  }

  strsplit(text, "_", fixed = TRUE)[[1]][1]
}

infer_cell_line <- function(condition, replicate_id) {
  condition <- as.character(condition)
  replicate_id <- as.character(replicate_id)
  condition_prefix <- strsplit(condition, "||", fixed = TRUE)[[1]][1]

  if (nzchar(condition_prefix)) {
    if (grepl("^SUM-?159", condition_prefix)) {
      prefix_cell_line <- normalize_cell_line_label(condition_prefix)
      if (prefix_cell_line != "SUM-159") {
        return(prefix_cell_line)
      }
    }

    match <- regexec("^(.+?)_([0-9]N|[0-9]+)(_|$)", condition_prefix)
    hit <- regmatches(condition_prefix, match)[[1]]
    if (length(hit) >= 2) {
      return(normalize_cell_line_label(hit[2]))
    }
    if (grepl("||", condition, fixed = TRUE)) {
      prefix_cell_line <- normalize_cell_line_label(condition_prefix)
      if (!(grepl("^SUM-?159", condition_prefix) && prefix_cell_line == "SUM-159")) {
        return(prefix_cell_line)
      }
    }
  }

  normalize_cell_line_label(replicate_id)
}

canonicalize_condition <- function(value) {
  text <- as.character(value)
  parts <- strsplit(text, "||", fixed = TRUE)[[1]]
  if (length(parts) > 1) {
    text <- parts[length(parts)]
  }

  normalized <- gsub("[-_]", " ", trimws(text))
  lower <- tolower(normalized)

  known <- c(
    control = "Control",
    glucose = "Glucose",
    oxygen = "Oxygen",
    phosphate = "Phosphate"
  )
  if (lower %in% names(known)) {
    return(unname(known[[lower]]))
  }

  rev_match <- regexpr("\\brev[a-z0-9]+\\b", lower, perl = TRUE)
  if (rev_match[1] > 0) {
    token <- regmatches(lower, rev_match)
    return(paste0("Rev", toupper(sub("^rev", "", token))))
  }

  if (!nzchar(normalized)) {
    return("Unknown")
  }
  paste0(toupper(substr(normalized, 1, 1)), substr(normalized, 2, nchar(normalized)))
}

stage_from_id <- function(sample_id) {
  text <- tolower(as.character(sample_id))
  if (grepl("harvest", text, fixed = TRUE)) {
    return("Harvest")
  }
  if (grepl("seed", text, fixed = TRUE)) {
    return("Seed")
  }
  if (grepl("mrca", text, fixed = TRUE)) {
    return("MRCA")
  }
  "Other"
}

add_plot_columns <- function(df) {
  missing <- setdiff(required_columns, names(df))
  if (length(missing)) {
    stop("neutral_probability CSV is missing required columns: ", paste(missing, collapse = ", "))
  }

  df$p_value <- as.numeric(df$p_value)
  df$cin_rate <- as.numeric(df$cin_rate)
  df$delta_pass <- as.numeric(df$delta_pass)

  if (anyNA(df$p_value)) {
    stop("p_value contains missing or non-numeric values.")
  }
  if (any(df$p_value < 0 | df$p_value > 1)) {
    stop("p_value must be between 0 and 1.")
  }

  positive_p <- df$p_value[df$p_value > 0]
  zero_floor <- if (length(positive_p)) min(positive_p) / 2 else 1e-300
  df$p_for_log <- ifelse(df$p_value <= 0, zero_floor, df$p_value)
  df$neglog10_p <- -log10(df$p_for_log)
  df$cell_line <- mapply(infer_cell_line, df$condition, df$replicate_id, USE.NAMES = FALSE)
  df$condition_group <- vapply(df$condition, canonicalize_condition, character(1))
  df$start_stage <- vapply(df$ancestor, stage_from_id, character(1))
  df$end_stage <- vapply(df$replicate_id, stage_from_id, character(1))
  df$main_group_key <- paste(df$cell_line, df$condition_group, sep = " | ")
  df
}

ordered_unique <- function(x) {
  unique(as.character(x))
}

chr_lengths_bp <- c(
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

condition_display_order <- function(values) {
  preferred <- c("Control", "Oxygen", "Glucose", "Phosphate")
  extras <- setdiff(sort(unique(as.character(values))), preferred)
  c(preferred[preferred %in% values], extras)
}

build_sample_condition_lookup <- function(neutral_probability_df) {
  endpoint_lookup <- setNames(
    as.character(neutral_probability_df$condition_group),
    as.character(neutral_probability_df$replicate_id)
  )

  ancestor_rows <- neutral_probability_df[
    !is.na(neutral_probability_df$ancestor) &
      nzchar(as.character(neutral_probability_df$ancestor)),
    ,
    drop = FALSE
  ]
  ancestor_lookup <- tapply(
    as.character(ancestor_rows$condition_group),
    as.character(ancestor_rows$ancestor),
    function(values) {
      condition_values <- sort(unique(values[!is.na(values) & nzchar(values)]))
      if (length(condition_values) == 1L) {
        return(condition_values)
      }
      "Control"
    }
  )

  lookup <- unlist(ancestor_lookup, use.names = TRUE)
  lookup[names(endpoint_lookup)] <- endpoint_lookup
  lookup
}

ordered_condition_groups <- function(groups, group_to_condition) {
  condition_order <- condition_display_order(unname(group_to_condition[groups]))
  group_conditions <- unname(group_to_condition[groups])
  condition_rank <- match(group_conditions, condition_order)
  condition_rank[is.na(condition_rank)] <- length(condition_order) + seq_len(sum(is.na(condition_rank)))
  groups[order(condition_rank, seq_along(groups))]
}

build_x_axis <- function(df) {
  group_min <- tapply(df$p_for_log, df$main_group_key, min)
  globally_ordered_groups <- names(sort(group_min, method = "radix"))
  group_to_condition <- setNames(
    df$condition_group[match(globally_ordered_groups, df$main_group_key)],
    globally_ordered_groups
  )
  group_to_cell_line <- setNames(
    df$cell_line[match(globally_ordered_groups, df$main_group_key)],
    globally_ordered_groups
  )
  ordered_cell_lines <- unique(group_to_cell_line)
  ordered_groups <- unlist(
    lapply(
      ordered_cell_lines,
      function(cell_line) {
        cell_groups <- globally_ordered_groups[group_to_cell_line == cell_line]
        ordered_condition_groups(cell_groups, group_to_condition)
      }
    ),
    use.names = FALSE
  )

  keys <- character()
  labels <- character()

  for (group in ordered_groups) {
    group_df <- df[df$main_group_key == group, , drop = FALSE]
    lineage_min <- tapply(group_df$p_for_log, group_df$condition, min)
    ordered_lineages <- names(sort(lineage_min, method = "radix"))

    for (lineage in ordered_lineages) {
      lineage_df <- group_df[group_df$condition == lineage, , drop = FALSE]
      endpoints <- sort(unique(as.character(lineage_df$replicate_id)))
      new_keys <- paste(group, lineage, endpoints, sep = "\r")
      keys <- c(keys, new_keys)
      labels <- c(labels, endpoints)
    }
  }

  list(keys = keys, labels = labels)
}

plot_grouped_neglog10p <- function(df, output_path, width, height, res) {
  axis_info <- build_x_axis(df)
  x_keys <- axis_info$keys
  x_labels <- axis_info$labels
  x_pos <- seq_along(x_keys)
  names(x_pos) <- x_keys

  cell_lines <- ordered_unique(df$cell_line[order(df$cell_line, df$condition_group)])
  cell_palette <- grDevices::hcl.colors(max(3, length(cell_lines)), palette = "Dark 3")
  cell_colors <- setNames(cell_palette[seq_along(cell_lines)], cell_lines)

  group_condition <- df[!duplicated(df$main_group_key), c("main_group_key", "condition_group")]
  group_to_condition <- setNames(group_condition$condition_group, group_condition$main_group_key)
  x_groups <- vapply(strsplit(x_keys, "\r", fixed = TRUE), `[`, character(1), 1)
  x_conditions <- unname(group_to_condition[x_groups])
  condition_levels <- ordered_unique(x_conditions)
  condition_palette <- grDevices::hcl.colors(max(3, length(condition_levels)), palette = "Set 2")
  condition_colors <- setNames(condition_palette[seq_along(condition_levels)], condition_levels)

  grDevices::png(output_path, width = width, height = height, res = res)
  old_par <- par(no.readonly = TRUE)
  on.exit({
    par(old_par)
    grDevices::dev.off()
  }, add = TRUE)
  par(mar = c(12, 5, 7, 14), xpd = NA)

  ylim_top <- max(df$neglog10_p, -log10(0.01), na.rm = TRUE) * 1.12
  plot(
    NA,
    xlim = c(0.5, length(x_keys) + 0.5),
    ylim = c(0, ylim_top),
    xaxt = "n",
    xlab = "",
    ylab = "-log10(p)",
    bty = "l"
  )

  set.seed(42)
  for (cell_line in cell_lines) {
    cell_df <- df[df$cell_line == cell_line, , drop = FALSE]
    point_keys <- paste(cell_df$main_group_key, cell_df$condition, cell_df$replicate_id, sep = "\r")
    xs <- unname(x_pos[point_keys]) + rnorm(nrow(cell_df), mean = 0, sd = 0.08)
    points(
      xs,
      cell_df$neglog10_p,
      pch = 16,
      cex = 1.0,
      col = grDevices::adjustcolor(cell_colors[[cell_line]], alpha.f = 0.82)
    )
  }

  previous_group <- x_groups[1]
  for (i in seq_along(x_groups)) {
    if (i > 1 && x_groups[i] != previous_group) {
      abline(v = i - 0.5, col = "grey85", lty = "dotted", lwd = 1)
    }
    previous_group <- x_groups[i]
  }

  abline(h = -log10(0.05), col = "grey35", lty = "dashed", lwd = 1)
  abline(h = -log10(0.01), col = "grey20", lty = "dotted", lwd = 1)

  tick_step <- max(1, ceiling(length(x_labels) / 60))
  visible <- ((seq_along(x_labels) - 1) %% tick_step) == 0
  axis(
    side = 1,
    at = x_pos[visible],
    labels = FALSE,
    cex.axis = 0.58,
    tick = TRUE
  )
  label_y <- par("usr")[3] - 0.02 * diff(par("usr")[3:4])
  text(
    x = x_pos[visible],
    y = label_y,
    labels = x_labels[visible],
    srt = 90,
    adj = 1,
    cex = 0.58,
    col = condition_colors[x_conditions[visible]]
  )
  mtext("Endpoint replicate, grouped by cell line and condition", side = 1, line = 10)

  usr <- par("usr")
  strip_y0 <- usr[4] + 0.035 * diff(usr[3:4])
  strip_y1 <- usr[4] + 0.075 * diff(usr[3:4])
  for (i in seq_along(x_pos)) {
    rect(
      xleft = x_pos[i] - 0.5,
      ybottom = strip_y0,
      xright = x_pos[i] + 0.5,
      ytop = strip_y1,
      border = NA,
      col = grDevices::adjustcolor(condition_colors[[x_conditions[i]]], alpha.f = 0.9)
    )
  }

  legend(
    "topright",
    inset = c(-0.28, 0),
    legend = cell_lines,
    col = cell_colors[cell_lines],
    pch = 16,
    title = "Cell line",
    bty = "n",
    cex = 0.9
  )
  legend(
    x = mean(usr[1:2]),
    y = usr[4] + 0.28 * diff(usr[3:4]),
    legend = condition_levels,
    fill = condition_colors[condition_levels],
    title = "Condition strip",
    bty = "n",
    horiz = TRUE,
    xjust = 0.5,
    cex = 0.85
  )
}

read_karyotyped_sample_ids <- function(path) {
  ids <- readLines(path, warn = FALSE)
  unique(trimws(ids[nzchar(trimws(ids))]))
}

extract_passage_rank <- function(sample_id, passage_id) {
  candidates <- c(as.character(sample_id), as.character(passage_id))
  if (any(grepl("MRCA", candidates, fixed = TRUE))) {
    return(-1)
  }

  passage_match <- regexec("_A([0-9]+)", candidates)
  hits <- regmatches(candidates, passage_match)
  nums <- vapply(hits, function(hit) {
    if (length(hit) >= 2) {
      return(as.integer(hit[2]))
    }
    NA_integer_
  }, integer(1))

  if (any(!is.na(nums))) {
    return(min(nums, na.rm = TRUE))
  }

  Inf
}

build_karyotype_cell_metadata <- function(karyotypes_path,
                                          db_col_path,
                                          media_path,
                                          karyotyped_samples_path,
                                          neutral_probability_df) {
  source("R/db_utils.R")
  source("R/workflow_utils.R")

  karyotypes <- readRDS(karyotypes_path)
  db_col <- readRDS(db_col_path)
  media_raw <- readRDS(media_path)
  karyotyped_ids <- read_karyotyped_sample_ids(karyotyped_samples_path)

  tree_build <- build_connected_trees(db_col, karyotyped_ids)
  condition_build <- assign_conditions_to_trees(
    connected_trees = tree_build$connected_trees,
    db_samples = db_col$samples,
    media_tbl = media_raw
  )
  connected_trees <- condition_build$connected_trees

  tree_membership <- do.call(rbind, lapply(seq_along(connected_trees), function(i) {
    tree <- connected_trees[[i]]
    data.frame(
      passage_id = tree$tree_passage_ids,
      connected_tree = names(connected_trees)[i],
      stringsAsFactors = FALSE
    )
  }))

  condition_lookup <- do.call(rbind, lapply(connected_trees, function(tree) {
    tree$tree_conditions[, c("passage_id", "condition"), drop = FALSE]
  }))
  condition_lookup <- unique(condition_lookup)

  sample_meta <- unique(tree_build$sample_map[, c("samples", "passage_id"), drop = FALSE])
  names(sample_meta)[names(sample_meta) == "samples"] <- "sample_id"
  sample_meta$sample_id <- as.character(sample_meta$sample_id)
  sample_meta$passage_id <- as.character(sample_meta$passage_id)
  sample_meta$connected_tree <- tree_membership$connected_tree[
    match(sample_meta$passage_id, tree_membership$passage_id)
  ]
  sample_meta$condition_raw <- condition_lookup$condition[
    match(sample_meta$passage_id, condition_lookup$passage_id)
  ]
  neutral_condition_lookup <- build_sample_condition_lookup(neutral_probability_df)
  sample_meta$condition_group <- unname(neutral_condition_lookup[sample_meta$passage_id])
  missing_condition <- is.na(sample_meta$condition_group) | !nzchar(sample_meta$condition_group)
  if (any(missing_condition)) {
    missing_passages <- sort(unique(sample_meta$passage_id[missing_condition]))
    warning(
      "Dropping heatmap passages absent from neutral_probability.csv replicate_id values: ",
      paste(missing_passages, collapse = ", "),
      call. = FALSE
    )
    sample_meta <- sample_meta[!missing_condition, , drop = FALSE]
  }
  sample_meta$cell_line <- mapply(
    infer_cell_line,
    sample_meta$connected_tree,
    sample_meta$sample_id,
    USE.NAMES = FALSE
  )
  sample_meta$passage_rank <- vapply(
    seq_len(nrow(sample_meta)),
    function(i) extract_passage_rank(sample_meta$sample_id[i], sample_meta$passage_id[i]),
    numeric(1)
  )
  sample_meta$passage_label <- sample_meta$passage_id

  cell_tbl <- lapply(split(karyotypes, as.character(karyotypes$id)), function(rows) {
    kmat <- do.call(rbind, lapply(rows$karyotype, function(k) as.numeric(unlist(k))))
    if (ncol(kmat) < 22) {
      stop("Expected karyotype vectors with at least 22 chromosomes.")
    }
    kmat <- kmat[, 1:22, drop = FALSE]
    out <- data.frame(
      sample_id = rep(as.character(rows$id[1]), nrow(kmat)),
      cell_index = seq_len(nrow(kmat)),
      kmat,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    out
  })
  cell_tbl <- do.call(rbind, cell_tbl)
  names(cell_tbl)[3:ncol(cell_tbl)] <- paste0("chr", seq_len(22))

  out <- merge(sample_meta, cell_tbl, by = "sample_id", all.x = TRUE, sort = FALSE)
  out <- out[!is.na(out$chr1), , drop = FALSE]
  out$cell_label <- paste0(out$sample_id, "#", out$cell_index)
  out
}

draw_discrete_legend <- function(title, labels, colors, x, y_top, cex = 0.75, box_width = 0.8, box_height = 0.5) {
  if (!length(labels)) {
    return(invisible(NULL))
  }

  text(x, y_top + 0.5, labels = title, adj = c(0, 0), font = 2, cex = cex)
  for (i in seq_along(labels)) {
    y_mid <- y_top - (i - 1) * 0.75
    rect(
      xleft = x,
      ybottom = y_mid - box_height / 2,
      xright = x + box_width,
      ytop = y_mid + box_height / 2,
      col = colors[[labels[i]]],
      border = NA,
      xpd = NA
    )
    text(x + box_width + 0.2, y_mid, labels = labels[i], adj = c(0, 0.5), cex = cex, xpd = NA)
  }
}

draw_continuous_legend <- function(colors, zlim, x, y_bottom, y_top, ticks = pretty(zlim, n = 5), cex = 0.75) {
  y_seq <- seq(y_bottom, y_top, length.out = length(colors) + 1)
  for (i in seq_along(colors)) {
    rect(
      xleft = x,
      ybottom = y_seq[i],
      xright = x + 0.8,
      ytop = y_seq[i + 1],
      col = colors[i],
      border = NA,
      xpd = NA
    )
  }

  text(x, y_top + 0.8, labels = "Mean CN", adj = c(0, 0), font = 2, cex = cex, xpd = NA)
  tick_pos <- y_bottom + (ticks - zlim[1]) / diff(zlim) * (y_top - y_bottom)
  segments(x + 0.8, tick_pos, x + 1.0, tick_pos, xpd = NA)
  text(x + 1.2, tick_pos, labels = format(ticks, trim = TRUE), adj = c(0, 0.5), cex = cex, xpd = NA)
}

compute_cell_distance <- function(cell_mat, mode = "chrom_weighted") {
  mode <- match.arg(mode, c("chrom_weighted", "unweighted"))
  if (nrow(cell_mat) <= 1) {
    return(stats::dist(cell_mat))
  }

  if (mode == "unweighted") {
    return(stats::dist(cell_mat))
  }

  chr_names <- colnames(cell_mat)
  weights <- chr_lengths_bp[chr_names]
  if (any(is.na(weights))) {
    stop("Missing chromosome lengths for: ", paste(chr_names[is.na(weights)], collapse = ", "))
  }
  weights <- weights / sum(weights)
  scaled_mat <- sweep(cell_mat, 2, sqrt(weights), `*`)
  stats::dist(scaled_mat)
}

calinski_harabasz_score <- function(cell_mat, membership) {
  k <- length(unique(membership))
  n <- nrow(cell_mat)
  if (k <= 1 || k >= n) {
    return(NA_real_)
  }

  overall_center <- colMeans(cell_mat)
  cluster_ids <- sort(unique(membership))

  between_ss <- 0
  within_ss <- 0
  for (cluster_id in cluster_ids) {
    cluster_mat <- cell_mat[membership == cluster_id, , drop = FALSE]
    cluster_n <- nrow(cluster_mat)
    if (cluster_n == 0) {
      next
    }
    cluster_center <- colMeans(cluster_mat)
    between_ss <- between_ss + cluster_n * sum((cluster_center - overall_center) ^ 2)
    within_ss <- within_ss + sum(rowSums((cluster_mat - matrix(cluster_center, nrow = cluster_n, ncol = ncol(cell_mat), byrow = TRUE)) ^ 2))
  }

  if (within_ss <= 0) {
    return(NA_real_)
  }
  (between_ss / (k - 1)) / (within_ss / (n - k))
}

gap_statistic_score <- function(cell_mat, membership, B = 10L) {
  k <- length(unique(membership))
  n <- nrow(cell_mat)
  if (k <= 1 || k >= n) {
    return(NA_real_)
  }

  cluster_dispersion <- function(mat, groups) {
    total <- 0
    for (cluster_id in unique(groups)) {
      cluster_mat <- mat[groups == cluster_id, , drop = FALSE]
      if (nrow(cluster_mat) <= 1) {
        next
      }
      center <- colMeans(cluster_mat)
      total <- total + sum(sqrt(rowSums((cluster_mat - matrix(center, nrow = nrow(cluster_mat), ncol = ncol(mat), byrow = TRUE)) ^ 2)))
    }
    total
  }

  wk_obs <- cluster_dispersion(cell_mat, membership)
  if (!is.finite(wk_obs) || wk_obs <= 0) {
    return(NA_real_)
  }

  mins <- apply(cell_mat, 2, min)
  maxs <- apply(cell_mat, 2, max)
  ref_logs <- replicate(B, {
    ref_mat <- sapply(seq_len(ncol(cell_mat)), function(j) stats::runif(nrow(cell_mat), min = mins[j], max = maxs[j]))
    ref_mat <- as.matrix(ref_mat)
    ref_hc <- stats::hclust(stats::dist(ref_mat), method = "ward.D")
    ref_groups <- stats::cutree(ref_hc, k = k)
    wk_ref <- cluster_dispersion(ref_mat, ref_groups)
    log(wk_ref)
  })

  mean(ref_logs) - log(wk_obs)
}

evaluate_cluster_condition_enrichment <- function(cell_df,
                                                  cluster_membership,
                                                  condition_levels,
                                                  alpha = 0.05,
                                                  adjust_method = "BH") {
  total_cells <- nrow(cell_df)
  clusters <- sort(unique(as.integer(cluster_membership)))
  out <- do.call(rbind, lapply(clusters, function(cluster_id) {
    in_cluster <- cluster_membership == cluster_id
    cluster_size <- sum(in_cluster)

    do.call(rbind, lapply(condition_levels, function(condition_name) {
      is_condition <- as.character(cell_df$condition_group) == condition_name
      a <- sum(in_cluster & is_condition)
      b <- sum(in_cluster & !is_condition)
      c <- sum(!in_cluster & is_condition)
      d <- sum(!in_cluster & !is_condition)
      expected <- cluster_size * sum(is_condition) / total_cells
      enrichment_fold <- if (expected > 0) a / expected else NA_real_
      fisher_p <- stats::fisher.test(
        matrix(c(a, b, c, d), nrow = 2),
        alternative = "greater"
      )$p.value

      data.frame(
        cluster_id = cluster_id,
        condition = condition_name,
        n_cells_total = total_cells,
        n_cells_cluster = cluster_size,
        n_condition_total = sum(is_condition),
        n_condition_in_cluster = a,
        expected_in_cluster = expected,
        enrichment_fold = enrichment_fold,
        p_value = fisher_p,
        stringsAsFactors = FALSE
      )
    }))
  }))

  out$p_adj <- stats::p.adjust(out$p_value, method = adjust_method)
  out$is_overrepresented <- with(
    out,
    n_condition_in_cluster > expected_in_cluster & p_adj <= alpha
  )
  out
}

select_cluster_count <- function(distance_obj,
                                 hc,
                                 cell_mat,
                                 cluster_mode = "auto",
                                 cluster_score = "silhouette",
                                 k_min = 2L,
                                 k_max = 10L,
                                 fixed_k = NULL) {
  cluster_mode <- match.arg(cluster_mode, c("auto", "fixed"))
  cluster_score <- match.arg(cluster_score, c("silhouette", "calinski_harabasz", "gap"))
  n_cells <- attr(distance_obj, "Size")
  if (n_cells <= 1 || is.null(hc)) {
    summary_df <- data.frame(
      k = 1L,
      score = NA_real_,
      selected = TRUE,
      method = cluster_mode,
      score_method = cluster_score,
      stringsAsFactors = FALSE
    )
    return(list(k = 1L, summary = summary_df))
  }

  if (cluster_mode == "fixed") {
    if (is.null(fixed_k) || !length(fixed_k) || is.na(fixed_k)) {
      stop("fixed_k must be provided when cluster_mode='fixed'.")
    }
    chosen_k <- max(1L, min(as.integer(fixed_k[1]), n_cells))
    summary_df <- data.frame(
      k = chosen_k,
      score = NA_real_,
      selected = TRUE,
      method = cluster_mode,
      score_method = cluster_score,
      stringsAsFactors = FALSE
    )
    return(list(k = chosen_k, summary = summary_df))
  }

  max_candidate <- min(as.integer(k_max), n_cells - 1L)
  min_candidate <- min(as.integer(k_min), max_candidate)
  if (max_candidate < 2L || min_candidate < 2L || min_candidate > max_candidate) {
    summary_df <- data.frame(
      k = 1L,
      score = NA_real_,
      selected = TRUE,
      method = cluster_mode,
      score_method = cluster_score,
      stringsAsFactors = FALSE
    )
    return(list(k = 1L, summary = summary_df))
  }

  candidate_ks <- seq.int(min_candidate, max_candidate)
  scores <- vapply(candidate_ks, function(k) {
    membership <- stats::cutree(hc, k = k)
    if (cluster_score == "silhouette") {
      sil <- cluster::silhouette(membership, dmatrix = as.matrix(distance_obj))
      return(mean(sil[, "sil_width"]))
    }
    if (cluster_score == "calinski_harabasz") {
      return(calinski_harabasz_score(cell_mat, membership))
    }
    gap_statistic_score(cell_mat, membership)
  }, numeric(1))

  chosen_k <- candidate_ks[which.max(scores)]
  summary_df <- data.frame(
    k = candidate_ks,
    score = scores,
    selected = candidate_ks == chosen_k,
    method = cluster_mode,
    score_method = cluster_score,
    stringsAsFactors = FALSE
  )
  list(k = chosen_k, summary = summary_df)
}

plot_column_dendrogram <- function(hclust_obj, x_positions, y_bottom, y_top, line_col = "grey20") {
  if (is.null(hclust_obj) || length(x_positions) <= 1) {
    return(invisible(NULL))
  }

  dend <- stats::as.dendrogram(hclust_obj)
  x_map <- setNames(as.numeric(x_positions), labels(dend))
  height_max <- max(hclust_obj$height, na.rm = TRUE)

  if (!is.finite(height_max) || height_max <= 0) {
    height_max <- 1
  }

  draw_node <- function(node) {
    if (is.leaf(node)) {
      return(list(x = x_map[[attr(node, "label")]], y = y_bottom))
    }

    kids <- lapply(node, draw_node)
    child_x <- vapply(kids, `[[`, numeric(1), "x")
    child_y <- vapply(kids, `[[`, numeric(1), "y")
    node_height <- attr(node, "height")
    y_here <- y_bottom + (node_height / height_max) * (y_top - y_bottom)

    for (i in seq_along(kids)) {
      segments(child_x[i], child_y[i], child_x[i], y_here, col = line_col, xpd = NA)
    }
    segments(min(child_x), y_here, max(child_x), y_here, col = line_col, xpd = NA)

    list(x = mean(range(child_x)), y = y_here)
  }

  draw_node(dend)
}

plot_cell_line_heatmap <- function(cell_df,
                                   output_path,
                                   width,
                                   height,
                                   res,
                                   distance_mode,
                                   enrichment_alpha,
                                   cluster_mode,
                                   cluster_score,
                                   k_min,
                                   k_max) {
  if (!nrow(cell_df)) {
    return(invisible(NULL))
  }

  condition_levels <- condition_display_order(cell_df$condition_group)
  cell_df$condition_group <- factor(cell_df$condition_group, levels = condition_levels)
  cell_df <- cell_df[order(
    cell_df$condition_group,
    cell_df$passage_rank,
    cell_df$passage_label,
    cell_df$sample_id,
    cell_df$cell_index
  ), , drop = FALSE]

  chr_cols <- paste0("chr", seq_len(22))
  cell_mat <- as.matrix(cell_df[, chr_cols, drop = FALSE])
  rownames(cell_mat) <- cell_df$cell_label

  hc <- NULL
  cluster_membership <- rep(1L, nrow(cell_mat))
  selected_k <- 1L
    cluster_selection_summary <- data.frame(
      k = 1L,
      score = NA_real_,
      selected = TRUE,
      method = cluster_mode,
      score_method = cluster_score,
      stringsAsFactors = FALSE
    )
  if (nrow(cell_mat) > 1) {
    d <- compute_cell_distance(cell_mat, mode = distance_mode)
    hc <- stats::hclust(d, method = "ward.D")
    ord <- hc$order
    cell_df <- cell_df[ord, , drop = FALSE]
    cell_mat <- cell_mat[ord, , drop = FALSE]

    fixed_k <- if (identical(unique(as.character(cell_df$cell_line)), "MDA-231")) 4L else 6L
    k_choice <- select_cluster_count(
      distance_obj = d,
      hc = hc,
      cell_mat = cell_mat,
      cluster_mode = cluster_mode,
      cluster_score = cluster_score,
      k_min = k_min,
      k_max = k_max,
      fixed_k = fixed_k
    )
    selected_k <- k_choice$k
    cluster_selection_summary <- k_choice$summary
    cluster_membership <- stats::cutree(hc, k = selected_k)[ord]
  }

  enrichment_df <- evaluate_cluster_condition_enrichment(
    cell_df = cell_df,
    cluster_membership = cluster_membership,
    condition_levels = condition_levels,
    alpha = enrichment_alpha
  )

  heatmap_mat <- t(cell_mat)
  rownames(heatmap_mat) <- chr_cols
  colnames(heatmap_mat) <- cell_df$cell_label
  n_rows <- nrow(heatmap_mat)
  n_cols <- ncol(heatmap_mat)

  condition_palette <- grDevices::hcl.colors(max(3, length(condition_levels)), palette = "Set 2")
  condition_colors <- setNames(condition_palette[seq_along(condition_levels)], condition_levels)

  passage_levels <- unique(cell_df$passage_label)
  passage_palette <- grDevices::colorRampPalette(c("#8C5A2B", "#F5F1E8", "#111111"))(max(3, length(passage_levels)))
  passage_colors <- setNames(passage_palette[seq_along(passage_levels)], passage_levels)

  zlim <- range(heatmap_mat, na.rm = TRUE)
  if (!all(is.finite(zlim)) || diff(zlim) == 0) {
    zlim <- c(zlim[1] - 0.5, zlim[2] + 0.5)
  }
  fill_colors <- grDevices::colorRampPalette(c("#2166AC", "#F7F7BF", "#B2182B"))(100)

  map_to_fill <- function(values) {
    idx <- floor((values - zlim[1]) / diff(zlim) * (length(fill_colors) - 1)) + 1
    idx[!is.finite(idx)] <- 1
    idx <- pmax(1, pmin(length(fill_colors), idx))
    fill_colors[idx]
  }

  grDevices::png(output_path, width = width, height = height, res = res)
  old_par <- par(no.readonly = TRUE)
  on.exit({
    par(old_par)
    grDevices::dev.off()
  }, add = TRUE)

  n_condition_rows <- length(condition_levels)
  star_y_base <- n_rows + 1.8
  star_y_step <- 0.5
  dend_y0 <- star_y_base + n_condition_rows * star_y_step + 0.3
  dend_y1 <- dend_y0 + 2.4

  par(mar = c(7, 7, 10, 19), xpd = NA)
  plot(
    NA,
    xlim = c(0.5, n_cols + 12.5),
    ylim = c(0.5, dend_y1 + 0.8),
    xaxt = "n",
    yaxt = "n",
    xlab = "",
    ylab = "",
    bty = "n"
  )

  y_pos <- rev(seq_len(n_rows))
  x_pos <- seq_len(n_cols)
  plot_column_dendrogram(hc, x_pos, dend_y0, dend_y1)

  for (row_idx in seq_len(n_rows)) {
    y <- y_pos[row_idx]
    row_colors <- map_to_fill(heatmap_mat[row_idx, ])
    for (col_idx in seq_len(n_cols)) {
      x <- x_pos[col_idx]
      rect(x - 0.5, y - 0.5, x + 0.5, y + 0.5, col = row_colors[col_idx], border = "grey92")
    }
    text(0.3, y, labels = rownames(heatmap_mat)[row_idx], adj = 1, cex = 0.7)
  }

  cond_y0 <- n_rows + 0.15
  cond_y1 <- n_rows + 0.75
  pass_y0 <- n_rows + 0.85
  pass_y1 <- n_rows + 1.45

  for (col_idx in seq_len(n_cols)) {
    x <- x_pos[col_idx]
    cond_col <- condition_colors[[as.character(cell_df$condition_group[col_idx])]]
    pass_col <- passage_colors[[cell_df$passage_label[col_idx]]]
    rect(x - 0.5, cond_y0, x + 0.5, cond_y1, col = cond_col, border = NA)
    rect(x - 0.5, pass_y0, x + 0.5, pass_y1, col = pass_col, border = NA)
  }

  axis(side = 2, at = y_pos, labels = rownames(heatmap_mat), las = 2, cex.axis = 0.8)
  axis(side = 1, at = x_pos, labels = FALSE, tick = FALSE)
  text(0.3, (cond_y0 + cond_y1) / 2, labels = "Cond", adj = 1, cex = 0.7, font = 2)
  text(0.3, (pass_y0 + pass_y1) / 2, labels = "Pass", adj = 1, cex = 0.7, font = 2)

  sig_df <- enrichment_df[enrichment_df$is_overrepresented, , drop = FALSE]
  if (nrow(sig_df)) {
    for (i in seq_along(condition_levels)) {
      condition_name <- condition_levels[i]
      y_star <- star_y_base + (n_condition_rows - i) * star_y_step
      text(
        0.3,
        y_star,
        labels = condition_name,
        adj = 1,
        cex = 0.62,
        col = condition_colors[[condition_name]],
        font = 2
      )

      cond_sig <- sig_df[sig_df$condition == condition_name, , drop = FALSE]
      if (!nrow(cond_sig)) {
        next
      }

      for (j in seq_len(nrow(cond_sig))) {
        cluster_id <- cond_sig$cluster_id[j]
        xs <- x_pos[cluster_membership == cluster_id]
        if (!length(xs)) {
          next
        }
        text(
          mean(range(xs)),
          y_star,
          labels = "*",
          cex = 1.1,
          col = condition_colors[[condition_name]],
          font = 2
        )
      }
    }
  }

  cluster_breaks <- which(cluster_membership[-1] != cluster_membership[-n_cols])
  for (idx in cluster_breaks) {
    x_sep <- x_pos[idx] + 0.5
    segments(x_sep, 0.5, x_sep, dend_y1, lwd = 1.2, col = "grey55")
  }

  title(
    main = paste0(unique(cell_df$cell_line), ": cell-level karyotype heatmap"),
    sub = paste0(
      "Columns are observed cells clustered with ward.D using ",
      if (identical(distance_mode, "chrom_weighted")) "chromosome-length-weighted" else "unweighted",
      " distance and cut into ",
      selected_k,
      " clusters; rows are chromosomes; * marks cluster-condition enrichments"
    )
  )

  draw_discrete_legend(
    title = "Condition",
    labels = condition_levels,
    colors = condition_colors,
    x = n_cols + 2.2,
    y_top = dend_y1 - 0.2
  )
  draw_discrete_legend(
    title = "Passage",
    labels = passage_levels,
    colors = passage_colors,
    x = n_cols + 2.2,
    y_top = n_rows - 2.5,
    cex = 0.65
  )
  draw_continuous_legend(
    colors = fill_colors,
    zlim = zlim,
    x = n_cols + 9.0,
    y_bottom = max(1, n_rows * 0.15),
    y_top = max(6, n_rows * 0.9)
  )

  enrichment_df$cell_line <- unique(as.character(cell_df$cell_line))
  enrichment_df$distance_mode <- distance_mode
  enrichment_df$n_clusters <- selected_k
  cluster_selection_summary$cell_line <- unique(as.character(cell_df$cell_line))
  cluster_selection_summary$distance_mode <- distance_mode
  cluster_selection_summary$selected_k <- selected_k

  list(
    enrichment = enrichment_df,
    cluster_selection = cluster_selection_summary
  )
}

write_cell_line_heatmaps <- function(output_dir,
                                     karyotypes_path,
                                     db_col_path,
                                     media_path,
                                     karyotyped_samples_path,
                                     neutral_probability_df,
                                     width,
                                     height,
                                     res,
                                     distance_mode,
                                     enrichment_csv_name,
                                     enrichment_alpha,
                                     cluster_selection_csv_name,
                                     cluster_mode,
                                     cluster_score,
                                     k_min,
                                     k_max) {
  required_paths <- c(karyotypes_path, db_col_path, media_path, karyotyped_samples_path)
  if (!all(file.exists(required_paths))) {
    missing <- required_paths[!file.exists(required_paths)]
    warning(
      "Skipping heatmap generation because required cached inputs are missing: ",
      paste(missing, collapse = ", ")
    )
    return(invisible(NULL))
  }

  cell_profiles <- build_karyotype_cell_metadata(
    karyotypes_path = karyotypes_path,
    db_col_path = db_col_path,
    media_path = media_path,
    karyotyped_samples_path = karyotyped_samples_path,
    neutral_probability_df = neutral_probability_df
  )
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  enrichment_results <- list()
  cluster_selection_results <- list()
  for (cell_line in ordered_unique(cell_profiles$cell_line)) {
    cell_df <- cell_profiles[cell_profiles$cell_line == cell_line, , drop = FALSE]
    output_path <- file.path(
      output_dir,
      paste0(gsub("[^A-Za-z0-9._-]+", "_", cell_line), "_karyotype_heatmap.png")
    )
    heatmap_result <- plot_cell_line_heatmap(
      cell_df = cell_df,
      output_path = output_path,
      width = width,
      height = height,
      res = res,
      distance_mode = distance_mode,
      enrichment_alpha = enrichment_alpha,
      cluster_mode = cluster_mode,
      cluster_score = cluster_score,
      k_min = k_min,
      k_max = k_max
    )
    enrichment_results[[cell_line]] <- heatmap_result$enrichment
    cluster_selection_results[[cell_line]] <- heatmap_result$cluster_selection
    message("Saved cell-line heatmap to ", output_path)
  }

  enrichment_df <- do.call(rbind, enrichment_results)
  cluster_selection_df <- do.call(rbind, cluster_selection_results)
  if (!is.null(enrichment_df) && nrow(enrichment_df)) {
    enrichment_df$alpha <- enrichment_alpha
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(
      enrichment_df,
      file.path(output_dir, enrichment_csv_name),
      row.names = FALSE
    )
    message("Saved cluster-condition enrichment results to ", file.path(output_dir, enrichment_csv_name))
  }

  if (!is.null(cluster_selection_df) && nrow(cluster_selection_df)) {
    utils::write.csv(
      cluster_selection_df,
      file.path(output_dir, cluster_selection_csv_name),
      row.names = FALSE
    )
    message("Saved cluster-selection summary to ", file.path(output_dir, cluster_selection_csv_name))
  }
}

main <- function() {
  opts <- parse_args(commandArgs(trailingOnly = TRUE))
  df <- read.csv(opts$csv, stringsAsFactors = FALSE, check.names = FALSE)
  annotated <- add_plot_columns(df)

  dir.create(opts$out_dir, showWarnings = FALSE, recursive = TRUE)
  output_path <- file.path(opts$out_dir, opts$output_name)
  annotated_path <- file.path(opts$out_dir, opts$annotated_csv)

  write.csv(annotated, annotated_path, row.names = FALSE)
  plot_grouped_neglog10p(annotated, output_path, opts$width, opts$height, opts$res)
  write_cell_line_heatmaps(
    output_dir = file.path(opts$out_dir, opts$heatmap_subdir),
    karyotypes_path = opts$karyotypes_path,
    db_col_path = opts$db_col_path,
    media_path = opts$media_path,
    karyotyped_samples_path = opts$karyotyped_samples_path,
    neutral_probability_df = annotated,
    width = opts$heatmap_width,
    height = opts$heatmap_height,
    res = opts$res,
    distance_mode = opts$heatmap_distance,
    enrichment_csv_name = opts$heatmap_enrichment_csv,
    enrichment_alpha = opts$heatmap_enrichment_alpha,
    cluster_selection_csv_name = opts$heatmap_cluster_selection_csv,
    cluster_mode = opts$heatmap_cluster_mode,
    cluster_score = opts$heatmap_cluster_score,
    k_min = opts$heatmap_k_min,
    k_max = opts$heatmap_k_max
  )

  message("Saved grouped -log10(p) plot to ", output_path)
  message("Saved annotated plotting data to ", annotated_path)
}

called_as_script <- function() {
  file_args <- sub("^--file=", "", commandArgs(trailingOnly = FALSE))
  any(grepl("(^|/)plot_neutral_probability_neglog10p[.]R$", file_args))
}

if (called_as_script()) {
  main()
}

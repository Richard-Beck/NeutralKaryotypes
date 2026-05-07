#!/usr/bin/env Rscript

parse_args <- function(args) {
  opts <- list(
    csv = "data/neutral_probability.csv",
    out_dir = "figures",
    output_name = "08_replicate_grouped_neglog10p.png",
    annotated_csv = "neutral_probability_plot_annotated.csv",
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
  if (grepl("^SUM-159", text)) {
    ploidy_match <- regexec("SUM-159.*_([0-9]N)(_|$)", text)
    ploidy_hit <- regmatches(text, ploidy_match)[[1]]
    if (length(ploidy_hit) >= 2) {
      if (grepl("(^|_)NLS(_|$)", text)) {
        return(paste0("SUM-159_NLS_", ploidy_hit[2]))
      }
      return(paste0("SUM-159_", ploidy_hit[2]))
    }

    ploidy_match <- regexpr("SUM-159_[0-9]N", text)
    if (ploidy_match[1] > 0) {
      return(regmatches(text, ploidy_match))
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
    if (grepl("^SUM-159", condition_prefix)) {
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
      if (!(grepl("^SUM-159", condition_prefix) && prefix_cell_line == "SUM-159")) {
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

ordered_condition_groups <- function(groups, group_to_condition) {
  condition_order <- c("Control", "Oxygen", "Glucose", "Phosphate")
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

main <- function() {
  opts <- parse_args(commandArgs(trailingOnly = TRUE))
  df <- read.csv(opts$csv, stringsAsFactors = FALSE, check.names = FALSE)
  annotated <- add_plot_columns(df)

  dir.create(opts$out_dir, showWarnings = FALSE, recursive = TRUE)
  output_path <- file.path(opts$out_dir, opts$output_name)
  annotated_path <- file.path(opts$out_dir, opts$annotated_csv)

  write.csv(annotated, annotated_path, row.names = FALSE)
  plot_grouped_neglog10p(annotated, output_path, opts$width, opts$height, opts$res)

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

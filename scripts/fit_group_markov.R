# Script helpers for grouped Markov fitting. Safe to source from the workflow.

if (!requireNamespace("optparse", quietly = TRUE)) {
  stop("scripts/fit_group_markov.R requires the 'optparse' package.")
}

# Fit the grouped Markov grid and save it to disk.
# Inputs: grouped interval path, output path, and optional parameter grids.
# Output: data.frame of grouped fit results; also written to `output_path`.
fit_group_markov_main <- function(grouped_intervals_path = "core_data/grouped_intervals.Rds",
                                  output_path = "results/group_fit_df.Rds",
                                  pwgd_grid = c(0, 0.001, 0.01, 0.05),
                                  pmis_grid = 10^(seq(-8, -1, length.out = 19)),
                                  fit_mode = "pooled_group_v1",
                                  checkpoint_dir = NULL,
                                  selected_output_path = NULL,
                                  reuse_grid_path = NULL) {
  fit_mode <- match.arg(fit_mode, c("pooled_group_v1", "interval_shared_wgd_v1", "interval_independent_v1"))
  if (missing(output_path)) output_path <- if (fit_mode == "interval_independent_v1") "results/interval_independent/interval_fit_grid.Rds" else if (fit_mode == "interval_shared_wgd_v1") "results/interval_fit_grid.Rds" else "results/group_fit_df.Rds"
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)

  source("R/estimate_pmis.R")
  source("R/interval_fit_utils.R")
  if (fit_mode != "pooled_group_v1") {
    if (is.null(checkpoint_dir)) checkpoint_dir <- paste0(output_path, ".checkpoints")
    if (is.null(selected_output_path)) selected_output_path <- file.path(dirname(output_path), "interval_fits.Rds")
    if (identical(normalizePath(output_path, mustWork = FALSE), normalizePath(selected_output_path, mustWork = FALSE))) {
      stop("Grid and selected fit output paths must differ.")
    }
  }
  if (fit_mode == "pooled_group_v1" && !is.null(reuse_grid_path)) stop("Surface reuse is supported only for explicit interval fit modes.")

  grouped_intervals <- readRDS(grouped_intervals_path)
  group_fit_df <- estimate_group_pmis_grid(
    grouped_intervals,
    pmis_grid = pmis_grid,
    pwgd_grid = pwgd_grid,
    fit_mode = fit_mode,
    checkpoint_dir = checkpoint_dir,
    reuse_grid = if (is.null(reuse_grid_path)) NULL else readRDS(reuse_grid_path)
  )

  nk_fit_atomic_save(group_fit_df, output_path, overwrite = TRUE)
  if (fit_mode != "pooled_group_v1") {
    selected <- nk_select_interval_fits(group_fit_df)
    nk_validate_interval_fits(selected, grouped_intervals)
    nk_fit_atomic_save(selected, selected_output_path, overwrite = TRUE)
    message("Saved ", fit_mode, " fits and profiles to ", selected_output_path)
  }
  message("Saved ", fit_mode, " fit grid to ", output_path)
  group_fit_df
}

# Load cached grouped fit results for the workflow, computing them if absent.
# Inputs: grouped interval path, output path, and optional parameter grids.
# Output: cached grouped fit data.frame.
load_group_markov_results <- function(grouped_intervals_path = "core_data/grouped_intervals.Rds",
                                      output_path = "results/group_fit_df.Rds",
                                      pwgd_grid = c(0, 0.001, 0.01, 0.05),
                                      pmis_grid = 10^(seq(-8, -1, length.out = 19)),
                                      fit_mode = "pooled_group_v1",
                                      checkpoint_dir = NULL,
                                      selected_output_path = NULL,
                                      reuse_grid_path = NULL) {
  fit_mode <- match.arg(fit_mode, c("pooled_group_v1", "interval_shared_wgd_v1", "interval_independent_v1"))
  if (fit_mode != "pooled_group_v1") {
    if (missing(output_path)) output_path <- if (fit_mode == "interval_independent_v1") "results/interval_independent/interval_fit_grid.Rds" else "results/interval_fit_grid.Rds"
    # Revalidate semantic inputs on every invocation. Per-interval checkpoint
    # identities decide reuse; mere existence of an output never does.
    return(fit_group_markov_main(grouped_intervals_path, output_path, pwgd_grid, pmis_grid,
      fit_mode, checkpoint_dir, selected_output_path, reuse_grid_path))
  }
  if (!file.exists(output_path)) {
    fit_group_markov_main(
      grouped_intervals_path = grouped_intervals_path,
      output_path = output_path,
      pwgd_grid = pwgd_grid,
      pmis_grid = pmis_grid,
      fit_mode = fit_mode
    )
  }

  out <- readRDS(output_path)
  if ("fit_mode" %in% names(out) && any(out$fit_mode != "pooled_group_v1")) {
    stop("An interval grid cannot be loaded under pooled compatibility mode.")
  }
  out
}

if (sys.nframe() == 0L) {
  option_list <- list(
    optparse::make_option("--grouped_intervals_path", type = "character", default = "core_data/grouped_intervals.Rds"),
    optparse::make_option("--output_path", type = "character", default = NULL),
    optparse::make_option("--fit_mode", type = "character", default = "pooled_group_v1"),
    optparse::make_option("--checkpoint_dir", type = "character", default = NULL),
    optparse::make_option("--selected_output_path", type = "character", default = NULL),
    optparse::make_option("--reuse_grid_path", type = "character", default = NULL),
    optparse::make_option("--pwgd_grid", type = "character", default = "0,0.001,0.01,0.05"),
    optparse::make_option("--pmis_grid", type = "character", default = NULL)
  )
  opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))

  pwgd_grid <- as.numeric(strsplit(opt$pwgd_grid, ",", fixed = TRUE)[[1]])
  pmis_grid <- if (!is.null(opt$pmis_grid)) {
    as.numeric(strsplit(opt$pmis_grid, ",", fixed = TRUE)[[1]])
  } else {
    10^(seq(-8, -1, length.out = 19))
  }

  invisible(fit_group_markov_main(
    grouped_intervals_path = opt$grouped_intervals_path,
    output_path = if (!is.null(opt$output_path)) opt$output_path else if (opt$fit_mode == "interval_independent_v1") "results/interval_independent/interval_fit_grid.Rds" else if (opt$fit_mode == "interval_shared_wgd_v1") "results/interval_fit_grid.Rds" else "results/group_fit_df.Rds",
    pwgd_grid = pwgd_grid,
    pmis_grid = pmis_grid,
    fit_mode = opt$fit_mode,
    checkpoint_dir = opt$checkpoint_dir,
    selected_output_path = opt$selected_output_path,
    reuse_grid_path = opt$reuse_grid_path
  ))
}

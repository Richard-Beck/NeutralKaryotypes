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
                                  pmis_grid = 10^(seq(-8, -1, length.out = 19))) {
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)

  source("R/estimate_pmis.R")

  grouped_intervals <- readRDS(grouped_intervals_path)
  group_fit_df <- estimate_group_pmis_grid(
    grouped_intervals,
    pmis_grid = pmis_grid,
    pwgd_grid = pwgd_grid
  )

  saveRDS(group_fit_df, output_path)
  message("Saved group fit grid to ", output_path)
  group_fit_df
}

# Load cached grouped fit results for the workflow, computing them if absent.
# Inputs: grouped interval path, output path, and optional parameter grids.
# Output: cached grouped fit data.frame.
load_group_markov_results <- function(grouped_intervals_path = "core_data/grouped_intervals.Rds",
                                      output_path = "results/group_fit_df.Rds",
                                      pwgd_grid = c(0, 0.001, 0.01, 0.05),
                                      pmis_grid = 10^(seq(-8, -1, length.out = 19))) {
  if (!file.exists(output_path)) {
    fit_group_markov_main(
      grouped_intervals_path = grouped_intervals_path,
      output_path = output_path,
      pwgd_grid = pwgd_grid,
      pmis_grid = pmis_grid
    )
  }

  readRDS(output_path)
}

if (sys.nframe() == 0L) {
  option_list <- list(
    optparse::make_option("--grouped_intervals_path", type = "character", default = "core_data/grouped_intervals.Rds"),
    optparse::make_option("--output_path", type = "character", default = "results/group_fit_df.Rds"),
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

  fit_group_markov_main(
    grouped_intervals_path = opt$grouped_intervals_path,
    output_path = opt$output_path,
    pwgd_grid = pwgd_grid,
    pmis_grid = pmis_grid
  )
}

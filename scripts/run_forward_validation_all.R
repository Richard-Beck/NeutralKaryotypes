if (!requireNamespace("optparse", quietly = TRUE)) {
  stop("scripts/run_forward_validation_all.R requires the 'optparse' package.")
}

source("R/workflow_utils.R")
source("scripts/run_forward_validation.R")

# Sanitize group ids for filesystem-safe result names.
# Input: arbitrary group id string.
# Output: filesystem-safe string.
sanitize_id <- function(x) {
  gsub("[^A-Za-z0-9._-]+", "_", x)
}

# Run forward validation for every best-fit group in the grouped fit table.
# Inputs: grouped interval path, group-fit path, output directories, and simulation settings.
# Output: invisibly returns a character vector of full result paths.
run_forward_validation_all_main <- function(grouped_intervals_path = "core_data/grouped_intervals.Rds",
                                            group_fit_path = "results/group_fit_df.Rds",
                                            output_dir = "results",
                                            summary_output_dir = "result_summaries",
                                            bottleneck_size = 1000L,
                                            expansion_factor = 32L,
                                            n_reps = 30L,
                                            record_every = 50L,
                                            n_null_pairs = 100L,
                                            n_cores = default_forward_n_cores(),
                                            seed = 1L) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(summary_output_dir, recursive = TRUE, showWarnings = FALSE)

  group_fit_df <- readRDS(group_fit_path)
  best_group_pmis <- select_best_group_fits(group_fit_df)
  result_paths <- character(nrow(best_group_pmis))

  for (i in seq_len(nrow(best_group_pmis))) {
    forward_group_id <- best_group_pmis$group_id[i]
    stem <- file.path(output_dir, paste0("forward_validation_", sanitize_id(forward_group_id)))
    output_path <- paste0(stem, ".Rds")
    summary_output_path <- file.path(
      summary_output_dir,
      paste0("forward_validation_", sanitize_id(forward_group_id), "_summary.Rds")
    )

    run_forward_validation_main(
      grouped_intervals_path = grouped_intervals_path,
      group_fit_path = group_fit_path,
      forward_group_id = forward_group_id,
      output_path = output_path,
      summary_output_path = summary_output_path,
      bottleneck_size = bottleneck_size,
      expansion_factor = expansion_factor,
      n_reps = n_reps,
      record_every = record_every,
      n_null_pairs = n_null_pairs,
      n_cores = n_cores,
      seed = seed + 1000L * i
    )
    result_paths[i] <- output_path
  }

  invisible(result_paths)
}

if (sys.nframe() == 0L) {
  option_list <- list(
    optparse::make_option("--grouped_intervals_path", type = "character", default = "core_data/grouped_intervals.Rds"),
    optparse::make_option("--group_fit_path", type = "character", default = "results/group_fit_df.Rds"),
    optparse::make_option("--output_dir", type = "character", default = "results"),
    optparse::make_option("--summary_output_dir", type = "character", default = "result_summaries"),
    optparse::make_option("--bottleneck_size", type = "integer", default = 1000L),
    optparse::make_option("--expansion_factor", type = "integer", default = 32L),
    optparse::make_option("--n_reps", type = "integer", default = 30L),
    optparse::make_option("--record_every", type = "integer", default = 50L),
    optparse::make_option("--n_null_pairs", type = "integer", default = 100L),
    optparse::make_option("--n_cores", type = "integer", default = default_forward_n_cores()),
    optparse::make_option("--seed", type = "integer", default = 1L)
  )
  opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))

  run_forward_validation_all_main(
    grouped_intervals_path = opt$grouped_intervals_path,
    group_fit_path = opt$group_fit_path,
    output_dir = opt$output_dir,
    summary_output_dir = opt$summary_output_dir,
    bottleneck_size = opt$bottleneck_size,
    expansion_factor = opt$expansion_factor,
    n_reps = opt$n_reps,
    record_every = opt$record_every,
    n_null_pairs = opt$n_null_pairs,
    n_cores = opt$n_cores,
    seed = opt$seed
  )
}

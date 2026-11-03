if (!requireNamespace("optparse", quietly = TRUE)) {
  stop("scripts/run_forward_validation_all.R requires the 'optparse' package.")
}

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
                                            group_fit_path = forward_validation_default_paths(fit_mode)$group_fit_path,
                                            output_dir = forward_validation_default_paths(fit_mode)$output_dir,
                                            summary_output_dir = forward_validation_default_paths(fit_mode)$summary_output_dir,
                                            bottleneck_size = 1000L,
                                            expansion_factor = 32L,
                                            n_reps = 30L,
                                            record_every = 50L,
                                            n_null_pairs = 100L,
                                            n_cores = default_forward_n_cores(),
                                            seed = 1L,
                                            distance_metric = c("chrom_weighted_wasserstein", "wasserstein"),
                                            fit_mode = c("pooled_group_v1", "interval_shared_wgd_v1", "interval_independent_v1"),
                                            checkpoint_dir = file.path(output_dir, "forward_checkpoints_interval_v1"),
                                            checkpoint_batch_size = n_cores,
                                            aggregate_output_path = forward_validation_default_paths(fit_mode)$aggregate_output_path) {
  fit_mode <- match.arg(fit_mode)
  distance_metric <- match.arg(distance_metric)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(summary_output_dir, recursive = TRUE, showWarnings = FALSE)

  group_fit_df <- readRDS(group_fit_path)
  best_group_pmis <- if (fit_mode == "pooled_group_v1") select_best_group_fits(group_fit_df) else {
    select_forward_fits(group_fit_df, readRDS(grouped_intervals_path), fit_mode)
  }
  group_ids <- unique(as.character(best_group_pmis$group_id))
  if (fit_mode != "pooled_group_v1") group_ids <- sort(group_ids)
  if (anyDuplicated(sanitize_id(group_ids))) stop("Group IDs collide in forward output filenames.")
  result_paths <- character(length(group_ids))

  for (i in seq_along(group_ids)) {
    forward_group_id <- group_ids[i]
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
      seed = if (fit_mode == "pooled_group_v1") seed + 1000L * i else seed,
      distance_metric = distance_metric, fit_mode = fit_mode,
      checkpoint_dir = if (fit_mode != "pooled_group_v1") checkpoint_dir else NULL,
      checkpoint_batch_size = checkpoint_batch_size, update_aggregate = FALSE
    )
    result_paths[i] <- output_path
  }

  write_neutral_probability_csv(
    summary_dir = summary_output_dir,
    grouped_intervals_path = grouped_intervals_path,
    output_path = aggregate_output_path, fit_mode = fit_mode
  )

  invisible(result_paths)
}

if (sys.nframe() == 0L) {
  option_list <- list(
    optparse::make_option("--grouped_intervals_path", type = "character", default = "core_data/grouped_intervals.Rds"),
    optparse::make_option("--group_fit_path", type = "character", default = NULL),
    optparse::make_option("--fit_mode", type = "character", default = "pooled_group_v1"),
    optparse::make_option("--checkpoint_dir", type = "character", default = NULL),
    optparse::make_option("--checkpoint_batch_size", type = "integer", default = NULL),
    optparse::make_option("--aggregate_output_path", type = "character", default = NULL),
    optparse::make_option("--output_dir", type = "character", default = NULL),
    optparse::make_option("--summary_output_dir", type = "character", default = NULL),
    optparse::make_option("--bottleneck_size", type = "integer", default = 1000L),
    optparse::make_option("--expansion_factor", type = "integer", default = 32L),
    optparse::make_option("--n_reps", type = "integer", default = 30L),
    optparse::make_option("--record_every", type = "integer", default = 50L),
    optparse::make_option("--n_null_pairs", type = "integer", default = 100L),
    optparse::make_option("--n_cores", type = "integer", default = default_forward_n_cores()),
    optparse::make_option("--seed", type = "integer", default = 1L),
    optparse::make_option(
      "--distance_metric",
      type = "character",
      default = "chrom_weighted_wasserstein",
      help = "Posterior-predictive endpoint distance: chrom_weighted_wasserstein or wasserstein."
    )
  )
  opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))
  default_paths <- forward_validation_default_paths(opt$fit_mode)
  for (field in c("group_fit_path", "output_dir", "summary_output_dir", "aggregate_output_path")) {
    if (is.null(opt[[field]])) opt[[field]] <- default_paths[[field]]
  }

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
    seed = opt$seed,
    distance_metric = opt$distance_metric, fit_mode = opt$fit_mode,
    checkpoint_dir = if (is.null(opt$checkpoint_dir)) file.path(opt$output_dir, "forward_checkpoints_interval_v1") else opt$checkpoint_dir,
    checkpoint_batch_size = if (is.null(opt$checkpoint_batch_size)) opt$n_cores else opt$checkpoint_batch_size,
    aggregate_output_path = opt$aggregate_output_path
  )
}

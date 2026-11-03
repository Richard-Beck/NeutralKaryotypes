source("R/estimate_pmis.R")
source("R/interval_fit_utils.R")

fails <- function(expr) stopifnot(inherits(tryCatch(force(expr), error = identity), "error"))
same <- function(x, y, tolerance = 1e-12) stopifnot(isTRUE(all.equal(x, y, tolerance = tolerance, check.attributes = FALSE)))
make_interval <- function(id, end = c(2, 2, 2, 2), divisions = 4L) {
  list(connected_tree = "native_fixture", start = paste0("start_", id), end = paste0("end_", id),
       start_ids = paste0("sample_start_", id), end_ids = paste0("sample_end_", id),
       start_alias = NULL, n_collapsed_passages = 0L, condition = "fixture",
       elapsed_passages = 1L, conditions = c(fixture = 1L),
       start_karyotype = matrix(2, nrow = 4, ncol = 22),
       end_karyotype = matrix(rep(end, 22), nrow = length(end), ncol = 22), n_divisions = divisions)
}
groups <- list(g = list(make_interval("a"), make_interval("b", c(1, 3, 4, 4))),
               h = list(make_interval("c", c(2, 2, 3, 2))))
p <- c(0.001, 0.02, 0.08)
q <- c(0, 0.01, 0.05)
scratch <- tempfile("interval-fitting-")
dir.create(scratch)
grid <- nk_interval_fit_grid(groups, p, q, scratch)
fits <- nk_select_interval_fits(grid)
stopifnot(nrow(fits) == 3L, nrow(grid) == 27L, all(fits$uncertainty_status == "not_calibrated"),
          all(attr(grid, "interval_checkpoints")$cache_status == "computed"),
          all(attr(fits, "group_checkpoints")$cache_status == "computed"))

# The maintained interval kernel is exactly the old pooled objective's term.
old <- estimate_group_pmis_grid(groups, p, q)
for (j in seq_len(nrow(old))) {
  part <- grid[grid$group_id == old$group_id[j] & grid$pmis == old$pmis[j] & grid$pwgd == old$pwgd[j], ]
  same(sum(part$negll), old$negll[j])
}
singleton_old <- old[old$group_id == "h", ]
singleton_old <- singleton_old[order(singleton_old$negll, singleton_old$pwgd, singleton_old$pmis), ]
same(fits[fits$group_id == "h", c("pmis", "pwgd", "negll")], singleton_old[1, c("pmis", "pwgd", "negll")])

# Independent exhaustive reference over p1, p2 and shared q (not reducer code).
exhaustive <- expand.grid(p1 = p, p2 = p, q = q)
exhaustive$nll <- vapply(seq_len(nrow(exhaustive)), function(i) {
  row <- exhaustive[i, ]
  a <- groups$g[[1]]
  b <- groups$g[[2]]
  get_interval_nll_markov(row$p1, row$q, a$start_karyotype, a$end_karyotype, a$n_divisions) +
    get_interval_nll_markov(row$p2, row$q, b$start_karyotype, b$end_karyotype, b$n_divisions)
}, numeric(1))
same(unique(fits$group_negll[fits$group_id == "g"]), min(exhaustive$nll))
stopifnot(sum(fits$negll[fits$group_id == "g"]) <= min(old$negll[old$group_id == "g"]),
          length(unique(fits$pmis[fits$group_id == "g"])) == 2L)
same(attr(fits, "profiles")$g$pwgd_profiles$negll,
     vapply(sort(q), function(value) min(exhaustive$nll[exhaustive$q == value]), numeric(1)))

# Interval order and grid traversal cannot alter parameter or fit identity.
reordered <- list(h = groups$h, g = rev(groups$g))
again_grid <- nk_interval_fit_grid(reordered, rev(p), rev(q), scratch)
again <- nk_select_interval_fits(again_grid[nrow(again_grid):1L, ])
m <- match(fits$physical_interval_id, again$physical_interval_id)
same(fits[, c("pmis", "pwgd", "negll", "fit_id")], again[m, c("pmis", "pwgd", "negll", "fit_id")])
stopifnot(all(attr(again_grid, "interval_checkpoints")$cache_status == "reused"),
          all(attr(again, "group_checkpoints")$cache_status == "reused"))
validated <- nk_validate_interval_fits(fits, reordered)
stopifnot(identical(validated$interval_idx, nk_interval_roster(reordered)$interval_idx[match(validated$physical_interval_id, nk_interval_roster(reordered)$physical_interval_id)]))
same(nk_select_interval_fits(fits), fits)

# Neighbour updates invalidate the shared reduction, not unrelated surfaces.
changed <- groups
changed$g[[1]]$end_karyotype[1, 1] <- 3
changed_grid <- nk_interval_fit_grid(changed, p, q, scratch)
changed_fits <- nk_select_interval_fits(changed_grid)
stopifnot(sum(attr(changed_grid, "interval_checkpoints")$cache_status == "computed") == 1L,
          identical(attr(changed_fits, "group_checkpoints")$cache_status, c("computed", "reused")))
stopifnot(all(changed_fits$group_profile_id[changed_fits$group_id == "g"] != fits$group_profile_id[fits$group_id == "g"]))
same(changed_fits[changed_fits$group_id == "h", c("pmis", "pwgd", "fit_id")], fits[fits$group_id == "h", c("pmis", "pwgd", "fit_id")])

# Native identity excludes numeric division convention and incidental labels.
labelled <- groups
labelled$g[[1]]$display_label <- "updated manuscript name"
stopifnot(all(attr(nk_interval_fit_grid(labelled, p, q, scratch), "interval_checkpoints")$cache_status == "reused"))
divided <- groups
divided$g[[1]]$n_divisions <- 5L
same(nk_physical_interval_id(divided$g[[1]]), nk_physical_interval_id(groups$g[[1]]))
stopifnot(sum(attr(nk_interval_fit_grid(divided, p, q, scratch), "interval_checkpoints")$cache_status == "computed") == 1L)

# Exact ties stay visible and use the declared lowest-q/lowest-p policy.
tie_groups <- list(t = list(make_interval("tie", divisions = 0L)))
tie_grid <- nk_interval_fit_grid(tie_groups, rev(p), rev(q))
tie <- nk_select_interval_fits(tie_grid)
stopifnot(tie$pmis == min(p), tie$pwgd == min(q), tie$pmis_tie_count == length(p),
          tie$pwgd_tie_count == length(q), tie$pmis_flat, tie$pwgd_flat,
          tie$estimability_status == "tied_grid_optimum")
tie_old <- estimate_group_pmis_grid(tie_groups, p, q)
same(tie$negll, min(tie_old$negll))
no_support <- tie_grid
no_support$negll <- Inf
attr(no_support, "grid_sha256") <- NULL
unavailable <- nk_select_interval_fits(no_support)
stopifnot(unavailable$estimability_status == "no_finite_shared_wgd_profile", is.na(unavailable$pmis), is.na(unavailable$pwgd))

# Invalid identities/rosters/grids fail; positional or display matches never fill gaps.
fails(nk_interval_roster(list(g = c(groups$g, groups$g[1]))))
bad <- groups
bad$g[[1]]$start_ids <- NULL
fails(nk_interval_roster(bad))
bad <- groups
bad$g[[1]]$comparison_type <- "explicit_full_window"
fails(nk_interval_roster(bad))
bad <- groups
bad$g[[2]] <- groups$g[[1]]
bad$g[[2]]$comparison_type <- "explicit_full_window"
bad$g[[2]]$model_path_passage_ids <- c(bad$g[[2]]$start, "exact_middle", bad$g[[2]]$end)
fails(nk_interval_roster(bad))
fails(nk_select_interval_fits(old))
fails(nk_interval_fit_grid(groups, c(p, p[1]), q))
fails(nk_select_interval_fits(grid[-1L, ]))
fails(nk_validate_interval_fits(fits[-1L, ], groups))
fails(nk_validate_interval_fits(fits, changed))
bad <- fits
bad$pmis[1] <- 0.007
fails(nk_validate_interval_fits(bad))

# Corrupt/stale payloads are never silently reused or replaced.
checkpoint_path <- attr(grid, "interval_checkpoints")$checkpoint_path[1]
saved <- readRDS(checkpoint_path)
bad <- saved
bad$data$negll[1] <- bad$data$negll[1] + 1
saveRDS(bad, checkpoint_path)
fails(nk_interval_fit_grid(groups, p, q, scratch))
saveRDS(saved, checkpoint_path)
bad <- saved
bad$identity$schema <- "old_schema"
saveRDS(bad, checkpoint_path)
fails(nk_interval_fit_grid(groups, p, q, scratch))
saveRDS(saved, checkpoint_path)

# Runtime and kernel changes invalidate all dependent likelihood surfaces.
original_runtime <- nk_interval_runtime_identity
nk_interval_runtime_identity <- function() c(original_runtime(), list(test_runtime_change = TRUE))
stopifnot(all(attr(nk_interval_fit_grid(groups, p, q, scratch), "interval_checkpoints")$cache_status == "computed"))
nk_interval_runtime_identity <- original_runtime
original_kernel <- get_interval_nll_markov
get_interval_nll_markov <- function(...) original_kernel(...) + 0
stopifnot(all(attr(nk_interval_fit_grid(groups, p, q, scratch), "interval_checkpoints")$cache_status == "computed"))
get_interval_nll_markov <- original_kernel

# Maintained entrypoint saves full surfaces and selected profiles separately.
source("scripts/fit_group_markov.R")
input_path <- file.path(scratch, "inputs.Rds")
output_path <- file.path(scratch, "interval_fit_grid.Rds")
selected_path <- file.path(scratch, "interval_fits.Rds")
saveRDS(groups, input_path)
entry <- fit_group_markov_main(input_path, output_path, q, p, "interval_shared_wgd_v1", scratch, selected_path)
same(readRDS(selected_path)$fit_id, fits$fit_id)
stopifnot(nrow(readRDS(output_path)) == nrow(grid), length(attr(readRDS(selected_path), "profiles")) == 2L)
fails(load_group_markov_results(input_path, output_path, q, p, "pooled_group_v1"))
unlink(scratch, recursive = TRUE)
message("Interval kernel parity, exhaustive shared-WGD profiling, identities, ties, caches and entrypoint checks passed.")

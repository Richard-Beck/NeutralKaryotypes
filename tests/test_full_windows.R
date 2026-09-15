source("R/db_utils.R")

fails <- function(expr) stopifnot(inherits(tryCatch(force(expr), error = identity), "error"))
tree <- list(
  sampled_passage_ids = c("proxy", "mid", "late"),
  tree_edges = data.frame(passage_id = c("root", "proxy", "early", "mid", "late_pre", "late"),
                         passage_from = c(NA, "root", "root", "early", "mid", "late_pre")),
  tree_conditions = data.frame(passage_id = c("root", "proxy", "early", "mid", "late_pre", "late"),
                              condition = c("control", "control", "glucose", "glucose", "glucose", "glucose")),
  tree_samples = data.frame(passage_id = c("proxy", "mid", "late"), samples = c("p", "m", "l"))
)
original <- extract_karyotyped_intervals(tree)
stopifnot(length(original) == 2L)
request <- data.frame(connected_tree = "fixture", start = "proxy", end = "late", start_alias = "root")
expanded <- extract_karyotyped_intervals(tree, request)
stopifnot(identical(original, expanded[1:2]), length(expanded) == 3L,
          identical(expanded[[3]]$model_path_passage_ids, c("proxy", "root", "early", "mid", "late_pre", "late")),
          expanded[[3]]$n_collapsed_passages == 1L,
          sum(expanded[[3]]$conditions) == 6L,
          sum(vapply(original, function(x) sum(x$conditions), numeric(1))) == 7L)
bad <- request
bad$start_alias <- "early"
fails(extract_karyotyped_intervals(tree, bad))
bad$start_alias <- NA_character_
fails(extract_karyotyped_intervals(tree, bad))
bad <- request
bad$end <- "missing"
fails(extract_karyotyped_intervals(tree, bad))
fails(extract_karyotyped_intervals(tree, rbind(request, request)))
bad <- request
bad$end <- "mid"
fails(extract_karyotyped_intervals(tree, bad))
bad <- request
bad$connected_tree <- "unknown"
fails(build_simulation_intervals(list(fixture = tree), full_windows = bad))
direct <- tree
direct$sampled_passage_ids <- c("root", "mid", "late")
direct_request <- data.frame(start = "root", end = "late", start_alias = NA_character_)
stopifnot(length(extract_karyotyped_intervals(direct, direct_request)) == 3L)
direct$sampled_passage_ids <- c("early", "mid", "late")
fails(extract_karyotyped_intervals(direct,
  data.frame(start = "early", end = "late", start_alias = "root")))
karyotypes <- data.frame(id = c("p", "m", "l"))
karyotypes$karyotype <- I(rep(list(rep(2, 22)), 3))
base <- build_simulation_intervals(list(fixture = tree), karyotypes)
out <- build_simulation_intervals(list(fixture = tree), karyotypes, request)
stopifnot(identical(base, out[1:2]), out[[3]]$elapsed_passages == 6L,
          identical(dim(out[[3]]$start_karyotype), c(1L, 22L)),
          group_simulation_intervals(out)[[1]][[3]]$n_divisions == 30L)
message("Full-window default parity, proxy/path arithmetic and negative checks passed.")

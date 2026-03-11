edge_df_to_phylo <- function(df, id_col = "passage_id", parent_col = "passage_from", edge_length = NULL) {
  library(dplyr)
  library(ape)

  df <- df %>%
    dplyr::transmute(id = .data[[id_col]], parent = .data[[parent_col]]) %>%
    dplyr::distinct()

  stopifnot(!anyDuplicated(df$id))
  stopifnot(!any(is.na(df$id)))

  roots <- df$id[is.na(df$parent) | !(df$parent %in% df$id)]
  if (length(roots) != 1) stop("Expected exactly one root; found: ", paste(roots, collapse = ", "))

  child_map <- split(df$id[!is.na(df$parent)], df$parent[!is.na(df$parent)])
  child_map <- lapply(child_map, unique)

  quote_label <- function(x) {
    paste0("'", gsub("'", "", x, fixed = TRUE), "'")
  }

  build_newick <- function(node) {
    kids <- child_map[[node]]
    if (is.null(kids) || !length(kids)) {
      return(quote_label(node))
    }

    kids <- sort(kids)
    paste0(
      "(",
      paste(vapply(kids, build_newick, character(1)), collapse = ","),
      ")",
      quote_label(node)
    )
  }

  phy <- ape::read.tree(text = paste0(build_newick(roots[1]), ";"))
  phy$tip.label <- gsub("^'|'$", "", phy$tip.label)
  phy$tip.label <- gsub(" ", "_", phy$tip.label, fixed = TRUE)
  if (!is.null(phy$node.label)) {
    phy$node.label <- gsub("^'|'$", "", phy$node.label)
    phy$node.label <- gsub(" ", "_", phy$node.label, fixed = TRUE)
  }
  if (!is.null(edge_length)) {
    el <- edge_length(df)
    if (length(el) == nrow(phy$edge)) {
      phy$edge.length <- as.numeric(el)
    }
  }
  phy
}

plot_passage_tree <- function(tree, component_name = NULL) {
  library(dplyr)
  library(ggplot2)
  library(ggtree)

  edge_df <- tree$tree_edges
  phy <- edge_df_to_phylo(edge_df)

  sampled_passages <- unique(tree$sampled_passage_ids)
  sample_counts <- table(tree$tree_samples$passage_id)
  node_ids <- unique(c(as.character(edge_df$passage_id), as.character(edge_df$passage_from)))
  node_ids <- node_ids[!is.na(node_ids)]

  meta <- data.frame(
    label = node_ids,
    stringsAsFactors = FALSE
  ) %>%
    dplyr::distinct() %>%
    dplyr::mutate(
      karyotyped = label %in% sampled_passages,
      sample_count = as.integer(sample_counts[label])
    )

  meta$sample_count[is.na(meta$sample_count)] <- 0L

  p <- ggtree(phy, layout = "circular", size = 0.35)
  p$data <- dplyr::left_join(p$data, meta, by = "label")
  p$data$karyotyped[is.na(p$data$karyotyped)] <- FALSE
  p$data$sample_count[is.na(p$data$sample_count)] <- 0L
  p$data$highlight <- !is.na(p$data$label) & p$data$label %in% sampled_passages

  p <- p +
    geom_point2(
      data = ~ dplyr::filter(.x, !is.na(label)),
      aes(size = pmax(sample_count, 1L)),
      shape = 16,
      color = "grey70",
      alpha = 0.9
    ) +
    geom_point2(
      data = ~ dplyr::filter(.x, highlight),
      aes(shape = "karyotyped"),
      color = "#c83e4d",
      size = 3.4
    ) +
    scale_size_continuous(
      range = c(1.3, 3.2),
      breaks = c(1, 2, 3, 5),
      name = "Samples"
    ) +
    scale_shape_discrete("") +
    theme_void() +
    ggtitle(if (is.null(component_name)) "Connected passage tree" else paste("Connected passage tree:", component_name))

  p
}

experiment_lineage_plot <- function(lineage_df_subset){
  library(dplyr)
  library(ggplot2)
  library(ggtree)

  phy <- edge_df_to_phylo(lineage_df_subset)

  meta <- lineage_df_subset %>%
    dplyr::transmute(
      label = as.character(passage_id),
      g = ifelse(is.finite(g), g, NA_real_),
      karyotyped = as.logical(karyotyped),
      has_flow = as.logical(has_flow)
    )

  p <- ggtree(phy, layout="circular", size=0.3)
  p$data <- dplyr::left_join(p$data, meta, by="label")

  p <- p +
    geom_point2(
      data = ~ dplyr::filter(.x, karyotyped),
      aes(shape = "karyotyped"), color = "red", size = 3.2
    ) +
    geom_point2(
      data = ~ dplyr::filter(.x, has_flow),
      aes(shape = "has_flow"), color = "green", size = 3.2
    ) +
    geom_point2(
      aes(color = g),
      shape = 16, size = 1.4
    ) +
    scale_color_viridis_c("growth\nrate", na.value = "grey70", trans = "log") +
    theme_void() +
    scale_shape_discrete("") +
    ggtitle(paste("Hypoxia experiment data:", unique(lineage_df_subset$label_value)))
  return(p)
}

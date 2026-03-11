load_db_vars <- function(filepath) {
  if (!file.exists(filepath)) {
    stop(sprintf("Error: File '%s' does not exist.", filepath))
  }
  
  lines <- tryCatch(readLines(filepath, warn = FALSE), 
                    error = function(e) stop("Error: Failed to read file '", filepath, "'. ", e$message))
  
  if (length(lines) == 0) {
    stop(sprintf("Error: File '%s' is empty.", filepath))
  }
  
  if (any(!grepl("^[A-Za-z_][A-Za-z0-9_]*=.+$", lines))) {
    stop("Error: All lines must be in the format KEY=value, with valid variable names.")
  }
  
  kv_pairs <- strsplit(lines, "=", fixed = TRUE)
  keys <- vapply(kv_pairs, `[`, "", 1)
  values <- vapply(kv_pairs, `[`, "", 2)
  vars <- setNames(values, keys)
  
  required_keys <- c("HOST", "DBNAME", "USER", "PASSWORD")
  missing_keys <- setdiff(required_keys, names(vars))
  if (length(missing_keys) > 0) {
    stop("Error: Missing required keys: ", paste(missing_keys, collapse = ", "))
  }
  
  return(vars)
}

get_karyotyping <- function(cloneIds){
  require(DBI)
  require(RMariaDB)
  dbvars <- load_db_vars("db_creds.txt")
  db <- dbConnect(
    MariaDB(),
    host = dbvars["HOST"],
    user = dbvars["USER"],
    password = dbvars["PASSWORD"],
    dbname = dbvars["DBNAME"]
  )
  
  origin_clause <- paste0("'", cloneIds, "'", collapse = ", ")
  q <- paste0("SELECT * FROM Perspective WHERE origin IN (", origin_clause, ") AND whichPerspective='GenomePerspective'")
  
  rs <- dbSendQuery(db, q)
  res <- dbFetch(rs)
  dbClearResult(rs) 
  dbDisconnect(db)
  
  kvecs <- lapply(res$profile, function(p) {
    raw_vec <- p
    readBin(raw_vec, what = "double", n = length(raw_vec) / 8, endian = "big")
  })
  
  df <- tibble::tibble(
    id = res$origin,
    karyotype = kvecs
  )
  df
}

collapse_db <- function(x) {
  required_cols <- c("id", "event", "passaged_from_id1")
  missing_cols <- setdiff(required_cols, names(x))
  if (length(missing_cols) > 0) {
    stop("collapse_db() missing required columns: ",
         paste(missing_cols, collapse = ", "))
  }

  if (anyDuplicated(x$id)) {
    stop("collapse_db() requires unique values in x$id.")
  }

  x <- as.data.frame(x, stringsAsFactors = FALSE)
  x$id <- as.character(x$id)
  x$event <- as.character(x$event)
  x$passaged_from_id1 <- as.character(x$passaged_from_id1)
  x$passaged_from_id1[is.na(x$passaged_from_id1) | x$passaged_from_id1 == ""] <- NA_character_

  valid_events <- c("seeding", "harvest")
  bad_event_rows <- !(x$event %in% valid_events) | is.na(x$event)
  if (any(bad_event_rows)) {
    warning(
      "collapse_db() stripped ", sum(bad_event_rows),
      " rows with unsupported event values."
    )
    x <- x[!bad_event_rows, , drop = FALSE]
  }

  missing_parent_rows <- !is.na(x$passaged_from_id1) & !(x$passaged_from_id1 %in% x$id)
  if (any(missing_parent_rows)) {
    warning(
      "collapse_db() stripped ", sum(missing_parent_rows),
      " rows whose parent id was absent from x$id."
    )
    x <- x[!missing_parent_rows, , drop = FALSE]
  }

  root_rows <- is.na(x$passaged_from_id1)
  is_passage_row <- x$event == "seeding" | root_rows
  valid_passage_ids <- x$id[is_passage_row]
  parent_lookup <- setNames(x$passaged_from_id1, x$id)

  resolve_passage_id <- function(id, parent_map, valid_ids) {
    cur <- id
    visited <- character()
    while (!is.null(cur) && !is.na(cur)) {
      if (cur %in% visited) {
        stop("collapse_db() detected a cycle while resolving passage ids.")
      }
      if (cur %in% valid_ids) return(cur)
      visited <- c(visited, cur)
      cur <- unname(parent_map[cur])
    }
    NA_character_
  }

  sample_passage_id <- vapply(
    x$id,
    resolve_passage_id,
    character(1),
    parent_map = parent_lookup,
    valid_ids = valid_passage_ids
  )

  unresolved_rows <- is.na(sample_passage_id)
  if (any(unresolved_rows)) {
    warning(
      "collapse_db() stripped ", sum(unresolved_rows),
      " rows that could not be resolved to any upstream passage_id."
    )
    x <- x[!unresolved_rows, , drop = FALSE]
    sample_passage_id <- sample_passage_id[!unresolved_rows]
    root_rows <- root_rows[!unresolved_rows]
    is_passage_row <- is_passage_row[!unresolved_rows]
    valid_passage_ids <- valid_passage_ids[valid_passage_ids %in% x$id]
    parent_lookup <- setNames(x$passaged_from_id1, x$id)
  }

  samples <- data.frame(
    passage_id = sample_passage_id,
    samples = x$id,
    stringsAsFactors = FALSE
  )

  passaging_rows <- x[is_passage_row, , drop = FALSE]
  parent_passage_id <- vapply(
    passaging_rows$passaged_from_id1,
    resolve_passage_id,
    character(1),
    parent_map = parent_lookup,
    valid_ids = valid_passage_ids
  )
  parent_passage_id[is.na(passaging_rows$passaged_from_id1)] <- NA_character_

  passaging <- data.frame(
    passage_id = passaging_rows$id,
    passage_from = unname(parent_passage_id),
    stringsAsFactors = FALSE
  )

  list(
    passaging = passaging,
    samples = samples
  )
}

recover_lineage <- function(row, data) {
  current_id <- row$last
  lineage <- c(current_id)
  
  while (!is.na(current_id) && current_id != row$first) {
    parent <- data$passaged_from_id1[which(data$id == current_id)]
    if (length(parent) == 0 || is.na(parent)) {
      stop(paste("Could not trace from", row$last, "to", row$first))
    }
    lineage <- c(parent, lineage)
    current_id <- parent
  }
  
  if (current_id != row$first) {
    stop(paste("Could not find", row$first, "starting from", row$last))
  }
  
  seeding_lineage <- lineage[data$event[match(lineage, data$id)] == "seeding"]
  
  filtered_x <- do.call(rbind, lapply(seq_along(seeding_lineage), function(i) {
    dfi <- data[data$id %in% seeding_lineage[i] | data$passaged_from_id1 %in% seeding_lineage[i], ]
    dfi$adjPass <- stringr::str_pad(i, width = 2)
    dfi$date <- as.Date(dfi$date)
    dfi$num_date <- as.numeric(as.Date(dfi$date))
    dfi$num_date <- dfi$num_date - min(dfi$num_date)
    dfi$intercept <- NaN
    dfi$g <- NaN
    if (nrow(dfi) < 2) return(dfi)
    if(sum(!is.na(dfi$correctedCount))<2) return(dfi)
    fit <- lm(log(pmax(1, dfi$correctedCount)) ~ dfi$num_date)
    dfi$intercept <- exp(coef(fit)[1])
    dfi$g <- coef(fit)[2]
    dfi
  }))
  
  filtered_x$label_value    <- row$label
  filtered_x$sublabel_value <- row$label2
  return(filtered_x)
}


fill_lineage_gaps <- function(ploidy_substr,df,x){
  g <- x[x$cellLine=="SUM-159",c("id","passaged_from_id1")]
  
  
  ## perhaps split 2N 4N.
  karyotyped_lineages <- lapply(cloneIds[grepl(ploidy_substr,cloneIds)],function(id){
    ids <- c()
    while(!is.na(id)){
      ids <- c(id,ids)
      id <- g$passaged_from_id1[g$id==id]
    }
    ids
  })
  
  init_ids <- sapply(filters,'[[',"first")
  init_ids <- init_ids[grepl(ploidy_substr,init_ids)]
  
  imaged_lineages <- lapply(init_ids,function(id){
    ids <- c()
    while(!is.na(id)){
      ids <- c(id,ids)
      id <- g$passaged_from_id1[g$id==id]
    }
    ids
  })
  
  all_lineages <- c(imaged_lineages,karyotyped_lineages)
  n_lineages <- length(all_lineages)
  
  uids <- table(unlist(all_lineages))
  n_remove <- length(uids[uids==n_lineages])-1
  
  all_lineages <- lapply(all_lineages, function(i){
    i[-c(1:n_remove)]
  })
  
  x <- x[x$id%in%unlist(all_lineages),]
  x <- x[!x$id%in%df$id,]
  x$label_value <- gsub("_","",ploidy_substr)
  x
}

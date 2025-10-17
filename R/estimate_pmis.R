states <- 0:8
eps <- 1e-9


library(Rcpp)

# Compile the C++ file once (creates run_karyotype_neutral)
Rcpp::sourceCpp("ksim2.cpp")

# From a cell-by-22 integer matrix -> per-chromosome histograms (rows=states, cols=chr)
hist_mat_from_matrix <- function(M){
  do.call(cbind, lapply(1:ncol(M), function(j){
    tab <- table(factor(M[,j], levels=states))
    as.numeric(tab)
  }))
}

hist_mat_from_named_counts <- function(tbl, nchr = 22, states = 0:8, sep = "\\.") {
  nm <- names(tbl); if (is.null(nm)) stop("Input must have names (karyotype strings).")
  parts <- strsplit(nm, sep)
  stopifnot(all(lengths(parts) == nchr))
  M <- do.call(rbind, lapply(parts, function(x) as.integer(x)))  # (#types x nchr)
  w <- as.numeric(tbl)                                           # weights per type
  
  H <- matrix(0, nrow = length(states), ncol = nchr,
              dimnames = list(states = states, chr = paste0("chr", seq_len(nchr))))
  for (j in seq_len(nchr)) {
    f <- factor(M[, j], levels = states)
    v <- tapply(w, f, sum)                                       # NA for empty levels
    H[, j] <- ifelse(is.na(v), 0, v)
  }
  H
}

matrix_to_named_counts <- function(M, sep=".") {
  nm <- apply(M, 1, paste, collapse=sep)
  as.numeric(table(factor(nm)))
  # but keep names
  tbl <- table(nm)
  tbl
}

gen_pop <- function(H,h=0.1,states=0:8,n_init=10000){
  p0 <- do.call(rbind,lapply(states,function(j){
    sapply(1:22,function(cr){
      sum(exp(-(j-H[,cr])^2/h^2))
    })
  }))
  for(i in 1:ncol(p0)) p0[,i] <- p0[,i]/sum(p0[,i])
  
  pop0 <- do.call(cbind,lapply(1:ncol(H),function(i){
    sample(states,n_init,replace = T,prob=p0[,i])
  }))  
}

# NLL of observed counts H_obs under simulated probs P_sim (columns = chr)
nll_multinomial <- function(H_obs, P_sim, eps=1e-9){
  P <- pmax(P_sim, eps); P <- sweep(P, 2, colSums(P), "/")
  -sum(H_obs * log(P))
}

get_nll <- function(p_mis,K0,H_obs,n_steps,max_pop,cull_keep,seed, record_every=5){
  pop0 <- matrix_to_named_counts(gen_pop(round(K0)))
  res <- run_karyotype_neutral(
    initial_counts_named = pop0,
    rate = 1.0,
    p_misseg = p_mis,
    dt = 0.1,
    n_steps = n_steps,
    max_pop = max_pop,
    cull_keep = cull_keep,
    record_every = record_every,
    seed = seed
  )
  last_tbl <- tail(res,1)[[1]]
  H_sim <- hist_mat_from_named_counts(last_tbl, nchr=22)
  P_sim <- sweep(H_sim, 2, colSums(H_sim), "/")
  nll_multinomial(H_obs, P_sim, eps)
}

get_pop <- function(p_mis,K0,nsteps,max_pop,cull_keep,seed, record_every=5){
  pop0 <- matrix_to_named_counts(gen_pop(round(K0)))
  res <- run_karyotype_neutral(
    initial_counts_named = pop0,
    rate = 1.0,
    p_misseg = p_mis,
    dt = 0.1,
    n_steps = n_steps,
    max_pop = max_pop,
    cull_keep = cull_keep,
    record_every = record_every,
    seed = seed
  )
  last_tbl <- tail(res,1)[[1]]
  last_tbl
}
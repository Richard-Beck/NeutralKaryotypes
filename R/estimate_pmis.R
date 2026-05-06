states <- 0:8
eps <- 1e-9


library(Rcpp)

# Compile the C++ simulator on first use.
# Input: none.
# Output: invisibly returns TRUE after `run_karyotype_neutral()` is available.
ensure_ksim_compiled <- local({
  compiled <- FALSE
  function() {
    if (!compiled || !exists("run_karyotype_neutral", mode = "function")) {
      Rcpp::sourceCpp("ksim2.cpp")
      compiled <<- TRUE
    }
    invisible(TRUE)
  }
})

# Convert a cell-by-chromosome matrix into copy-number histograms.
# Inputs: `M` numeric matrix, `states` allowed copy-number states.
# Output: matrix with rows = states and columns = chromosomes.
hist_mat_from_matrix <- function(M, states = 0:8){
  do.call(cbind, lapply(1:ncol(M), function(j){
    tab <- table(factor(M[,j], levels=states))
    as.numeric(tab)
  }))
}

# Convert named whole-karyotype counts into per-chromosome histograms.
# Inputs: `tbl` named counts, `nchr` number of chromosomes, `states` support.
# Output: matrix with rows = states and columns = chromosomes.
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

# Collapse a karyotype matrix to named whole-karyotype counts.
# Input: `M` cell-by-chromosome matrix.
# Output: named `table` of unique karyotype strings.
matrix_to_named_counts <- function(M, sep=".") {
  nm <- apply(M, 1, paste, collapse=sep)
  as.numeric(table(factor(nm)))
  # but keep names
  tbl <- table(nm)
  tbl
}

# Generate an initial synthetic population from marginal histograms.
# Inputs: `H` histogram matrix, kernel width `h`, state support, and sample size.
# Output: simulated cell-by-chromosome matrix.
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

# Compute multinomial negative log-likelihood per chromosome.
# Inputs: observed histogram `H_obs`, predicted probabilities `P_sim`, and `eps`.
# Output: scalar negative log-likelihood.
nll_multinomial <- function(H_obs, P_sim, eps=1e-9){
  P <- pmax(P_sim, eps); P <- sweep(P, 2, colSums(P), "/")
  -sum(H_obs * log(P))
}

# Raise a square matrix to a non-negative integer power.
# Inputs: square matrix `A`, integer exponent `n`.
# Output: matrix `A^n`.
mat_pow <- function(A, n) {
  if (length(n) != 1L || n < 0 || n != as.integer(n)) stop("n must be a single non-negative integer")
  n <- as.integer(n)
  if (n == 0L) return(diag(nrow(A)))
  out <- diag(nrow(A))
  base <- A
  while (n > 0L) {
    if (n %% 2L == 1L) out <- out %*% base
    base <- base %*% base
    n <- n %/% 2L
  }
  out
}

# Build the one-division copy-number transition matrix.
# Inputs: missegregation rate `pmis`, WGD rate `pwgd`, and state support.
# Output: square transition matrix over copy-number states.
make_cn_transition <- function(pmis, pwgd = 0, states = 0:8) {
  states <- as.integer(states)
  if (any(states != seq(min(states), max(states)))) stop("states must be consecutive integers")
  if (pmis < 0) stop("pmis must be >= 0")
  if (pwgd < 0) stop("pwgd must be >= 0")
  if (max(states) * pmis + pwgd > 1) {
    stop("Need max(states) * pmis + pwgd <= 1 so that diagonal probabilities stay non-negative")
  }
  S <- length(states)
  Q <- matrix(0, S, S, dimnames = list(from = states, to = states))
  for (i in seq_along(states)) {
    x <- states[i]
    px <- x * pmis
    Q[i, i] <- 1 - px - pwgd
    if ((x - 1L) %in% states) Q[i, match(x - 1L, states)] <- Q[i, match(x - 1L, states)] + px / 2
    if ((x + 1L) %in% states) Q[i, match(x + 1L, states)] <- Q[i, match(x + 1L, states)] + px / 2
    if ((2L * x) %in% states) Q[i, match(2L * x, states)] <- Q[i, match(2L * x, states)] + pwgd / 2
  }
  Q
}

# Propagate one initial copy number through repeated divisions.
# Inputs: initial state `init_cn`, division count `n_div`, or precomputed `Q`.
# Output: list with raw probabilities, conditional probabilities, survival, and `Q`.
cn_distribution <- function(init_cn, n_div, pmis = NULL, pwgd = 0, Q = NULL, states = 0:8) {
  states <- as.integer(states)
  if (!(init_cn %in% states)) {
    stop("init_cn must be one of: ", paste(states, collapse = ", "))
  }
  if (is.null(Q)) Q <- make_cn_transition(pmis = pmis, pwgd = pwgd, states = states)
  v0 <- setNames(rep(0, length(states)), states)
  v0[as.character(init_cn)] <- 1
  v_raw <- drop(v0 %*% mat_pow(Q, n_div))
  names(v_raw) <- states
  surv <- sum(v_raw)
  v_cond <- if (surv > 0) v_raw / surv else rep(NA_real_, length(v_raw))
  names(v_cond) <- states
  list(raw = v_raw, cond = v_cond, survival = surv, Q = Q)
}

# Predict per-chromosome endpoint probabilities from a starting karyotype matrix.
# Inputs: `K0` start matrix, `n_div`, `pmis`, `pwgd`, and allowed `states`.
# Output: matrix of predicted marginal copy-number probabilities.
predict_marginal_cn_probs <- function(K0, n_div, pmis, pwgd = 0, states = 0:8) {
  K0 <- round(as.matrix(K0))
  if (!nrow(K0) || !ncol(K0)) {
    stop("predict_marginal_cn_probs() requires a non-empty K0 matrix.")
  }
  if (any(!(K0 %in% states))) {
    stop("predict_marginal_cn_probs() found copy-number states outside the allowed set: ",
         paste(setdiff(sort(unique(as.integer(K0))), states), collapse = ", "))
  }

  Q <- make_cn_transition(pmis = pmis, pwgd = pwgd, states = states)
  Qn <- mat_pow(Q, n_div)
  cond_by_state <- apply(Qn, 1, function(v_raw) {
    surv <- sum(v_raw)
    if (surv > 0) v_raw / surv else rep(NA_real_, length(v_raw))
  })
  cond_by_state <- t(cond_by_state)
  P <- matrix(0, nrow = length(states), ncol = ncol(K0),
              dimnames = list(states = states, chr = paste0("chr", seq_len(ncol(K0)))))

  for (j in seq_len(ncol(K0))) {
    init_counts <- table(factor(K0[, j], levels = states))
    init_weights <- as.numeric(init_counts) / sum(init_counts)
    P[, j] <- drop(t(init_weights) %*% cond_by_state)
  }

  sweep(P, 2, colSums(P), "/")
}

# Score one observed interval under the Markov approximation.
# Inputs: start matrix `K0`, end matrix `KT`, division count, and model parameters.
# Output: scalar negative log-likelihood.
get_interval_nll_markov <- function(p_mis, p_wgd = 0, K0, KT, n_div, states = 0:8, eps = 1e-9) {
  K0 <- round(as.matrix(K0))
  KT <- round(as.matrix(KT))
  if (any(!(KT %in% states))) {
    stop("get_interval_nll_markov() found observed copy-number states outside the allowed set: ",
         paste(setdiff(sort(unique(as.integer(KT))), states), collapse = ", "))
  }
  H_obs <- hist_mat_from_matrix(KT, states = states)
  P_pred <- predict_marginal_cn_probs(K0 = K0, n_div = n_div, pmis = p_mis, pwgd = p_wgd, states = states)
  nll_multinomial(H_obs, P_pred, eps = eps)
}

# Score a whole interval group by summing interval likelihoods.
# Inputs: `intervals` list and model parameters.
# Output: scalar group negative log-likelihood.
get_group_nll_markov <- function(p_mis, p_wgd = 0, intervals, states = 0:8, eps = 1e-9) {
  sum(vapply(intervals, function(interval) {
    get_interval_nll_markov(
      p_mis = p_mis,
      p_wgd = p_wgd,
      K0 = interval$start_karyotype,
      KT = interval$end_karyotype,
      n_div = interval$n_divisions,
      states = states,
      eps = eps
    )
  }, numeric(1)))
}

# Evaluate the Markov objective on a grid of `p_mis` and `p_wgd`.
# Inputs: grouped intervals and parameter grids.
# Output: data.frame of one row per group/parameter combination.
estimate_group_pmis_grid <- function(grouped_intervals,
                                     pmis_grid = 10^(seq(-8, -1, length.out = 19)),
                                     pwgd_grid = 0,
                                     states = 0:8,
                                     eps = 1e-9) {
  if (is.null(names(grouped_intervals))) {
    names(grouped_intervals) <- as.character(seq_along(grouped_intervals))
  }

  param_grid <- expand.grid(
    pmis = pmis_grid,
    pwgd = pwgd_grid,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  do.call(rbind, lapply(seq_along(grouped_intervals), function(i) {
    group_id <- names(grouped_intervals)[i]
    group_intervals <- grouped_intervals[[i]]
    data.frame(
      group_id = group_id,
      pmis = param_grid$pmis,
      pwgd = param_grid$pwgd,
      negll = vapply(seq_len(nrow(param_grid)), function(ii) {
        get_group_nll_markov(
          p_mis = param_grid$pmis[ii],
          p_wgd = param_grid$pwgd[ii],
          intervals = group_intervals,
          states = states,
          eps = eps
        )
      }, numeric(1)),
      n_members = length(group_intervals),
      stringsAsFactors = FALSE
    )
  }))
}

# Predict ploidy distributions implied by fitted marginal copy-number probabilities.
# Inputs: `group_intervals`, fitted `p_mis`, fitted `p_wgd`.
# Output: data.frame of ploidy support and probabilities by interval.
predict_group_ploidy_distribution <- function(group_intervals, p_mis, p_wgd = 0, states = 0:8) {
  do.call(rbind, lapply(seq_along(group_intervals), function(i) {
    interval <- group_intervals[[i]]
    P_pred <- predict_marginal_cn_probs(
      K0 = interval$start_karyotype,
      n_div = interval$n_divisions,
      pmis = p_mis,
      pwgd = p_wgd,
      states = states
    )
    sum_probs <- c(1)
    sum_support <- 0
    for (j in seq_len(ncol(P_pred))) {
      chr_probs <- P_pred[, j]
      next_probs <- rep(0, length(sum_probs) + length(chr_probs) - 1L)
      for (k in seq_along(chr_probs)) {
        next_probs[k:(k + length(sum_probs) - 1L)] <- next_probs[k:(k + length(sum_probs) - 1L)] +
          chr_probs[k] * sum_probs
      }
      sum_probs <- next_probs
    }
    ploidy_support <- seq(0, by = 1 / ncol(P_pred), length.out = length(sum_probs))
    data.frame(
      interval_id = i,
      ploidy = ploidy_support,
      prob = sum_probs,
      stringsAsFactors = FALSE
    )
  }))
}

# Extract observed ploidy values from group endpoint karyotypes.
# Input: `group_intervals` list.
# Output: data.frame with one observed ploidy value per endpoint cell.
observed_group_ploidy_distribution <- function(group_intervals) {
  do.call(rbind, lapply(seq_along(group_intervals), function(i) {
    interval <- group_intervals[[i]]
    KT <- round(as.matrix(interval$end_karyotype))
    data.frame(
      interval_id = i,
      ploidy = rowMeans(KT),
      stringsAsFactors = FALSE
    )
  }))
}

# Legacy stochastic likelihood wrapper around the C++ simulator.
# Inputs: start matrix `K0`, observed histograms `H_obs`, simulator controls, and `p_mis`.
# Output: scalar negative log-likelihood.
get_nll <- function(p_mis,K0,H_obs,n_steps,max_pop,cull_keep,seed, record_every=5){
  ensure_ksim_compiled()
  pop0 <- matrix_to_named_counts(gen_pop(round(K0)))
  res <- run_karyotype_neutral(
    initial_counts_named = pop0,
    rate = 1.0,
    p_misseg = p_mis,
    p_wgd = 0,
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

# Run the legacy stochastic simulator and return the final named population.
# Inputs: start matrix `K0`, simulator controls, and `p_mis`.
# Output: named counts table from the final simulation time point.
get_pop <- function(p_mis,K0,nsteps,max_pop,cull_keep,seed, record_every=5){
  ensure_ksim_compiled()
  pop0 <- matrix_to_named_counts(gen_pop(round(K0)))
  res <- run_karyotype_neutral(
    initial_counts_named = pop0,
    rate = 1.0,
    p_misseg = p_mis,
    p_wgd = 0,
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

# Bootstrap an initial named population from observed start cells.
# Inputs: start matrix `K0`, `bottleneck_size`, optional RNG seed.
# Output: named whole-karyotype count table.
sample_initial_population <- function(K0, bottleneck_size = 2000, seed = NULL) {
  K0 <- round(as.matrix(K0))
  if (!is.null(seed)) set.seed(seed)
  idx <- sample(seq_len(nrow(K0)), size = bottleneck_size, replace = TRUE)
  matrix_to_named_counts(K0[idx, , drop = FALSE])
}

# Expand named whole-karyotype counts back into a cell-by-chromosome matrix.
# Inputs: `tbl` named count vector, `nchr`, and separator regex.
# Output: numeric matrix with one row per simulated cell.
named_counts_to_matrix <- function(tbl, nchr = 22, sep = "\\.") {
  if (!length(tbl)) {
    return(matrix(numeric(), nrow = 0, ncol = nchr))
  }
  out <- do.call(rbind, lapply(seq_along(tbl), function(i) {
    matrix(as.integer(strsplit(names(tbl)[i], sep)[[1]]),
           nrow = as.integer(tbl[i]),
           ncol = nchr,
           byrow = TRUE)
  }))
  mode(out) <- "numeric"
  out
}

# Run one forward simulation for a single interval under fitted parameters.
# Inputs: one `interval`, fitted `p_mis`/`p_wgd`, and simulation controls.
# Output: list containing trajectory, final counts, and final cell matrix.
simulate_interval_forward <- function(interval,
                                      p_mis,
                                      p_wgd,
                                      bottleneck_size = 2000,
                                      expansion_factor = 32,
                                      cull_keep = NULL,
                                      divide_all_each_step = TRUE,
                                      seed = 1,
                                      record_every = 25) {
  ensure_ksim_compiled()
  if (is.null(cull_keep)) {
    cull_keep <- 1 / expansion_factor
  }
  max_pop <- as.integer(round(bottleneck_size / cull_keep))
  initial_counts <- sample_initial_population(interval$start_karyotype, bottleneck_size = bottleneck_size, seed = seed)
  res <- run_karyotype_neutral(
    initial_counts_named = initial_counts,
    rate = 1.0,
    p_misseg = p_mis,
    p_wgd = p_wgd,
    dt = 1.0,
    n_steps = as.integer(interval$n_divisions),
    max_pop = max_pop,
    cull_keep = cull_keep,
    record_every = record_every,
    divide_all_each_step = divide_all_each_step,
    seed = seed
  )
  final_counts <- tail(res, 1)[[1]]
  final_matrix <- named_counts_to_matrix(final_counts, nchr = ncol(interval$start_karyotype))
  list(
    interval = interval,
    p_mis = p_mis,
    p_wgd = p_wgd,
    bottleneck_size = bottleneck_size,
    max_pop = max_pop,
    cull_keep = cull_keep,
    divide_all_each_step = divide_all_each_step,
    trajectory = res,
    final_counts = final_counts,
    final_matrix = final_matrix
  )
}

# Repeat forward simulation across all intervals in one fitted group.
# Inputs: `group_intervals`, one-row `fit_row`, and replication controls.
# Output: nested list indexed by interval then replicate.
simulate_group_forward <- function(group_intervals,
                                   fit_row,
                                   bottleneck_size = 2000,
                                   expansion_factor = 32,
                                   n_reps = 3,
                                   divide_all_each_step = TRUE,
                                   seed = 1,
                                   record_every = 25) {
  lapply(seq_along(group_intervals), function(i) {
    lapply(seq_len(n_reps), function(rep_idx) {
      simulate_interval_forward(
        interval = group_intervals[[i]],
        p_mis = fit_row$pmis,
        p_wgd = fit_row$pwgd,
        bottleneck_size = bottleneck_size,
        expansion_factor = expansion_factor,
        divide_all_each_step = divide_all_each_step,
        seed = seed + rep_idx + 1000L * i,
        record_every = record_every
      )
    })
  })
}

# Compare observed and simulated populations with a Wasserstein permutation test.
# Inputs: observed and simulated karyotype matrices plus permutation settings.
# Output: list with Wasserstein statistic, permutation p-value, and null stats.
wasserstein_permutation_test <- function(observed_matrix,
                                         simulated_matrix,
                                         n_perm = 200,
                                         seed = 1) {
  if (!requireNamespace("transport", quietly = TRUE)) {
    stop("wasserstein_permutation_test() requires the 'transport' package.")
  }

  observed_matrix <- round(as.matrix(observed_matrix))
  simulated_matrix <- round(as.matrix(simulated_matrix))
  if (!nrow(observed_matrix) || !nrow(simulated_matrix)) {
    return(list(statistic = NA_real_, p_value = NA_real_, perm_stats = numeric()))
  }
  if (ncol(observed_matrix) != ncol(simulated_matrix)) {
    stop("Observed and simulated matrices must have the same number of columns.")
  }

  pooled <- rbind(observed_matrix, simulated_matrix)
  n_obs <- nrow(observed_matrix)
  n_sim <- nrow(simulated_matrix)

  make_wpp <- function(M) {
    samstr <- table(apply(M, 1, paste, collapse = "."))
    coords <- matrix(as.numeric(unlist(strsplit(names(samstr), split = "[.]"))),
                     ncol = ncol(M),
                     byrow = TRUE)
    mass <- as.numeric(samstr)
    mass <- mass / sum(mass)
    transport::wpp(coords, mass = mass)
  }

  obs_wpp <- make_wpp(observed_matrix)
  sim_wpp <- make_wpp(simulated_matrix)
  stat_obs <- transport::wasserstein(obs_wpp, sim_wpp)

  set.seed(seed)
  perm_stats <- replicate(n_perm, {
    idx_obs <- sample(seq_len(nrow(pooled)), size = n_obs, replace = FALSE)
    idx_sim <- setdiff(seq_len(nrow(pooled)), idx_obs)
    perm_obs <- pooled[idx_obs, , drop = FALSE]
    perm_sim <- pooled[idx_sim, , drop = FALSE]
    transport::wasserstein(make_wpp(perm_obs), make_wpp(perm_sim))
  })

  list(
    statistic = stat_obs,
    p_value = mean(perm_stats >= stat_obs),
    perm_stats = perm_stats
  )
}

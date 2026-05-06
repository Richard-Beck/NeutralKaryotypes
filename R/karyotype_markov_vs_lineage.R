# Pure-R implementation of:
# 1) single-chromosome viability-conditioned Markov chain with per-copy misseg
# 2) single-lineage karyotype simulator tracking the first daughter only
# 3) wrappers to compare the two approaches

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

kary_to_string <- function(x) paste(x, collapse = ".")

make_cn_transition <- function(pmis, states = 1:8) {
  states <- as.integer(states)
  if (any(states != seq(min(states), max(states)))) stop("states must be consecutive integers")
  if (pmis < 0) stop("pmis must be >= 0")
  if (max(states) * pmis > 1) stop("Need max(states) * pmis <= 1 so that diagonal probabilities stay non-negative")
  S <- length(states)
  Q <- matrix(0, S, S, dimnames = list(from = states, to = states))
  for (i in seq_along(states)) {
    x <- states[i]
    px <- x * pmis
    Q[i, i] <- 1 - px
    if ((x - 1L) %in% states) Q[i, match(x - 1L, states)] <- Q[i, match(x - 1L, states)] + px / 2
    if ((x + 1L) %in% states) Q[i, match(x + 1L, states)] <- Q[i, match(x + 1L, states)] + px / 2
  }
  Q
}

cn_distribution <- function(init_cn, n_steps, pmis = NULL, Q = NULL, states = 1:8) {
  states <- as.integer(states)
  if (!(init_cn %in% states)) stop("init_cn must be one of: ", paste(states, collapse = ", "))
  if (is.null(Q)) Q <- make_cn_transition(pmis = pmis, states = states)
  v0 <- setNames(rep(0, length(states)), states)
  v0[as.character(init_cn)] <- 1
  v_raw <- drop(v0 %*% mat_pow(Q, n_steps))
  names(v_raw) <- states
  surv <- sum(v_raw)
  v_cond <- if (surv > 0) v_raw / surv else rep(NA_real_, length(v_raw))
  names(v_cond) <- states
  list(raw = v_raw, cond = v_cond, survival = surv, Q = Q)
}

sample_cn_chain <- function(init_cn, n_steps, n, pmis = NULL, Q = NULL, states = 1:8) {
  d <- cn_distribution(init_cn = init_cn, n_steps = n_steps, pmis = pmis, Q = Q, states = states)
  if (!is.finite(d$survival) || d$survival <= 0) return(integer(0))
  sample(as.integer(states), n, replace = TRUE, prob = d$cond)
}

simulate_lineage <- function(kary0, n_div, pmis, viable = 1:8) {
  k <- as.integer(kary0)
  viable <- as.integer(viable)
  if (pmis < 0) stop("pmis must be >= 0")
  if (max(viable) * pmis > 1) stop("Need max(viable) * pmis <= 1 so that per-chromosome event probabilities stay <= 1")
  if (any(!(k %in% viable))) return(NULL)
  if (n_div == 0L || pmis == 0) return(k)

  div <- 0L
  while (div < n_div) {
    p_chr <- k * pmis
    q_event <- 1 - prod(1 - p_chr)
    if (q_event <= 0) return(k)

    wait <- rgeom(1L, q_event)
    if (div + wait >= n_div) return(k)
    div <- div + wait + 1L

    hit <- rbinom(length(k), 1L, p_chr)
    while (!any(hit)) hit <- rbinom(length(k), 1L, p_chr)
    idx <- which(hit == 1L)
    delta <- sample(c(-1L, 1L), length(idx), replace = TRUE)
    k[idx] <- k[idx] + delta

    if (any(!(k %in% viable))) return(NULL)
  }
  k
}

simulate_population_lineage <- function(kary0, n_div, pmis, M, viable = 1:8) {
  sims <- replicate(M, simulate_lineage(kary0, n_div, pmis, viable), simplify = FALSE)
  keep <- Filter(Negate(is.null), sims)
  if (!length(keep)) return(list(samples = list(), counts = structure(integer(0), names = character(0)), survival = 0))
  labs <- vapply(keep, kary_to_string, FUN.VALUE = character(1))
  list(samples = keep, counts = sort(table(labs), decreasing = TRUE), survival = length(keep) / M)
}

simulate_population_markov <- function(kary0, n_div, pmis, M, states = 1:8) {
  Q <- make_cn_transition(pmis, states)
  samp <- do.call(cbind, lapply(kary0, function(x) sample_cn_chain(x, n_div, M, Q = Q, states = states)))
  labs <- apply(samp, 1L, kary_to_string)
  list(samples = samp, counts = sort(table(labs), decreasing = TRUE), Q = Q)
}

exact_joint_distribution <- function(kary0, n_div, pmis, states = 1:8) {
  dists <- lapply(kary0, function(x) cn_distribution(x, n_div, pmis = pmis, states = states)$cond)
  grids <- expand.grid(rep(list(as.integer(states)), length(kary0)), KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  P <- rep(1, nrow(grids))
  for (j in seq_along(kary0)) P <- P * dists[[j]][match(grids[[j]], as.integer(states))]
  nm <- apply(as.matrix(grids), 1L, kary_to_string)
  names(P) <- nm
  P[order(P, decreasing = TRUE)]
}

compare_populations <- function(counts1, counts2, B = 5000) {
  nm <- union(names(counts1), names(counts2))
  x <- setNames(integer(length(nm)), nm)
  y <- setNames(integer(length(nm)), nm)
  x[names(counts1)] <- as.integer(counts1)
  y[names(counts2)] <- as.integer(counts2)
  p1 <- x / sum(x)
  p2 <- y / sum(y)
  tv <- 0.5 * sum(abs(p1 - p2))
  tst <- if (length(nm) >= 2L && sum(x) > 0L && sum(y) > 0L) suppressWarnings(chisq.test(rbind(lineage = x, markov = y), simulate.p.value = TRUE, B = B)) else NULL
  list(total_variation = tv, chisq = tst, p1 = p1, p2 = p2)
}

set.seed(1)
pmis <- 0.001
kary0 <- c(2, 2)
N <- 250
M <- 50000

pop_lineage <- simulate_population_lineage(kary0, N, pmis, M)
M_surv <- sum(pop_lineage$counts)
pop_markov <- simulate_population_markov(kary0, N, pmis, M_surv)
cmp <- compare_populations(pop_lineage$counts, pop_markov$counts, B = 2000)

cat("Initial karyotype:", kary_to_string(kary0), "\n")
cat("pmis:", pmis, "\n")
cat("Tracked divisions:", N, "\n")
cat("Requested lineage replicates:", M, "\n")
cat("Lineage survival fraction:", pop_lineage$survival, "\n")
cat("Retained lineage samples:", M_surv, "\n\n")

cat("Single-chromosome conditional distribution for CN=2 after N divisions:\n")
print(round(cn_distribution(2, N, pmis = pmis)$cond, 6))
cat("\nTop lineage frequencies:\n")
print(round(head(sort(cmp$p1, decreasing = TRUE), 10), 6))
cat("\nTop markov frequencies:\n")
print(round(head(sort(cmp$p2, decreasing = TRUE), 10), 6))
cat("\nTotal variation distance:\n")
print(cmp$total_variation)
cat("\nChi-squared test comparing the two empirical distributions:\n")
print(cmp$chisq)
if (length(kary0) <= 6L) {
  cat("\nExact joint distribution from method 1 (top states):\n")
  print(round(head(exact_joint_distribution(kary0, N, pmis), 10), 6))
}

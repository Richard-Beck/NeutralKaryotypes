#!/usr/bin/env Rscript

# Usage: Rscript R/analyse_lineage.R <index> <n_cores>

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Error: Please provide the object index and number of cores. \nUsage: Rscript R/analyse_lineage.R <index> <n_cores>")
}

idx <- as.integer(args[1])
n_cores <- as.integer(args[2])

# Directory setup
outputDir <- "data/fits"
if (!dir.exists(outputDir)) {
  dir.create(outputDir)
}

outputPops <- "data/pops"
if (!dir.exists(outputPops)) {
  dir.create(outputPops)
}

# Load the main data
if (!file.exists("fit_objects.Rds")) {
  stop("fit_objects.Rds not found in the current directory.")
}
fit_objects <- readRDS("fit_objects.Rds")

# Select the specific object
target_obj <- fit_objects[[idx]]

# Extract parameters
K0 <- target_obj$K0
KT <- target_obj$KT
H0 <- target_obj$H0
HT <- target_obj$HT

# Setup Parallel Cluster
library(parallel)
cl <- makeCluster(n_cores)

# Source the estimation function
if (!file.exists("R/estimate_pmis.R")) {
  stop("R/estimate_pmis.R not found.")
}
source("R/estimate_pmis.R")
clusterEvalQ(cl, source("R/estimate_pmis.R"))

# --- Step 1: Estimate PMIS (Likelihood calculation) ---

pmis <- 10^(seq(-8, -1, length.out=19))
df <- expand.grid(pmis=pmis, rep=1:3)

# Time parameters
Ndoublings <- target_obj$delta_pass * 5
dt_growth <- log(2)
total_time <- dt_growth * Ndoublings
dt_step <- 0.1
n_steps <- round(total_time / dt_step)

# FIX: Added "n_steps" to the export list
clusterExport(cl, c("K0", "HT", "df", "n_steps"))

# Calculate Negative Log Likelihood
res <- parSapplyLB(cl, 1:nrow(df), function(i) {
  get_nll(
    p_mis = df$pmis[i],
    K0 = K0,
    H_obs = HT,
    n_steps = n_steps, 
    max_pop = 0.3*10^6,
    cull_keep = 1/(2^5),
    seed = i
  )
})

resdf <- df
resdf$negll <- res

# Update object with likelihood results
target_obj$resdf <- resdf

# --- Step 2: Generate Populations ---

# Find the p_mis corresponding to the maximum likelihood (min negll)
best_fit_index <- which.min(resdf$negll)
best_pmis <- resdf$pmis[best_fit_index]

Nreps <- 150

# FIX: Export "best_pmis" so workers can see it
clusterExport(cl, c("best_pmis"))

# Generate populations using the best estimated p_mis
pops <- parLapplyLB(cl, 1:Nreps, function(i) {
  get_pop(
    p_mis = best_pmis,
    K0 = K0,
    nsteps = n_steps,
    max_pop = 0.3*10^6,
    cull_keep = 1/(2^5),
    seed = i
  )
})

# Update object with generated populations
target_obj$pops <- pops

out_file <- file.path(outputPops, paste0(target_obj$id_end, ".Rds"))
saveRDS(pops, out_file)

stopCluster(cl)
out_file <- file.path(outputDir, paste0(target_obj$id_end, ".Rds"))

#saveRDS(target_obj, out_file)
# --- Step 3: Null Hypothesis Testing (Transport Distance) ---

test_null <- function(pop, KT) {
  library(transport)
  
  sampleToInt <- function(pp, Nsam=20) {
    samstr <- sample(names(pp), size = Nsam, prob = pp)
    samstr <- table(samstr)
    # Ensure consistent matrix dimensions for 22 chromosomes/bins
    coords <- matrix(as.numeric(unlist(strsplit(names(samstr), split="[.]"))), ncol=22, byrow = TRUE)
    list(coords = coords, counts = samstr)
  }
  
  Nsam <- nrow(KT)
  
  # Process Target (KT)
  samstr <- table(apply(KT, 1, paste, collapse="."))
  smm <- list(
    coords = matrix(as.numeric(unlist(strsplit(names(samstr), split="[.]"))), ncol=22, byrow = TRUE),
    counts = samstr
  )
  wtst <- wpp(smm$coords, mass = smm$counts)
  
  # 1. Internal distances (Null distribution)
  dnull <- sapply(1:length(pop), function(i) {
    sapply(1:length(pop), function(j) {
      if (i >= j) return(NaN)
      
      popi <- sampleToInt(pop[[i]], Nsam)
      popj <- sampleToInt(pop[[j]], Nsam)
      
      wi <- wpp(popi$coords, mass = popi$counts)
      wj <- wpp(popj$coords, mass = popj$counts)
      
      wasserstein(wi, wj)
    })
  })
  dnull <- dnull[is.finite(dnull)]
  
  # 2. Test distances (Populations vs Target)
  dtst <- sapply(1:length(pop), function(j) {
    popj <- sampleToInt(pop[[j]], Nsam)
    
    wj <- wpp(popj$coords, mass = popj$counts)
    wasserstein(wj, wtst)
  })
  
  tst_stat <- median(dtst)
  empp <- mean(dnull > tst_stat)
  
  list(dnull = dnull, dtst = dtst, empp = empp)
}

# Run test (run locally as it's not computationally heavy enough to require cluster overhead)
res_test <- test_null(target_obj$pops, round(KT))

# Update object with test results
target_obj$test_null_res <- res_test

# --- Step 4: Save Result ---

# Construct filename
out_file <- file.path(outputDir, paste0(target_obj$id_end, ".Rds"))
target_obj$pops <- NULL ## otherwise outputs are huge
saveRDS(target_obj, out_file)

cat(sprintf("Finished processing object %d. Saved to %s\n", idx, out_file))
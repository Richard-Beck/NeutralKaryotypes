#!/usr/bin/env Rscript

# --- Configuration ---
n_cores_local <- parallel::detectCores() - 2  
outputDir <- "data/fits"
outputPops <- "data/pops"

if (!dir.exists(outputDir)) dir.create(outputDir, recursive = TRUE)
if (!dir.exists(outputPops)) dir.create(outputPops, recursive = TRUE)

# --- Load Data ---
if (!file.exists("fit_objects.Rds")) stop("fit_objects.Rds not found.")
fit_objects <- readRDS("fit_objects.Rds")

library(parallel)
library(transport)
library(pbapply) # Load the progress bar package

# Configure pbapply to print to stdout (so it shows up in RStudio/Terminal)
pboptions(type = "timer", style = 3, char = "=")

# --- Main Processing Loop ---
for (idx in seq_along(fit_objects)) {
  
  target_obj <- fit_objects[[idx]]
  if (is.null(target_obj)) next 
  
  out_file <- file.path(outputDir, paste0(target_obj$id_end, ".Rds"))
  
  if (file.exists(out_file)) {
    message(sprintf("[%d/%d] Skipping %s (Output exists)", idx, length(fit_objects), target_obj$id_end))
    next
  }
  
  message(sprintf("[%d/%d] Processing %s...", idx, length(fit_objects), target_obj$id_end))
  
  # Extract parameters
  K0 <- target_obj$K0
  KT <- target_obj$KT
  HT <- target_obj$HT
  
  # --- Setup Cluster ---
  cl <- makeCluster(n_cores_local)
  
  # Initialize workers with C++ compilation
  clusterEvalQ(cl, {
    source("R/estimate_pmis.R")
  })
  
  # --- Step 1: Estimate PMIS (Optimized for Speed) ---
  # Search 19 points, 1 replicate each (sufficient for large N)
  pmis <- 10^(seq(-8, -1, length.out=19))
  df <- expand.grid(pmis=pmis, rep=1) 
  
  Ndoublings <- target_obj$delta_pass * 5
  dt_growth <- log(2)
  total_time <- dt_growth * Ndoublings
  dt_step <- 0.1
  n_steps <- round(total_time / dt_step)
  
  clusterExport(cl, c("K0", "HT", "df", "n_steps"))
  
  message("  Step 1: Grid Search Likelihood...")
  # Replaced parSapplyLB with pbsapply
  res <- pbsapply(1:nrow(df), function(i) {
    get_nll(
      p_mis = df$pmis[i],
      K0 = K0,
      H_obs = HT,
      n_steps = n_steps, 
      max_pop = 1e5, # 100k Population for speed/accuracy balance
      cull_keep = 1/(2^5),
      seed = i
    )
  }, cl = cl) # Pass cluster here
  
  resdf <- df
  resdf$negll <- res
  target_obj$resdf <- resdf
  
  # --- Step 2: Generate Populations ---
  best_fit_index <- which.min(resdf$negll)
  best_pmis <- resdf$pmis[best_fit_index]
  Nreps <- 100 # 100 Replicates for p-value resolution ~0.01
  
  clusterExport(cl, c("best_pmis"))
  
  message("  Step 2: Generating Null Populations...")
  # Replaced parLapplyLB with pblapply
  pops <- pblapply(1:Nreps, function(i) {
    get_pop(
      p_mis = best_pmis,
      K0 = K0,
      nsteps = n_steps,
      max_pop = 1e5,
      cull_keep = 1/(2^5),
      seed = i
    )
  }, cl = cl) # Pass cluster here
  
  stopCluster(cl)
  
  target_obj$pops <- pops
  saveRDS(pops, file.path(outputPops, paste0(target_obj$id_end, ".Rds")))

  # --- Step 3: Null Hypothesis Testing (Hybrid Approach) ---
  message("  Step 3: Calculating Wasserstein Distances...")
  
  test_null <- function(pop, KT) {
    sampleToInt <- function(pp, Nsam=20) {
      if (length(pp) == 0) return(NULL)
      samstr <- sample(names(pp), size = Nsam, prob = pp, replace = TRUE)
      samstr <- table(samstr)
      coords <- matrix(as.numeric(unlist(strsplit(names(samstr), split="[.]"))), ncol=22, byrow = TRUE)
      list(coords = coords, counts = samstr)
    }
    
    Nsam <- nrow(KT) # Dynamically matches the number of cells in the specific flask
    samstr <- table(apply(KT, 1, paste, collapse="."))
    smm <- list(
      coords = matrix(as.numeric(unlist(strsplit(names(samstr), split="[.]"))), ncol=22, byrow = TRUE),
      counts = samstr
    )
    wtst <- wpp(smm$coords, mass = smm$counts)
    
    # 1. Internal distances
    dnull <- sapply(1:length(pop), function(i) {
      sapply(1:length(pop), function(j) {
        if (i >= j) return(NaN)
        popi <- sampleToInt(pop[[i]], Nsam)
        popj <- sampleToInt(pop[[j]], Nsam)
        if (is.null(popi) || is.null(popj)) return(NaN)
        wi <- wpp(popi$coords, mass = popi$counts)
        wj <- wpp(popj$coords, mass = popj$counts)
        wasserstein(wi, wj)
      })
    })
    dnull <- dnull[is.finite(dnull)]
    
    # 2. Test distances
    dtst <- sapply(1:length(pop), function(j) {
      popj <- sampleToInt(pop[[j]], Nsam)
      if (is.null(popj)) return(NaN)
      wj <- wpp(popj$coords, mass = popj$counts)
      wasserstein(wj, wtst)
    })
    
    tst_stat <- median(dtst, na.rm=TRUE)
    empp <- mean(dnull > tst_stat, na.rm=TRUE)
    
    list(dnull = dnull, dtst = dtst, empp = empp)
  }
  
  # A. Test the merged population (Optional, but good for reference)
  res_merged <- test_null(target_obj$pops, round(target_obj$KT))
  
  # B. Test EACH replicate independently against the shared Joint Null Model
  res_indiv <- lapply(target_obj$KT_list, function(kt_indiv) {
    test_null(target_obj$pops, round(kt_indiv))
  })
  
  target_obj$test_null_res <- list(
    merged = res_merged,
    individuals = res_indiv
  )
  
  # --- Step 4: Save & Clean ---
  target_obj$pops <- NULL 
  saveRDS(target_obj, out_file)
  
  message(sprintf("Finished %s. Saved to %s\n", target_obj$id_end, out_file))
}

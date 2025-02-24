#' Simulate Data and Compare LDWAS and GWASbyCor
#'
#' This function simulates data, runs LDWAS and GWASbyCor, and compares their performance.
#'
#' @param n Number of individuals.
#' @param m Number of SNPs.
#' @param num_causal Number of causal SNPs.
#' @param effect_size Effect size of causal SNPs.
#' @param num_covariates Number of covariates.
#' @param num_replicates Number of replicates.
#' @return A data frame containing the results for each replicate.
#' @export
simulate_and_compare <- function(n = 100, m = 1000, num_causal = 10, effect_size = 0.5, num_covariates = 2, num_replicates = 30) {
  results <- data.frame(
    replicate = integer(),
    method = character(),
    power = numeric(),
    fpr = numeric(),
    time = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (i in 1:num_replicates) {
    # Simulate data
    X <- matrix(sample(0:2, n * m, replace = TRUE), nrow = n)  # Genotype data
    C <- matrix(rnorm(n * num_covariates), nrow = n)  # Covariate data
    causal_snps <- sample(1:m, num_causal)  # Randomly select causal SNPs
    y <- rowSums(X[, causal_snps] * effect_size) + rnorm(n)  # Phenotype data
    
    # SNP information
    snp_info <- data.frame(SNP = paste0("rs", 1:m), chr = rep(1:10, each = m / 10), pos = 1:m)
    
    # Run LDWAS
    start_time <- Sys.time()
    ldwas_results <- run_LDWAS(y, X, C, snp_info)
    ldwas_time <- as.numeric(Sys.time() - start_time, units = "secs")
    
    # Run GWASbyCor
    start_time <- Sys.time()
    gwasbycor_results <- GWASbyCor(y, X, C)
    gwasbycor_time <- as.numeric(Sys.time() - start_time, units = "secs")
    
    # Calculate power and FPR for LDWAS
    ldwas_significant <- ldwas_results$significant_snps$SNP
    ldwas_power <- sum(ldwas_significant %in% paste0("rs", causal_snps)) / num_causal
    ldwas_fpr <- sum(!ldwas_significant %in% paste0("rs", causal_snps)) / (m - num_causal)
    
    # Calculate power and FPR for GWASbyCor
    gwasbycor_significant <- gwasbycor_results$significant_snps$SNP
    gwasbycor_power <- sum(gwasbycor_significant %in% paste0("rs", causal_snps)) / num_causal
    gwasbycor_fpr <- sum(!gwasbycor_significant %in% paste0("rs", causal_snps)) / (m - num_causal)
    
    # Store results
    results <- rbind(results, data.frame(
      replicate = i,
      method = "LDWAS",
      power = ldwas_power,
      fpr = ldwas_fpr,
      time = ldwas_time
    ))
    results <- rbind(results, data.frame(
      replicate = i,
      method = "GWASbyCor",
      power = gwasbycor_power,
      fpr = gwasbycor_fpr,
      time = gwasbycor_time
    ))
  }
  
  return(results)
}

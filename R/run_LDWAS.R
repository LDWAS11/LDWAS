run_LDWAS <- function(y, X, C, snp_info, num_pcs = 10, threshold = 1e-6, alpha = 0.05) {
  # Perform PCA and exclude dependent PCs
  pcs <- PCA_Cofactors(X, C, threshold = threshold)

  # Use the first `num_pcs` PCs (or fewer if some were excluded)
  selected_pcs <- pcs[, 1:min(num_pcs, ncol(pcs)), drop = FALSE]

  # Combine covariates and selected PCs
  covariates <- cbind(C, selected_pcs)

  # Perform GWAS
  p_values <- GWAS_GLM(y, X, covariates)

  # Calculate genome-wide threshold
  gwas_threshold <- genome_wide_threshold(p_values, alpha = alpha)

  # Identify significant SNPs
  significant_snps <- get_significant_snps(p_values, gwas_threshold, snp_info)

  # Calculate MAF
  maf <- calculate_maf(X)

  # Analyze MAF of significant SNPs
  maf_analysis <- analyze_maf(significant_snps, maf)

  # Save Manhattan plot to a file
  png("manhattan_plot.png", width = 1200, height = 800)
  manhattan_plot(p_values, snp_info, threshold = gwas_threshold)
  dev.off()

  # Save QQ plot to a file
  png("qq_plot.png", width = 800, height = 800)
  qq_plot(p_values)
  dev.off()

  # Return all results
  results <- list(
    p_values = p_values,
    threshold = gwas_threshold,
    significant_snps = significant_snps,
    maf_analysis = maf_analysis,
    manhattan_plot = "manhattan_plot.png",
    qq_plot = "qq_plot.png"
  )

  return(results)
}

<<<<<<< HEAD
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

=======
#' Run LDWAS Analysis
#'
#' This function performs GWAS using GLM while incorporating principal components (PCs) as covariates and excluding PCs that are linearly dependent on user-provided covariates. It also generates a Manhattan plot, QQ plot, genome-wide threshold, list of associated SNPs, and MAF analysis.
#'
#' @param y A numeric vector of phenotypes (n x 1).
#' @param X A numeric matrix of genotypes (n x m).
#' @param C A numeric matrix of covariates (n x t).
#' @param snp_info A data frame containing SNP information (e.g., SNP ID, chromosome, position).
#' @param num_pcs The number of principal components to extract (default: 10).
#' @param threshold A numeric value specifying the tolerance for linear dependence (default: 1e-6).
#' @param alpha A numeric value for the significance level (default: 0.05).
#' @return A list containing:
#'   - p_values: A numeric vector of p-values (1 x m) for each marker.
#'   - threshold: The genome-wide significance threshold.
#'   - significant_snps: A data frame of significantly associated SNPs.
#'   - maf_analysis: A summary of MAF analysis for significant SNPs.
#'   - manhattan_plot: A Manhattan plot.
#'   - qq_plot: A QQ plot.
#' @examples
#' y <- rnorm(100)  # Simulate phenotype data
#' X <- matrix(sample(0:2, 100 * 1000, replace = TRUE), nrow = 100)  # Simulate genotype data
#' C <- matrix(rnorm(100 * 2), nrow = 100)  # Simulate covariate data
#' snp_info <- data.frame(SNP = paste0("rs", 1:1000), chr = rep(1:10, each = 100), pos = 1:1000)  # Simulate SNP info
#' results <- run_LDWAS(y, X, C, snp_info, num_pcs = 10)  # Run the function
#' @export
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
  
  # Generate Manhattan plot
  manhattan_plot <- manhattan_plot(p_values, snp_info, threshold = gwas_threshold)
  
  # Generate QQ plot
  qq_plot <- qq_plot(p_values)
  
>>>>>>> d72c2305b517db46d74d7c36ff9f1ebbcec181cc
  # Return all results
  results <- list(
    p_values = p_values,
    threshold = gwas_threshold,
    significant_snps = significant_snps,
    maf_analysis = maf_analysis,
<<<<<<< HEAD
    manhattan_plot = "manhattan_plot.png",
    qq_plot = "qq_plot.png"
  )

=======
    manhattan_plot = manhattan_plot,
    qq_plot = qq_plot
  )
  
>>>>>>> d72c2305b517db46d74d7c36ff9f1ebbcec181cc
  return(results)
}

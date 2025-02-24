#' Calculate Minor Allele Frequency (MAF)
#'
#' This function calculates the Minor Allele Frequency (MAF) for each SNP.
#'
#' @param X A numeric matrix of genotypes (n x m), where genotypes are coded as 0, 1, or 2.
#' @return A numeric vector of MAF values (1 x m).
#' @examples
#' X <- matrix(sample(0:2, 100 * 1000, replace = TRUE), nrow = 100)  # Simulate genotype data
#' maf <- calculate_maf(X)
#' print(maf)
calculate_maf <- function(X) {
  # Check if X is a numeric matrix
  if (!is.matrix(X) || !is.numeric(X)) {
    stop("X must be a numeric matrix.")
  }
  
  # Check if X contains valid genotype values (0, 1, or 2)
  if (any(!X %in% c(0, 1, 2))) {
    stop("X must contain genotype values coded as 0, 1, or 2.")
  }
  
  # Calculate MAF
  maf <- colMeans(X) / 2
  maf <- pmin(maf, 1 - maf)  # MAF is the smaller of the two allele frequencies
  return(maf)
}

#' Analyze MAF of Significant SNPs
#'
#' This function analyzes the Minor Allele Frequency (MAF) of significantly associated SNPs.
#'
#' @param significant_snps A data frame of significantly associated SNPs. Must contain a column named "SNP".
#' @param maf A numeric vector of MAF values (1 x m).
#' @return A summary of MAF analysis, including mean, median, min, and max MAF.
#' @examples
#' X <- matrix(sample(0:2, 100 * 1000, replace = TRUE), nrow = 100)  # Simulate genotype data
#' p_values <- runif(1000)  # Simulate p-values
#' snp_info <- data.frame(SNP = paste0("rs", 1:1000), chr = rep(1:10, each = 100), pos = 1:1000)
#' threshold <- genome_wide_threshold(p_values)
#' significant_snps <- get_significant_snps(p_values, threshold, snp_info)
#' maf <- calculate_maf(X)
#' maf_analysis <- analyze_maf(significant_snps, maf)
#' print(maf_analysis)
analyze_maf <- function(significant_snps, maf) {
  # Check if significant_snps is a data frame and contains a "SNP" column
  if (!is.data.frame(significant_snps) || !"SNP" %in% colnames(significant_snps)) {
    stop("significant_snps must be a data frame with a column named 'SNP'.")
  }
  
  # Check if maf is a numeric vector
  if (!is.numeric(maf)) {
    stop("maf must be a numeric vector.")
  }
  
  # Handle case where there are no significant SNPs
  if (nrow(significant_snps) == 0) {
    return(list(
      mean_maf = NA,
      median_maf = NA,
      min_maf = NA,
      max_maf = NA
    ))
  }
  
  # Add MAF to significant SNPs
  significant_snps$maf <- maf[significant_snps$SNP]
  
  # Summarize MAF analysis
  summary <- list(
    mean_maf = mean(significant_snps$maf),
    median_maf = median(significant_snps$maf),
    min_maf = min(significant_snps$maf),
    max_maf = max(significant_snps$maf)
  )
  
  return(summary)
}

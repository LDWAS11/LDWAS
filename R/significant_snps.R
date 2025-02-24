#' Identify Significant SNPs
#'
#' This function identifies SNPs that are significantly associated with the phenotype.
#'
#' @param p_values A numeric vector of p-values (1 x m).
#' @param threshold A numeric value for the significance threshold.
#' @param snp_info A data frame containing SNP information (e.g., SNP ID, chromosome, position).
#' @return A data frame of significantly associated SNPs.
#' @examples
#' p_values <- runif(1000)  # Simulate p-values
#' snp_info <- data.frame(SNP = paste0("rs", 1:1000), chr = rep(1:10, each = 100), pos = 1:1000)
#' threshold <- genome_wide_threshold(p_values)
#' significant_snps <- get_significant_snps(p_values, threshold, snp_info)
#' print(significant_snps)
#' @export
get_significant_snps <- function(p_values, threshold, snp_info) {
  # Identify significant SNPs
  significant_indices <- which(p_values < threshold)
  significant_snps <- snp_info[significant_indices, ]
  significant_snps$p_value <- p_values[significant_indices]
  return(significant_snps)
}

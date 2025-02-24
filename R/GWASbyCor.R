#' Perform GWAS Using Correlation
#'
#' This function performs GWAS by testing the correlation between each SNP and the phenotype.
#'
#' @param y A numeric vector of phenotypes (n x 1).
#' @param X A numeric matrix of genotypes (n x m).
#' @param C A numeric matrix of covariates (n x t).
#' @return A list containing:
#'   - p_values: A numeric vector of p-values (1 x m) for each marker.
#'   - significant_snps: A data frame of significantly associated SNPs.
#' @examples
#' y <- rnorm(100)  # Simulate phenotype data
#' X <- matrix(sample(0:2, 100 * 1000, replace = TRUE), nrow = 100)  # Simulate genotype data
#' C <- matrix(rnorm(100 * 2), nrow = 100)  # Simulate covariate data
#' results <- GWASbyCor(y, X, C)
#' @export
GWASbyCor <- function(y, X, C) {
  p_values <- numeric(ncol(X))  # Initialize a vector to store p-values
  
  for (i in 1:ncol(X)) {
    # Fit a linear model to test the correlation between SNP and phenotype
    model <- lm(y ~ X[, i] + C)
    p_values[i] <- summary(model)$coefficients[2, 4]  # Extract p-value for the SNP
  }
  
  # Identify significant SNPs (e.g., using a Bonferroni threshold)
  threshold <- 0.05 / ncol(X)
  significant_snps <- data.frame(
    SNP = paste0("rs", 1:ncol(X))[p_values < threshold],
    p_value = p_values[p_values < threshold]
  )
  
  return(list(
    p_values = p_values,
    significant_snps = significant_snps
  ))
}

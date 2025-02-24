#' Perform GWAS using GLM
#'
#' This function performs genome-wide association studies (GWAS) using a generalized linear model (GLM).
#'
#' @param y A numeric vector of phenotypes (n x 1).
#' @param X A numeric matrix of genotypes (n x m).
#' @param C A numeric matrix of covariates (n x t).
#' @param family A family function for the GLM (default is gaussian()).
#' @param num_cores Number of cores for parallel processing (default is 1).
#' @return A numeric vector of p-values (1 x m) for each marker.
#' @examples
#' y <- rnorm(100)  # Simulate phenotype data
#' X <- matrix(sample(0:2, 100 * 1000, replace = TRUE), nrow = 100)  # Simulate genotype data
#' C <- matrix(rnorm(100 * 2), nrow = 100)  # Simulate covariate data
#' p_values <- GWAS_GLM(y, X, C)  # Run the function
#' @export
GWAS_GLM <- function(y, X, C, family = gaussian(), num_cores = 1) {
  # Input validation
  if (!is.numeric(y) || length(y) != nrow(X)) stop("Invalid phenotype vector 'y'.")
  if (!is.matrix(X) || !is.numeric(X)) stop("Invalid genotype matrix 'X'.")
  if (!is.matrix(C) || !is.numeric(C)) stop("Invalid covariate matrix 'C'.")
  if (nrow(X) != nrow(C)) stop("Number of rows in 'X' and 'C' must match.")
  if (any(is.na(y)) || any(is.nan(y)) || any(is.infinite(y))) stop("Phenotype vector 'y' contains NA, NaN, or Inf values.")
  if (any(is.na(X)) || any(is.nan(X)) || any(is.infinite(X))) stop("Genotype matrix 'X' contains NA, NaN, or Inf values.")
  if (any(is.na(C)) || any(is.nan(C)) || any(is.infinite(C))) stop("Covariate matrix 'C' contains NA, NaN, or Inf values.")
  
  # Convert C to a data frame for proper formula handling
  C_df <- as.data.frame(C)
  colnames(C_df) <- paste0("C", 1:ncol(C))
  
  # Parallel processing (if num_cores > 1)
  if (num_cores > 1) {
    library(parallel)
    cl <- makeCluster(num_cores)
    clusterExport(cl, c("y", "C_df", "family", "glm"))
    
    p_values <- parLapply(cl, 1:ncol(X), function(i) {
      model <- glm(y ~ X[, i] + ., data = C_df, family = family)
      return(summary(model)$coefficients[2, 4])
    })
    
    stopCluster(cl)
    p_values <- unlist(p_values)
  } else {
    # Sequential processing
    p_values <- numeric(ncol(X))
    for (i in 1:ncol(X)) {
      model <- glm(y ~ X[, i] + ., data = C_df, family = family)
      p_values[i] <- summary(model)$coefficients[2, 4]
    }
  }
  
  return(p_values)
}

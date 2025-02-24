#' Perform PCA and exclude dependent PCs
#'
#' This function performs principal component analysis (PCA) on genotype data and excludes PCs that are linearly dependent on covariates.
#'
#' @param X A numeric matrix of genotypes (n x m).
#' @param C A numeric matrix of covariates (n x t).
#' @param threshold A numeric value specifying the tolerance for linear dependence (default: 1e-6).
#' @return A numeric matrix of selected PCs (n x k), where k is the number of non-dependent PCs.
#' @examples
#' X <- matrix(sample(0:2, 100 * 1000, replace = TRUE), nrow = 100)  # Simulate genotype data
#' C <- matrix(rnorm(100 * 2), nrow = 100)  # Simulate covariate data
#' pcs <- PCA_Cofactors(X, C)  # Run the function
#' @export
PCA_Cofactors <- function(X, C, threshold = 1e-6) {
  
  # Check dimensions of inputs
  cat("Dimensions of X:", dim(X), "\n")
  cat("Dimensions of C:", dim(C), "\n")
  
  # Ensure that the number of rows in X and C match
  if (nrow(X) != nrow(C)) stop("The number of rows in X and C must match.")
  
  # Remove constant or zero columns from X
  constant_cols <- apply(X, 2, function(col) length(unique(col)) == 1)
  zero_cols <- apply(X, 2, function(col) all(col == 0))
  problematic_cols <- constant_cols | zero_cols
  X <- X[, !problematic_cols, drop = FALSE]  # Ensure X remains a matrix
  
  # Debugging: Print dimensions after removing problematic columns
  cat("Dimensions of X after removing constant/zero columns:", dim(X), "\n")
  
  # Check for NA, NaN, or infinite values in X
  if (any(is.na(X))) stop("Matrix 'X' contains NA values.")
  if (any(is.nan(X))) stop("Matrix 'X' contains NaN values.")
  if (any(is.infinite(X))) stop("Matrix 'X' contains infinite values.")
  
  # Perform PCA using SVD for better numerical stability
  X_scaled <- scale(X)  # Center and scale the genotype matrix
  pca <- svd(X_scaled)  # Singular Value Decomposition (SVD)
  pcs <- pca$u %*% diag(pca$d)  # Principal components
  
  # Debugging: Print dimensions of PCs
  cat("Dimensions of PCs before exclusion:", dim(pcs), "\n")
  
  # Check for linear dependence between PCs and covariates
  dependent_pcs <- numeric(0)
  for (i in 1:ncol(pcs)) {
    model <- lm(pcs[, i] ~ C)
    if (sum(model$residuals^2) < threshold) {
      dependent_pcs <- c(dependent_pcs, i)
    }
  }
  
  # Debugging: Print dependent PCs
  cat("Dependent PCs:", dependent_pcs, "\n")
  
  # Exclude dependent PCs, but ensure at least one PC is retained
  if (length(dependent_pcs) > 0) {
    pcs <- pcs[, -dependent_pcs, drop = FALSE]
  }
  
  # If all PCs are excluded, retain the first PC
  if (ncol(pcs) == 0) {
    warning("All PCs were dependent on covariates. Retaining the first PC.")
    pcs <- pcs[, 1, drop = FALSE]
  }
  
  # Debugging: Print final dimensions of PCs
  cat("Final dimensions of PCs:", dim(pcs), "\n")
  
  return(pcs)
}

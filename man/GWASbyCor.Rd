% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GWASbyCor.R
\name{GWASbyCor}
\alias{GWASbyCor}
\title{Perform GWAS Using Correlation}
\usage{
GWASbyCor(y, X, C)
}
\arguments{
\item{y}{A numeric vector of phenotypes (n x 1).}

\item{X}{A numeric matrix of genotypes (n x m).}

\item{C}{A numeric matrix of covariates (n x t).}
}
\value{
A list containing:
\itemize{
\item p_values: A numeric vector of p-values (1 x m) for each marker.
\item significant_snps: A data frame of significantly associated SNPs.
}
}
\description{
This function performs GWAS by testing the correlation between each SNP and the phenotype.
}
\examples{
y <- rnorm(100)  # Simulate phenotype data
X <- matrix(sample(0:2, 100 * 1000, replace = TRUE), nrow = 100)  # Simulate genotype data
C <- matrix(rnorm(100 * 2), nrow = 100)  # Simulate covariate data
results <- GWASbyCor(y, X, C)
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/qq_plot.R
\name{qq_plot}
\alias{qq_plot}
\title{Generate a QQ Plot}
\usage{
qq_plot(p_values)
}
\arguments{
\item{p_values}{A numeric vector of p-values (1 x m).}
}
\value{
A QQ plot.
}
\description{
This function generates a QQ plot for GWAS results.
}
\examples{
p_values <- runif(1000)  # Simulate p-values
qq_plot(p_values)
}

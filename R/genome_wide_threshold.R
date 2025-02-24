#' Calculate Genome-Wide Significance Threshold
#'
#' This function calculates a genome-wide significance threshold using Bonferroni correction.
#'
#' @param p_values A numeric vector of p-values (1 x m).
#' @param alpha A numeric value for the significance level (default: 0.05).
#' @return A numeric value representing the genome-wide threshold.
#' @examples
#' p_values <- runif(1000)  # Simulate p-values
#' threshold <- genome_wide_threshold(p_values)
#' print(threshold)
#' @export
genome_wide_threshold <- function(p_values, alpha = 0.05) {
  # Bonferroni correction
  threshold <- alpha / length(p_values)
  return(threshold)
}

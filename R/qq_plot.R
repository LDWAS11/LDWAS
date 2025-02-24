#' Generate a QQ Plot
#'
#' This function generates a QQ plot for GWAS results.
#'
#' @param p_values A numeric vector of p-values (1 x m).
#' @return A QQ plot.
#' @examples
#' p_values <- runif(1000)  # Simulate p-values
#' qq_plot(p_values)
#' @export
qq_plot <- function(p_values) {
  # Calculate expected and observed p-values
  observed <- -log10(sort(p_values))
  expected <- -log10(ppoints(length(p_values)))
  
  # Generate the plot
  ggplot2::ggplot(data.frame(expected = expected, observed = observed), ggplot2::aes(x = expected, y = observed)) +
    ggplot2::geom_point() +
    ggplot2::geom_abline(slope = 1, intercept = 0, color = "red") +
    ggplot2::labs(
      x = "Expected -log10(p-value)",
      y = "Observed -log10(p-value)",
      title = "QQ Plot"
    ) +
    ggplot2::theme_minimal()
}

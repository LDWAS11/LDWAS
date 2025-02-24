#' Generate a Manhattan Plot
#'
#' This function generates a Manhattan plot for GWAS results.
#'
#' @param p_values A numeric vector of p-values (1 x m).
#' @param snp_info A data frame containing SNP information (e.g., chromosome and position).
#' @param threshold A numeric value for the significance threshold (default: 5e-8).
#' @return A Manhattan plot.
#' @examples
#' p_values <- runif(1000)  # Simulate p-values
#' snp_info <- data.frame(chr = rep(1:10, each = 100), pos = 1:1000)  # Simulate SNP info
#' manhattan_plot(p_values, snp_info)
#' @export
manhattan_plot <- function(p_values, snp_info, threshold = 5e-8) {
  # Combine p-values and SNP information
  data <- data.frame(
    chr = snp_info$chr,
    pos = snp_info$pos,
    p = -log10(p_values)
  )
  
  # Add a cumulative position for plotting
  data <- data %>%
    dplyr::group_by(chr) %>%
    dplyr::mutate(cum_pos = cumsum(pos))
  
  # Generate the plot
  ggplot2::ggplot(data, ggplot2::aes(x = cum_pos, y = p, color = as.factor(chr))) +
    ggplot2::geom_point() +
    ggplot2::geom_hline(yintercept = -log10(threshold), linetype = "dashed", color = "red") +
    ggplot2::labs(
      x = "Chromosome",
      y = "-log10(p-value)",
      title = "Manhattan Plot"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "none")
}

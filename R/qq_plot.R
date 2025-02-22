qq_plot_gwas <- function(p_values) {
  if (is.null(p_values) || length(p_values) == 0) {
    stop("p_values must be a non-empty numeric vector.")
  }

  # Remove NA values
  p_values <- na.omit(p_values)

  # Number of p-values
  n <- length(p_values)

  # Expected p-values from a uniform(0,1) distribution
  p.uni <- -log10((1:n) / (n + 1))  # Theoretical quantiles

  # Observed p-values sorted in ascending order
  p.obs <- -log10(sort(p_values, decreasing = FALSE))  # Observed quantiles

  # Create the QQ-plot with hollow circles
  plot(p.uni, p.obs,
       xlab = "Expected -log10(P)",
       ylab = "Observed -log10(P)",
       main = "QQ-Plot of GWAS P-values",
       pch = 1, col = "black",  # pch = 1 makes hollow circles
       cex = 1.2)  # Adjust size if needed

  # Add the reference line (y = x)
  abline(a = 0, b = 1, col = "red", lwd = 2)
}





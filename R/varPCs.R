varPCs <- function(X, variance_threshold = 0.80, remove_NAs = FALSE) {
  # Ensure X is a matrix or data frame
  if (!is.matrix(X) && !is.data.frame(X)) {
    stop("X must be a matrix or data frame.")
  }

  # Ensure that variance_threshold is between 0 and 1
  if (variance_threshold <= 0 || variance_threshold >= 1) {
    stop("variance_threshold should be between 0 and 1.")
  }

  # Remove ID column if it's non-numeric (assuming it's the first column)
  if (!is.numeric(X[, 1])) {
    X <- X[, -1]  # Remove first column (ID column)
  }

  # Ensure all SNP columns are numeric and contain only 0, 1, 2, or NA
  if (!all(apply(X, 2, function(col) all(na.omit(col) %in% c(0, 1, 2))))) {
    stop("SNP columns must contain only 0, 1, or 2 (excluding NAs).")
  }

  # Handle missing values if remove_NAs is TRUE
  if (any(is.na(X))) {
    if (remove_NAs) {
      warning("Missing values in genotype (X) detected. Removing rows with NA values.")
      X <- X[complete.cases(X), ]
    } else {
      stop("Genotype matrix (X) contains missing values. Set remove_NAs = TRUE to remove them.")
    }
  }

  # Perform PCA on the genotype data
  pca_result <- prcomp(X)

  # Calculate the proportion of variance explained by each PC
  variance_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2)

  # Calculate the cumulative variance explained
  cumulative_variance <- cumsum(variance_explained)

  # Find the number of PCs that explain at least the given variance threshold
  num_PCs <- which(cumulative_variance >= variance_threshold)[1]

  return(num_PCs)
}

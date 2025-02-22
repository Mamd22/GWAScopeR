coVar_Filter <- function(X, C, num_PCs = 10, variance_threshold = 0.80, remove_NAs = FALSE) {
  # If num_PCs is provided as a numeric value, use it directly
  if (is.numeric(num_PCs)) {
    num_PCs <- num_PCs
  } else if (is.character(num_PCs) && num_PCs == "varPCs") {
    # If num_PCs is specified as "varPCs", use varPCs function to determine the number of PCs
    num_PCs <- varPCs(X, variance_threshold)
  }

  # Remove ID column from X if it's non-numeric (assuming it's the first column)
  if (!is.numeric(X[, 1])) {
    X <- X[, -1]  # Remove first column (ID column)
  }

  # Ensure X is numeric and check for missing values
  if (any(is.na(X))) {
    if (remove_NAs) {
      warning("Missing values in genotype (X) detected. Removing rows with NA values.")
      valid_indices <- complete.cases(X)
      X <- X[valid_indices, ]
      if (!is.null(C)) {
        C <- C[valid_indices, ]
      }
    } else {
      stop("Genotype data (X) contains missing values. Set remove_NAs = TRUE to remove them.")
    }
  }

  # Perform PCA on X (genotype matrix)
  pca_result <- prcomp(X)

  # Ensure num_PCs does not exceed the number of available PCs
  if (num_PCs > ncol(pca_result$x)) {
    warning("num_PCs exceeds the available number of PCs, using maximum available PCs.")
    num_PCs <- ncol(pca_result$x)
  }

  # Get the first 'num_PCs' principal components from PCA result
  PCs <- as.data.frame(pca_result$x[, 1:num_PCs])

  # Check if C (covariate matrix) has an ID column and remove it (if it's non-numeric)
  if (!is.null(C) && !is.numeric(C[, 1])) {
    C <- C[, -1]  # Remove first column (ID column)
  }

  # Ensure C is numeric and check for missing values
  if (any(is.na(C))) {
    if (remove_NAs) {
      warning("Missing values in covariates (C) detected. Removing rows with NA values.")
      valid_indices <- complete.cases(C)
      C <- C[valid_indices, ]
      X <- X[valid_indices, ]
    } else {
      stop("Covariate data (C) contains missing values. Set remove_NAs = TRUE to remove them.")
    }
  }

  # Combine PCs (from PCA) and existing covariates (C)
  C_combined <- cbind(C, PCs)

  # Remove linearly dependent PCs (if any)
  keep_PCs <- !apply(C_combined, 2, function(col) {
    tryCatch({
      lm(col ~ ., data = C_combined)  # Check linear dependence by running regression
      FALSE  # If no error, the column is not linearly dependent
    }, error = function(e) {
      TRUE  # If error occurs, column is linearly dependent
    })
  })

  # Filter out linearly dependent PCs and return the cleaned covariate matrix
  C_filtered <- C_combined[, keep_PCs]

  return(C_filtered)
}

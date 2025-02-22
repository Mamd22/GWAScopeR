GWAS_GLM <- function(y, X, SNP_names = NULL, C = NULL, return_full_results = FALSE) {
  # If y is a data frame, remove the first column if it's non-numeric (ID column)
  if (is.data.frame(y) && !is.numeric(y[, 1])) {
    y <- y[, -1]  # Remove ID column
  }

  # Ensure y is a numeric vector
  y <- as.numeric(y)

  # Check if the first column of X is non-numeric (ID column)
  if (!is.numeric(X[, 1])) {
    X <- X[, -1]  # Remove first column (ID column)
  }

  # Number of individuals and markers
  n <- length(y)
  m <- ncol(X)

  # Store p-values
  p_values <- numeric(m)
  snp_identifiers <- if (is.null(SNP_names)) colnames(X) else SNP_names

  # Loop over each SNP
  for (i in 1:m) {
    snp <- X[, i]  # Extract SNP column

    # If SNP has no variation (all values are the same), assign p-value 1
    if (length(unique(snp)) == 1) {
      p_values[i] <- 1
      next  # Skip further processing for this SNP
    }

    # Create model formula
    if (!is.null(C) && ncol(C) > 0) {
      # Ensure covariates (C) are included correctly in the formula
      data <- data.frame(y = y, snp = snp, C)
      formula <- as.formula(paste("y ~ snp +", paste(colnames(C), collapse = " + ")))
    } else {
      data <- data.frame(y = y, snp = snp)
      formula <- y ~ snp
    }

    # Fit GLM model
    model <- tryCatch({
      glm(formula, data = data, family = gaussian())
    }, error = function(e) {
      warning(paste("Error fitting model for SNP", i, ":", e$message))
      return(NULL)
    })

    if (!is.null(model)) {
      model_summary <- summary(model)

      # Check if "snp" is in the coefficients
      if ("snp" %in% rownames(model_summary$coefficients)) {
        p_values[i] <- model_summary$coefficients["snp", 4]  # Get the p-value
      } else {
        warning(paste("No SNP effect found for SNP", i))
        p_values[i] <- NA  # If SNP term is missing, assign NA
      }
    } else {
      p_values[i] <- NA  # If model fitting failed, assign NA
    }
  }

  # If return_full_results is TRUE, return SNPs and their p-values
  if (return_full_results) {
    gwas_results <- data.frame(
      SNP = snp_identifiers,
      P_Value = p_values
    )
    return(gwas_results)
  } else {
    # If return_full_results is FALSE, return only the p-value vector
    return(p_values)
  }
}

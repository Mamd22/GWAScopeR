manhattan_plot <- function(map_data, pval_vector = NULL, gwas_data = NULL,
                           threshold_type = c("fixed"), alpha = 0.01,
                           show_threshold_labels = TRUE, num_qtns = 10,
                           user_qtn_positions = NULL) {

  # If a vector of P-values is given
  if (!is.null(pval_vector)) {
    if (length(pval_vector) != nrow(map_data)) {
      stop("The length of pval_vector must match the number of SNPs in map_data.")
    }
    gwas_data <- data.frame(SNP = map_data$SNP, P = as.numeric(pval_vector))
  }

  # If gwas_data is provided, handle different column names for P-values
  if (!is.null(gwas_data)) {
    if ("P_Value" %in% colnames(gwas_data)) {
      gwas_data <- gwas_data %>% rename(P = P_Value)
    }

    # Ensure P-values are numeric
    gwas_data$P <- as.numeric(gwas_data$P)

    # Remove NAs
    gwas_data <- na.omit(gwas_data)

    # Merge with SNP mapping data
    data <- merge(gwas_data, map_data, by = "SNP")
  } else {
    stop("gwas_data is required or pval_vector must be provided.")
  }

  # Convert Chromosome to a properly ordered factor
  data <- data %>%
    mutate(Chromosome = factor(Chromosome, levels = sort(unique(as.numeric(as.character(Chromosome)))))) %>%
    arrange(Chromosome, Position)

  # Convert P-values to -log10(P) for plotting
  data$logP <- -log10(data$P)

  # Calculate genome-wide thresholds based on user input
  num_snps <- nrow(data)
  thresholds <- list()  # Store all threshold values
  threshold_labels <- list()  # Store threshold labels for the legend

  if ("bonferroni" %in% threshold_type) {
    thresholds$bonferroni <- alpha / num_snps
    threshold_labels$bonferroni <- paste("Bonferroni: ", round(-log10(thresholds$bonferroni), 2))
  }
  if ("inflation" %in% threshold_type) {
    chisq_vals <- qchisq(1 - data$P, df = 1)  # Compute chi-square values
    lambda <- median(chisq_vals) / qchisq(0.5, df = 1)  # Compute genomic inflation factor
    adjusted_pvals <- pchisq(chisq_vals / lambda, df = 1, lower.tail = FALSE)
    thresholds$inflation <- max(adjusted_pvals[adjusted_pvals < alpha], na.rm = TRUE)  # Inflation threshold
    threshold_labels$inflation <- paste("Inflation: ", round(-log10(thresholds$inflation), 2))
  }
  if ("fdr" %in% threshold_type) {
    fdr_vals <- p.adjust(data$P, method = "BH")
    thresholds$fdr <- max(data$P[fdr_vals < alpha], na.rm = TRUE)  # Maximum significant p-value under FDR control
    threshold_labels$fdr <- paste("FDR: ", round(-log10(thresholds$fdr), 2))
  }
  if ("fixed" %in% threshold_type) {
    thresholds$fixed <- alpha
    threshold_labels$fixed <- paste("Fixed: ", round(-log10(thresholds$fixed), 2))
  }

  # Assign X-axis positions
  data$Index <- seq_len(nrow(data))

  # Define colors for alternating chromosomes
  colors <- rep(c("black", "gray60"), length.out = length(unique(data$Chromosome)))

  # Identify the top 'num_qtns' SNPs based on the highest p-values (smallest -log10(p))
  top_qtns <- data[order(data$logP, decreasing = TRUE),][1:num_qtns, ]

  # Plot using ggplot
  plot <- ggplot(data, aes(x = Index, y = logP, color = Chromosome)) +
    geom_point(size = 1.5) +
    scale_color_manual(values = colors)  # Color for chromosomes

  # Add threshold lines for each threshold type
  for (threshold_name in names(thresholds)) {
    threshold_value <- thresholds[[threshold_name]]
    threshold_color <- ifelse(threshold_name == "fixed", "blue",
                              ifelse(threshold_name == "bonferroni", "green3",
                                     ifelse(threshold_name == "fdr", "red",
                                            ifelse(threshold_name == "inflation", "purple", "black"))))

    plot <- plot + geom_hline(yintercept = -log10(threshold_value),
                              linetype = "dashed",
                              color = threshold_color)  # Color for threshold lines
  }

  # Add vertical dashed lines for top QTN positions
  plot <- plot + geom_vline(xintercept = top_qtns$Index, linetype = "dashed", color = "deepskyblue")

  # Add vertical dashed lines for user-provided QTN positions (if any)
  if (!is.null(user_qtn_positions)) {
    user_qtn_positions <- user_qtn_positions[user_qtn_positions %in% data$Index]
    plot <- plot + geom_vline(xintercept = user_qtn_positions, linetype = "dashed", color = "orange")
  }

  # Add threshold labels if requested
  if (show_threshold_labels) {
    y_offset <- 0  # Start with no vertical offset
    for (threshold_name in names(thresholds)) {
      threshold_value <- thresholds[[threshold_name]]
      label <- threshold_labels[[threshold_name]]
      threshold_color <- ifelse(threshold_name == "fixed", "blue",
                                ifelse(threshold_name == "bonferroni", "green",
                                       ifelse(threshold_name == "fdr", "red",
                                              ifelse(threshold_name == "inflation", "purple", "black"))))

      # Apply vertical offset if the labels overlap
      plot <- plot + annotate("text", x = nrow(data) * 0.95,
                              y = -log10(threshold_value) + y_offset,
                              label = label, color = threshold_color, hjust = 0)

      # Adjust the offset for the next label
      y_offset <- y_offset + 0.3  # Increase the offset for subsequent labels
    }
  }

  # Customize the plot
  plot <- plot +
    labs(x = "Chromosome", y = "-log10(P-value)", title = "Manhattan Plot") +
    theme_minimal() +
    theme(legend.position = "none") +
    scale_x_continuous(breaks = data %>%
                         group_by(Chromosome) %>%
                         summarize(mid = mean(Index)) %>%
                         pull(mid),
                       labels = levels(data$Chromosome))

  print(plot)
}

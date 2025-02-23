\name{manhattan_plot}
\alias{manhattan_plot}
\title{Generate a Manhattan Plot for GWAS Data}
\description{
  This function generates a Manhattan plot for visualizing the results of a Genome-Wide Association Study (GWAS).
  The plot displays the -log10 of the p-values for each SNP, with horizontal threshold lines indicating significance levels.
  Vertical lines can be added to mark the top Quantitative Trait Nucleotides (QTNs) or user-provided QTN positions.
}
\usage{
manhattan_plot(map_data, pval_vector = NULL, gwas_data = NULL,
                threshold_type = c("fixed"), alpha = 0.01,
                show_threshold_labels = TRUE, num_qtns = 10,
                user_qtn_positions = NULL)
}
\arguments{
  \item{map_data}{
    A data frame containing SNP mapping information with columns `SNP`, `Chromosome`, and `Position`.
    The rows correspond to individual SNPs.
  }
  \item{pval_vector}{
    A numeric vector of p-values corresponding to the SNPs in `map_data`. The length of `pval_vector` must
    match the number of rows in `map_data`. Default is `NULL`.
  }
  \item{gwas_data}{
    A data frame containing GWAS results, including a column with SNP identifiers and p-values. If provided,
    this will be used instead of `pval_vector`. The column name for p-values should be `P_Value`. Default is `NULL`.
  }
  \item{threshold_type}{
    A character vector specifying the type of threshold to use for significance lines. Available options:
    "fixed", "bonferroni", "inflation", "fdr". Default is `c("fixed")`.
  }
  \item{alpha}{
    A numeric value representing the significance threshold for the p-values. Default is `0.01`.
  }
  \item{show_threshold_labels}{
    A logical value indicating whether threshold labels should be shown on the plot. Default is `TRUE`.
  }
  \item{num_qtns}{
    A numeric value indicating the number of top Quantitative Trait Nucleotides (QTNs) to highlight on the plot.
    Default is `10`.
  }
  \item{user_qtn_positions}{
    A numeric vector of user-defined QTN positions to highlight on the plot. Default is `NULL`.
  }
}
\details{
  This function generates a Manhattan plot where the x-axis represents the SNPs along the chromosomes, and the y-axis
  represents the -log10 of the p-values. The function allows you to specify multiple threshold types to display
  significance lines, including Bonferroni correction, genomic inflation, false discovery rate (FDR), and a fixed alpha value.
  The plot can also highlight the top SNPs based on their p-values or user-defined positions.
}
\value{
  A `ggplot2` object representing the Manhattan plot.
}
\references{

}
\author{
  Michelle Marcano-Delgado
}
\note{

}
\seealso{
  \code{\link{varPCs}}, \code{\link{prcomp}}, \code{\link{p.adjust}}, \code{\link{qchisq}}
}
\examples{
# Load required packages
library(dplyr)
library(ggplot2)

# Simulated data
set.seed(123)
y <- rnorm(100)  # Phenotype values
X <- matrix(sample(c(0, 1, 2), 100 * 500, replace = TRUE), ncol = 500)  # Genotype matrix
C <- data.frame(
  ID = 1:100,  # Unique identifier for each individual
  Cov1 = rnorm(100),
  Cov2 = rnorm(100),
  Cov3 = rnorm(100),
  Cov4 = rnorm(100),
  Cov5 = rnorm(100)
)

map_data <- data.frame(
  SNP = paste0("SNP", 1:500),                   # SNP identifiers
  Chromosome = sample(1:10, 500, replace = TRUE),  # Randomly assign SNPs to 10 chromosomes
  Position = sample(1:1e6, 500, replace = TRUE)    # Random positions in base pairs
)

# Ensure GWAS_GLM function is available
if (exists("GWAS_GLM")) {
  p_values <- GWAS_GLM(y, X, C=C)  # Run GWAS

  # Generate Manhattan plot
  manhattan_plot(map_data, pval_vector = p_values, threshold_type = c("fixed", "bonferroni"), alpha = 0.01)
} else {
  warning("GWAS_GLM function not found. Ensure GWAS_GLM is correctly defined before running this example.")
}
}

\keyword{GWAS}
\keyword{Manhattan Plot}

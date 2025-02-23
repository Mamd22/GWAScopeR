\name{GWAS_GLM}
\alias{GWAS_GLM}
\title{Perform Genome-Wide Association Study (GWAS) using GLM}
\description{
This function performs a Genome-Wide Association Study (GWAS) using a Generalized Linear Model (GLM).
It tests the association between SNPs and a quantitative phenotype, optionally adjusting for covariates.
}
\usage{
GWAS_GLM(y, X, SNP_names = NULL, C = NULL, return_full_results = FALSE)
}
\arguments{
  \item{y}{
    A numeric vector or data frame of phenotype values. If a data frame, the first column is assumed to be an ID column and will be removed if non-numeric.
  }
  \item{X}{
    A numeric matrix or data frame of SNP genotypes, where rows represent individuals and columns represent SNPs. The first column is removed if it is non-numeric (assumed ID column).
  }
  \item{SNP_names}{
    A character vector of SNP names. If NULL, column names of X are used.
  }
  \item{C}{
    A numeric matrix or data frame of covariates (e.g., principal components). Default is NULL (no covariates).
  }
  \item{return_full_results}{
    Logical. If TRUE, returns a data frame with SNP names and p-values. If FALSE, returns only a numeric vector of p-values. Default is FALSE.
  }
}
\details{
The function fits a Generalized Linear Model (GLM) with a Gaussian family for each SNP, adjusting for covariates if provided. SNPs with no variation are assigned a p-value of 1. If model fitting fails, the p-value is set to NA.
}
\value{
  If `return_full_results = TRUE`, returns a data frame with columns:
  \itemize{
    \item \code{SNP} - SNP identifier.
    \item \code{P_Value} - P-value for SNP association.
  }
  If `return_full_results = FALSE`, returns a numeric vector of p-values.
}
\references{
}
\author{
  Michelle Marcano-Delgado
}
\note{
  This function assumes the phenotype data and SNP data are numeric and that covariates are optional.
}

\seealso{
  \code{\link{}}
}

\examples{
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
  Cov5 = rnorm(100))

# Run GWAS
p_values <- GWAS_GLM(y, X, C=C)
head(p_values)

# Get full results with SNP names
snp_names <- paste0("SNP", 1:500)
results <- GWAS_GLM(y, X, SNP_names = snp_names, C = C, return_full_results = TRUE)
head(results)
}

\keyword{GWAS}
\keyword{GLM}
\keyword{genetic association}

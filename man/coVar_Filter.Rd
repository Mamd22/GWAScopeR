\name{coVar_Filter}
\alias{coVar_Filter}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Filter Covariates by Removing Linearly Dependent Principal Components
}
\description{
This function filters covariates by removing linearly dependent principal components (PCs) from a genotype matrix (X)
and a covariate matrix (C). The number of PCs to consider can either be specified directly or determined using the `varPCs` function.
It performs Principal Component Analysis (PCA) on the genotype data and checks for linear dependence in the combined matrix
of covariates and PCs.
}
\usage{
coVar_Filter(X, C, num_PCs = 10, variance_threshold = 0.80, remove_NAs = FALSE)
}
\arguments{
  \item{X}{
    A numeric matrix or data frame representing the genotype data (SNPs). The rows represent individuals, and columns represent SNPs.
    The first column, if non-numeric, will be removed assuming it's an ID column.
  }
  \item{C}{
    A numeric matrix or data frame representing covariates. The first column, if non-numeric, will be removed assuming it's an ID column.
  }
  \item{num_PCs}{
    Either a numeric value indicating the number of principal components to retain or a character string `"varPCs"`.
    If `"varPCs"`, the function uses the `varPCs` function to determine the number of PCs that explain the specified variance threshold.
    Default is 10.
  }
  \item{variance_threshold}{
    A numeric value between 0 and 1 indicating the proportion of the variance to explain using the principal components.
    Only used if `num_PCs = "varPCs"`. Default is 0.80 (80% of the variance).
  }
  \item{remove_NAs}{
    A logical value. If `TRUE`, rows containing missing values (NA) in either `X` or `C` will be removed.
    If `FALSE`, the function will stop if missing values are detected. Default is `FALSE`.
  }
}
\details{
This function combines the PCs obtained from the genotype data (`X`) with the covariate matrix (`C`), then removes
linearly dependent columns using linear regression to detect dependence. The returned matrix contains only the covariates
and PCs that are not linearly dependent on each other.
}
\value{
A numeric matrix or data frame of covariates (C) with linearly dependent PCs removed. The returned matrix includes both the original covariates and the retained PCs.
}
\references{
See the `varPCs` function for more details on determining the number of principal components to retain.
}
\author{
Michelle Marcano-Delgado
}
\note{
This function is useful for filtering out PCs that may not provide independent variance, ensuring that the remaining PCs are orthogonal to each other.
}

\seealso{
\code{\link{varPCs}}, \code{\link{prcomp}}
}
\examples{
## Example 1: Filter covariates by removing linearly dependent PCs (keeping top 10 PCs)
filtered_covariates <- coVar_Filter(genotype_matrix, covariate_matrix, num_PCs = 10)

## Example 2: Automatically determine the number of PCs to retain based on the variance threshold
filtered_covariates <- coVar_Filter(genotype_matrix, covariate_matrix, num_PCs = "varPCs", variance_threshold = 0.80)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.

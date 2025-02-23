\name{varPCs}
\alias{varPCs}
\title{Determine the Number of Principal Components (PCs) Explaining a Given Variance Threshold}
\description{
  This function performs Principal Component Analysis (PCA) on a genotype matrix (X) and determines the
  number of principal components that explain at least a specified proportion of the variance. The function
  handles missing values based on the `remove_NAs` parameter.
}
\usage{
varPCs(X, variance_threshold = 0.80, remove_NAs = FALSE)
}
\arguments{
  \item{X}{
    A numeric matrix or data frame where columns represent SNPs (genotypes) and rows represent individuals.
    The values in the matrix should be coded as 0, 1, or 2 (for different genotypes) with potential NA values.
    The first column, if non-numeric, will be removed assuming it's an ID column.
  }
  \item{variance_threshold}{
    A numeric value between 0 and 1 indicating the proportion of the variance to explain using the principal components.
    Default is 0.80 (80% of the variance).
  }
  \item{remove_NAs}{
    A logical value. If `TRUE`, rows containing missing values (NA) will be removed from the matrix `X`.
    If `FALSE`, the function will stop if missing values are detected. Default is `FALSE`.
  }
}
\details{
  This function performs PCA on the genotype data matrix `X` and calculates the proportion of variance explained
  by each principal component. It then returns the number of PCs that cumulatively explain the variance up to the
  `variance_threshold` specified.
}
\value{
  An integer indicating the number of principal components (PCs) that explain at least the specified variance threshold.
}
\references{
  \code{\link{prcomp}} for performing PCA.
}
\author{
  Michelle Marcano-Delgado
}
\note{
  The input matrix `X` must have columns representing SNPs and rows representing individuals. The first column, if non-numeric, is assumed to be an ID column and is removed. Missing values can be handled based on the `remove_NAs` parameter.
}
\seealso{
  \code{\link{prcomp}}, \code{\link{complete.cases}}
}
\examples{
  # Example 1: Calculate the number of PCs that explain at least 80% of the variance
  # Assuming 'genotype_matrix' is a data frame or matrix of SNP genotypes
   set.seed(123)
   genotype_matrix <- matrix(sample(c(0, 1, 2), 100 * 500, replace = TRUE), ncol = 500)  # Genotype matrix
  num_PCs <- varPCs(genotype_matrix, variance_threshold = 0.80)
  num_PCs

  # Example 2: Remove NAs before calculating the number of PCs
  num_PCs_no_NAs <- varPCs(genotype_matrix, variance_threshold = 0.80, remove_NAs = TRUE)
}
\keyword{PCA}
\keyword{Principal Components}
\keyword{Variance Threshold}

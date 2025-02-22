\name{plot_cumulative_variance}
\alias{plot_cumulative_variance}
\title{Plot Cumulative Variance Explained by Principal Components}
\description{
  This function performs Principal Component Analysis (PCA) on a genotype matrix and generates a plot showing the cumulative variance explained by the principal components.
  The plot includes a dashed line representing the given variance threshold, and the number of principal components (PCs) required to exceed this threshold is annotated.
}
\usage{
plot_cumulative_variance(X, variance_threshold = 0.80, remove_NAs = FALSE)
}
\arguments{
  \item{X}{
    A matrix or data frame containing genotype data, where rows represent samples and columns represent SNPs.
    The first column is assumed to be an ID column (if non-numeric, it will be removed), and all SNP columns should contain numeric values 0, 1, or 2, with optional NA values.
  }
  \item{variance_threshold}{
    A numeric value between 0 and 1 specifying the cumulative variance threshold for determining the number of PCs to retain. Default is 0.80.
  }
  \item{remove_NAs}{
    A logical value indicating whether to remove rows with missing values (NA). If `TRUE`, rows with missing values will be removed. Default is `FALSE`.
  }
}
\details{
  This function performs PCA on the input genotype matrix and calculates the cumulative variance explained by the principal components.
  A plot is generated with the cumulative variance on the y-axis and the number of principal components on the x-axis.
  A dashed line is drawn at the specified variance threshold, and the number of principal components required to reach the threshold is highlighted.
}
\value{
  A numeric value indicating the number of principal components needed to explain at least the specified variance threshold.
}
\references{

}
\author{
  Michelle Marcano-Delgado
}
\note{
  If there are missing values in the genotype matrix and `remove_NAs` is set to `FALSE`, the function will stop with an error.
}
\seealso{
  \code{\link{prcomp}}, \code{\link{complete.cases}}, \code{\link{ggplot2}}
}
\examples{
  # Example 1: Plot cumulative variance with a variance threshold of 80%
  plot_cumulative_variance(X, variance_threshold = 0.80)

  # Example 2: Remove rows with NA values and use a 90% variance threshold
  plot_cumulative_variance(X, variance_threshold = 0.90, remove_NAs = TRUE)
}
\keyword{PCA}
\keyword{Variance}
\keyword{Genotype Data}

\name{qq_plot_gwas}
\alias{qq_plot_gwas}
\title{QQ-Plot for GWAS P-values}
\description{
  This function generates a QQ-plot for GWAS p-values, comparing the observed \eqn{-log10(p)} values with the expected values under a uniform(0,1) distribution.
}
\usage{
qq_plot_gwas(p_values)
}
\arguments{
  \item{p_values}{
    A numeric vector of p-values. Must be non-empty and contain valid numeric values.
  }
}
\details{
  This function creates a QQ-plot by plotting the observed \eqn{-log10(p)} values against the expected \eqn{-log10(p)} values assuming a uniform distribution. A red reference line (y = x) is added to help assess the deviation of observed values from the expected distribution.
}
\value{
  The function generates and displays a QQ-plot showing the comparison between observed and expected \eqn{-log10(p)} values.
}
\references{

}
\author{
  Michelle Marcano-Delgado
}
\note{
  The input vector of p-values must be non-empty and should contain valid numeric values. Missing values will be removed automatically.
}
\seealso{
  \code{\link{plot}}, \code{\link{sort}}, \code{\link{abline}}
}
\examples{
  # Example 1: Simulated p-values and generate a QQ-plot
  set.seed(123)
  p_vals <- runif(1000, min = 1e-8, max = 1)  # Simulated p-values
  qq_plot_gwas(p_vals)
}
\keyword{GWAS}
\keyword{QQ-Plot}
\keyword{P-values}

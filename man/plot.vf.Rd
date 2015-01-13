\name{plot.vf}
\alias{plot.vf}
\title{ Plots fitted vectors onto an ordination diagram }
\description{
Add vector fitting arrows to an existing ordination plot.
}
\usage{
\method{plot}{vf}(x, pval = 1, cex = 0.8, ascale = 0.9, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ an object of S3 class \code{vf}, created by \code{vf()} }
  \item{pval}{ optional, critical p-value for choosing variables to plot }
  \item{cex}{ text size }
  \item{ascale}{ optional, proportion of plot area to use when calculating arrow length }
  \item{\dots}{ optional, other graphics parameters }
}
\value{
  Adds arrows to an existing ordination plot. Only arrows with a p-value less than pval are added. By default, all variables are shown.
  }
\author{ Sarah Goslee, Sarah.Goslee@ars.usda.gov }

\seealso{ \code{\link{vf}}}
\examples{

\dontrun{
# Example of multivariate analysis using built-in iris dataset
data(iris)
iris <- iris[seq(1, 150, by=3),]
iris.md <- distance(iris[,1:4], "mahal")

# Minimum-stress 2-dimensional nonmetric multidimensional scaling configuration
# Uses small number of separate ordinations (5) to increase speed of example.
# Use more for final analysis.
iris.nmds <- nmds(iris.md, mindim=2, maxdim=2, nits=3)
iris.nmin <- nmds.min(iris.nmds)

# Plot NMDS result with symbols denoting species
plot(iris.nmin, pch=as.numeric(iris[,5]))

# Fit vectors for the main variables to the NMDS configuration
iris.vf <- vf(iris.nmin, iris[,1:4], nperm=10)
plot(iris.vf, col="blue")
}

}
\keyword{ multivariate }
\keyword{ hplot }

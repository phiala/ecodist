\name{nmds.min}
\alias{nmds.min}
\title{ Find minimum stress configuration }
\description{
  Finds minimum stress configuration from output of nmds()
}
\usage{
nmds.min(x, dims = 0)
}
\arguments{
  \item{x}{ output from nmds() }
  \item{dims}{ desired dimensionality of result. If dims == 0 (the default) then the overall minimum stress configuration is chosen. }
}
\value{
  Configuration of minimum stress ordination (dataframe of
coordinates). The stress and r^2 for the minimum stress configuration are stored as attributes.
}
\author{ Sarah Goslee, Sarah.Goslee@ars.usda.gov }

\seealso{ \code{\link{nmds}} }
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

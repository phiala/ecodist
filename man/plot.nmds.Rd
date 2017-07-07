\name{plot.nmds}
\alias{plot.nmds}
\title{ Plot information about NMDS ordination }
\description{
Graphical display of stress and r2 for NMDS ordination along number of dimensions.
}
\usage{
\method{plot}{nmds}(x, plot = TRUE, xlab = "Dimension", ...)
}

\arguments{
  \item{x}{ an object of S3 class \code{nmds}, created by \code{nmds()} }
  \item{plot}{ optional, if TRUE a figure is produced }
  \item{xlab}{ optional, label for x axis of graph }
  \item{\dots}{ optional, other graphics parameters }
}
\value{
  Produces a two-panel plot with stress and r2 for ordination by number of dimensions. Points show the mean value; lines indicate minimum and maximum.
  }
\author{ Sarah Goslee, Sarah.Goslee@ars.usda.gov }

\seealso{ \code{\link{nmds}}}
\examples{

data(iris)
iris.d <- dist(iris[,1:4])

### nmds() is timeconsuming, so this was generated
### in advance and saved.
### set.seed(1234)
### iris.nmds <- nmds(iris.d, nits=20, mindim=1, maxdim=4)
### save(iris.nmds, file="ecodist/data/iris.nmds.rda")
data(iris.nmds)

# examine fit by number of dimensions
plot(iris.nmds)

# choose the best two-dimensional solution to work with
iris.nmin <- min(iris.nmds, dims=2)

}
\keyword{ multivariate }
\keyword{ hplot }

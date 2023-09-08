\name{plot.vf}
\alias{plot.vf}
\title{ Plots fitted vectors onto an ordination diagram }
\description{
Add vector fitting arrows to an existing ordination plot.
}
\usage{
\method{plot}{vf}(x, pval = NULL, r = NULL, cex = 0.8, ascale = 0.9, \dots)
}

\arguments{
  \item{x}{ an object of S3 class \code{vf}, created by \code{vf()} }
  \item{pval}{ optional, critical p-value for choosing variables to plot }
  \item{r}{ optional, minimum Mantel r for choosing variables to plot }
  \item{cex}{ text size }
  \item{ascale}{ optional, proportion of plot area to use when calculating arrow length }
  \item{\dots}{ optional, other graphics parameters }
}
\value{
  Adds arrows to an existing ordination plot. Only arrows with a p-value less than pval are added. By default, all variables are shown.
  }
\author{ Sarah Goslee }

\seealso{ \code{\link{vf}}}
\examples{

# Example of multivariate analysis using built-in iris dataset
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

# fit the data to the ordination as vectors
### vf() is timeconsuming, so this was generated
### in advance and saved.
### set.seed(1234)
### iris.vf <- vf(iris.nmin, iris[,1:4], nperm=1000)
### save(iris.vf, file="ecodist/data/iris.vf.rda")
data(iris.vf)
plot(iris.nmin, col=as.numeric(iris$Species), pch=as.numeric(iris$Species), main="NMDS")
plot(iris.vf)

}
\keyword{ multivariate }
\keyword{ hplot }

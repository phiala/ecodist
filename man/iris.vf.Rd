\name{iris.vf}
\alias{iris.vf}
\docType{data}
\title{Example for vector fitting on ordination}
\description{
  An object of class vf for use in the examples for \code{\link{nmds}} and  \code{\link{vf}}. Many of the functions in \code{ecodist} take a long time to run, so prepared examples have been included.
}
\usage{data(iris.vf)}
\format{
See \code{\link{vf}} for current format specification.
}
\author{ Sarah Goslee }
\seealso{ \code{\link{nmds}},  \code{\link{vf}} }
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

# rotate the configuration to maximize variance
iris.rot <- princomp(iris.nmin)$scores

# rotation preserves distance apart in ordination space
cor(dist(iris.nmin), dist(iris.rot))

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
\keyword{datasets}

\name{rotate2d}
\alias{rotate2d}
\title{
  Rotate a 2D ordination.
}
\description{
  Rotates a two-dimensional ordination configuration to place the direction indicated along the horizontal axis.
}
\usage{
rotate2d(ord, x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ord}{
    A matrix or data frame with two columns, or a vf object, containing the points of an ordination configuration.
}
  \item{x}{
    The coordinates of a point in the ordination space. }
}
\details{
    The configuration ord is rotated so that the vector defined by c(0, 0), and x is along the horizontal axis. This can be useful for placing a specific variable, for instance from vf(), in a consistent direction across multiple ordinations. Doing so can facilitate interpretation.
}
\value{
    A rotated data frame of coordinates of the same size as ord and in the same order. If ord was produced by vf(), the complete vf object is returned.
}
\author{
Sarah Goslee
}

\seealso{
    \code{\link{vf}}, \code{\link{nmds}}
}
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

# rotate configuration so Sepal Width is along the horizontal axis

iris.nmin.rot <- rotate2d(iris.nmin, iris.vf[2, 1:2])
iris.vf.rot <- rotate2d(iris.vf, iris.vf[2, 1:2])

plot(iris.nmin.rot, col=as.numeric(iris$Species), pch=as.numeric(iris$Species), main="NMDS")
plot(iris.vf.rot)

}
\keyword{ multivariate }


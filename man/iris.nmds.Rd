\name{iris.nmds}
\alias{iris.nmds}
\docType{data}
\title{Example for nmds}
\description{
  An object of class nmds for use in the example for \code{\link{nmds}}. Many of the functions in \code{ecodist} take a long time to run, so prepared examples have been included.
}
\usage{data(iris.nmds)}


\format{
See \code{\link{nmds}} for current format specification.
}


\seealso{ \code{\link{nmds}} }

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


}
\keyword{datasets}

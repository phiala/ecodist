\name{min.nmds}
\alias{min.nmds}
\alias{nmds.min}
\title{ Find minimum stress configuration }
\description{
  Finds minimum stress configuration from output of nmds()
}
\usage{
\method{min}{nmds}(..., na.rm = FALSE, dims = 2)
nmds.min(x, dims = 2)
}
\arguments{
  \item{...}{ output from nmds() }
  \item{x}{ output from nmds() }
  \item{dims}{ desired dimensionality of result. If dims = 0 then the overall minimum stress configuration is chosen. By default, the best two-dimensional configuration is returned. }
  \item{na.rm}{Not used; there should be no NA values in a NMDS configuration.}
}
\value{
  Configuration of minimum stress ordination (dataframe of
coordinates). The stress and r^2 for the minimum stress configuration are stored as attributes.
}
\details{
	For back-compatibility, the \code{nmds.min} function remains.
}
\author{ Sarah Goslee }

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

# choose the best two-dimensional solution to work with
iris.nmin <- min(iris.nmds, dims=2)

}
\keyword{ multivariate }

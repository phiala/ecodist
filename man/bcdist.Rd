\name{bcdist}
\alias{bcdist}
\title{ Bray-Curtis distance }
\description{
  Returns the Bray-Curtis (also known as Sorenson, 1 - percent similarity) pairwise distances for the objects in the data. It is duplicated by functionality within \code{\link{distance}} but remains for backward compatibility and because it is substantially faster.
}
\usage{
bcdist(x, rmzero = FALSE)
}
\arguments{
  \item{x}{ matrix or data frame with rows as samples and columns as variables (such as 
species). Distances will be calculated for each pair of rows. }
  \item{rmzero}{ If rmzero=TRUE, empty rows will be removed from the data before 
distances are calculated. Otherwise, the distance between two empty rows is assumed to be 0 (the default). }
}
\value{
  This function returns a column-order lower-triangular distance matrix. The returned object has an attribute, Size, giving the number of objects, that is, nrow(x). The length of the vector that is returned is nrow(x)*(nrow(x)-1)/2. 
  }
\author{ Sarah Goslee, Sarah.Goslee@ars.usda.gov }

\seealso{ \code{\link{dist}},  \code{\link{distance}}}
\examples{
data(iris)
iris.bc <- bcdist(iris[,1:4])
# equivalent to
iris.bc2 <- distance(iris[,1:4], "bray-curtis")
}
\keyword{ multivariate }

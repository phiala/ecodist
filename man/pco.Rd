\name{pco}
\alias{pco}
\title{ Principal coordinates analysis }
\description{
  Principal coordinates analysis (classical scaling).
}
\usage{
pco(x, negvals = "zero", dround = 0)
}
\arguments{
  \item{x}{ a lower-triangular dissimilarity matrix. }
  \item{negvals}{ if = "zero" sets all negative eigenvalues to zero;
if = "rm" corrects for negative eigenvalues using method
1 of Legendre and Anderson 1999. }
  \item{dround}{ if greater than 0, attempts to correct for round-off error by
rounding to that number of places. }
}
\details{
  PCO (classical scaling, metric multidimensional scaling) is very similar to principal components analysis, but allows the use of any dissimilarity metric.
}
\value{
  \item{values }{eigenvalue for each component. This is a measure of the variance explained by each dimension.}
  \item{vectors }{eigenvectors. Each column contains the scores for that dimension.}
}
\author{ Sarah Goslee }

\seealso{ \code{\link{princomp}}, \code{\link{nmds}} }
\examples{
data(iris)
iris.d <- dist(iris[,1:4])
iris.pco <- pco(iris.d)

# scatterplot of the first two dimensions
plot(iris.pco$vectors[,1:2], col=as.numeric(iris$Species),
  pch=as.numeric(iris$Species), main="PCO", xlab="PCO 1", ylab="PCO 2")
}
\keyword{ multivariate }

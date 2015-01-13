\name{pmgram}
\alias{pmgram}
\title{ Partial Mantel correlogram }
\description{
 This function calculates simple and partial multivariate correlograms. 
}
\usage{
pmgram(data, space, partial, breaks, nclass, stepsize, resids = FALSE, nperm = 1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ lower-triangular dissimilarity matrix. This can be either an object of class dist (treated as one column) or a matrix or data frame with one or two columns, each of which is an independent lower-triangular dissimilarity in vector form. }
  \item{space}{ lower-triangular matrix of geographic distances. }
  \item{partial}{ optional, lower-triangular dissimilarity matrix of ancillary data. }
  \item{breaks}{ locations of class breaks. If specified, overrides nclass and stepsize. }
  \item{nclass}{ number of distance classes. If not specified, Sturge's rule will be used
to determine an appropriate number of classes. }
  \item{stepsize}{ width of each distance class. If not specified, nclass and the range of space.d will be used to calculate an appropriate default. }
  \item{resids}{ if resids=TRUE, will return the residuals for each distance class. Otherwise returns 0. }
  \item{nperm}{ number of permutations to use. If set to 0, the permutation test will be omitted. }
}
\details{
This function does four different analyses:
	If data has 1 column and partial is missing, 
	calculates a multivariate correlogram for data.

If data has 2 columns and partial is missing,
	calculates Mantel r between the two columns for each distance class.

If data has 1 column and partial exists,
	does a multivariate correlogram for the residuals,
	taking residuals over whole extent.

If data has 2 columns and partial exists,
	does a partial multivariate correlogram, 
	calculating partial for each distance class separately.

The statistic used is not the standard Mantel r. For one variable, using Euclidean distance, this metric converges on the familiar Moran au-
tocorrelation. Like the Moran autocorrelation function, this statistic usually falls between -1 and 1, but is not bounded by those limits. Unlike the Moran function, this correlogram can be used for multivariate data, and can be extended to partial tests.

}
\value{
  Returns a object of class mgram, which is a list containing two objects:
  mgram is a matrix with one row for each distance class and 4 columns:
  \item{lag }{midpoint of the distance class.}
  \item{ngroup }{number of distances in that class.}
  \item{piecer }{Mantel r value or appropriate statistic (see Details).}
  \item{pval }{two-sided p-value.}
  resids is a vector of the residuals (if calcuated) and can be accessed with the residuals() method.
}
\author{ Sarah Goslee, Sarah.Goslee@ars.usda.gov }

\seealso{ \code{\link{mgram}}, \code{\link{mantel}} }
\examples{

\dontrun{
# generate a simple surface
x <- matrix(1:10, nrow=10, ncol=10, byrow=FALSE)
y <- matrix(1:10, nrow=10, ncol=10, byrow=TRUE)
z1 <- x + 3*y
z2 <- 2*x - y

# look at patterns
par(mfrow=c(1,2))
image(z1)
image(z2)

# analyze the pattern of z across space
z1 <- as.vector(z1)
z2 <- as.vector(z2)
z1.d <- distance(z1, "eucl")
z2.d <- distance(z2, "eucl")

space <- cbind(as.vector(x), as.vector(y))
space.d <- distance(space, "eucl")

# take partial correlogram of z2 on the residuals of z1 ~ space.d
z.pmgram <- pmgram(z1.d, space.d, z2.d, nperm=0)
par(mfrow=c(1,1))
plot(z.pmgram, pval=0.1)
}

}
\keyword{ multivariate }

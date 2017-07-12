\name{pmgram}
\alias{pmgram}
\title{ Partial Mantel correlogram }
\description{
 This function calculates simple and partial multivariate correlograms.
}
\usage{
pmgram(data, space, partial, breaks, nclass, stepsize, resids = FALSE, nperm = 1000)
}

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
	calculates a Mantel cross-correlogram, calculating the Mantel r between
	the two columns for each distance class separately.

If data has 1 column and partial exists,
	calculates a partial multivariate correlogram based on residuals of data ~ partial.

If data has 2 columns and partial exists,
	does a partial Mantel cross-correlogram,
	calculating partial Mantel r for each distance class separately.

The Iwt statistic used for the multivariate correlograms is not the standard Mantel r. For one variable, using Euclidean distance, this metric converges on the familiar Moran autocorrelation. Like the Moran autocorrelation function, this statistic usually falls between -1 and 1, but is not bounded by those limits. Unlike the Moran function, this correlogram can be used for multivariate data, and can be extended to partial tests.

The comparisons in \code{vignette("dissimilarity", package="ecodist")} may help.
}
\value{
  Returns a object of class mgram, which is a list containing two objects:
  mgram is a matrix with one row for each distance class and 4 columns:
  \item{lag }{midpoint of the distance class.}
  \item{ngroup }{number of distances in that class.}
  \item{piecer or Iwt }{Mantel r value or appropriate statistic  (see Details).}
  \item{pval }{two-sided p-value.}
  resids is a vector of the residuals (if calculated) and can be accessed with the \code{residuals()} method.
}
\author{ Sarah Goslee }

\seealso{ \code{\link{mgram}}, \code{\link{mantel}}, \code{\link{residuals.mgram}}, \code{\link{plot.mgram}} }
\examples{

data(bump)

par(mfrow=c(1, 2))
image(bump, col=gray(seq(0, 1, length=5)))

z <- as.vector(bump)
x <- rep(1:25, times=25)
y <- rep(1:25, each=25)

X <- col(bump)
Y <- row(bump)
# calculate dissimilarities for data and space
geo.dist <- dist(cbind(as.vector(X), as.vector(Y)))
value.dist <- dist(as.vector(bump))

### pmgram() is time-consuming, so this was generated
### in advance and saved.
### set.seed(1234)
### bump.pmgram <- pmgram(value.dist, geo.dist, nperm=10000)

data(bump.pmgram)
plot(bump.pmgram)

#### Partial pmgram example

# generate a simple surface
# with complex nonlinear spatial pattern

x <- matrix(1:25, nrow=25, ncol=25, byrow=FALSE)
y <- matrix(1:25, nrow=25, ncol=25, byrow=TRUE)

# create z1 and z2 as functions of x, y
# and scale them to [0, 1]
z1 <- x + 3*y
z2 <- y - cos(x)

z1 <- (z1 - min(z1)) / (max(z1) - min(z1))
z2 <- (z2 - min(z2)) / (max(z2) - min(z2))

z12 <- (z1 + z2*2)/3

# look at patterns

layout(matrix(c(
1, 1, 2, 2,
1, 1, 2, 2,
3, 3, 4, 4, 
3, 3, 5, 5), nrow=4, byrow=TRUE))


image(z1, col=gray(seq(0, 1, length=20)), zlim=c(0,1))
image(z2, col=gray(seq(0, 1, length=20)), zlim=c(0,1))
image(z12, col=gray(seq(0, 1, length=20)), zlim=c(0,1))

# analyze the pattern of z across space
z1 <- as.vector(z1)
z2 <- as.vector(z2)
z12 <- as.vector(z12)
z1.d <- dist(z1)
z2.d <- dist(z2)
z12.d <- dist(z12)

space <- cbind(as.vector(x), as.vector(y))
space.d <- dist(space)

# take partial correlogram without effects of z1
### pmgram() is time-consuming, so this was generated
### in advance and saved.
### set.seed(1234)
### z.no <- pmgram(z12.d, space.d, nperm=1000, resids=FALSE)
### save(z.no, file="ecodist/data/z.no.rda")
data(z.no)
plot(z.no)


# take partial correlogram of z12 given z1
### pmgram() is time-consuming, so this was generated
### in advance and saved.
### set.seed(1234)
### z.z1 <- pmgram(z12.d, space.d, z2.d, nperm=1000, resids=FALSE)
### save(z.z1, file="ecodist/data/z.z1.rda")
data(z.z1)
plot(z.z1)

}
\keyword{ multivariate }

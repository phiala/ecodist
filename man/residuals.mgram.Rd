\name{residuals.mgram}
\alias{residuals.mgram}
\title{ Residuals of a Mantel correlogram }
\description{
  Extracts residuals from an S3 object of class \code{mgram} (only relevant for objects created by \code{pmgram{}}).
}
\usage{
\method{residuals}{mgram}(object, ...)
}

\arguments{
  \item{object}{ an object of class \code{mgram} }
  \item{\dots}{ additional arguments }
}

\value{
  vector of residuals.
}
\author{ Sarah Goslee }

\seealso{ \code{\link{pmgram}}, \code{\link{mgram}}  }
\examples{

#### Partial pmgram example

# generate a simple surface
# with complex nonlinear spatial pattern

x <- matrix(1:10, nrow=10, ncol=10, byrow=FALSE)
y <- matrix(1:10, nrow=10, ncol=10, byrow=TRUE)

# create z1 and z2 as functions of x, y
# and scale them to [0, 1]
z1 <- x + 3*y
z2 <- y - cos(x)

z1 <- (z1 - min(z1)) / (max(z1) - min(z1))
z2 <- (z2 - min(z2)) / (max(z2) - min(z2))

z12 <- (z1 + z2*2)/3

# analyze the pattern of z across space
z1 <- as.vector(z1)
z2 <- as.vector(z2)
z12 <- as.vector(z12)
z1.d <- dist(z1)
z2.d <- dist(z2)
z12.d <- dist(z12)

space <- cbind(as.vector(x), as.vector(y))
space.d <- dist(space)

# take partial correlogram of z12 given z1
z.z1 <- pmgram(z12.d, space.d, z2.d, nperm=0, resids=TRUE)
summary(residuals(z.z1))

}
\keyword{ multivariate }

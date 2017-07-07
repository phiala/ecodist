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
\author{ Sarah Goslee, Sarah.Goslee@ars.usda.gov }

\seealso{ \code{\link{mgram}} }
\examples{

#### Partial pmgram example

# generate a simple surface
x <- matrix(1:25, nrow=25, ncol=25, byrow=FALSE)
y <- matrix(1:25, nrow=25, ncol=25, byrow=TRUE)

# create z1 and z2 as functions of x, y
# and scale them to [0, 1]
z1 <- x + 3*y
z2 <- y - cos(x)

z1 <- (z1 - min(z1)) / (max(z1) - min(z1))
z2 <- (z2 - min(z2)) / (max(z2) - min(z2))

z12 <- (z1 + z2)/2

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
### pgram() is time-consuming, so this was generated
### in advance and saved.
### set.seed(1234)
### z.no <- pgram(z12.d, space.d, nperm=10000, resids=TRUE)
### save(z.no, file="ecodist/data/z.no.rda")
data(z.no)
plot(z.no)
summary(residuals(z.no))


# take partial correlogram of z12 given z1
### pgram() is time-consuming, so this was generated
### in advance and saved.
### set.seed(1234)
### z.z1 <- pmgram(z12.d, space.d, z2.d, nperm=10000, resids=TRUE)
### save(z.z1, file="ecodist/data/z.z1.rda")
data(z.z1)
plot(z.z1)
summary(residuals(z.z1))

}
\keyword{ multivariate }

\name{plot.mgram}
\alias{plot.mgram}
\title{ Plot a Mantel correlogram }
\description{
  Plot a Mantel correlogram from an object of S3 class \code{mgram}, using solid symbols for significant values.
}
\usage{
\method{plot}{mgram}(x, pval = 0.05, xlab = "Distance", ylab = "Mantel r", \dots)
}

\arguments{
  \item{x}{ an object of class \code{mgram} }
  \item{pval}{ cut-off level for statistical significance. }
  \item{xlab}{ x-axis label. }
  \item{ylab}{ y-axis label. }
  \item{\dots}{ optional, any additional graphics parameters. }
}

\value{
  draws a plot (graphics device must be active).
}
\author{ Sarah Goslee }

\seealso{ \code{\link{mgram}} }
\examples{

# generate a simple surface
x <- matrix(1:10, nrow=10, ncol=10, byrow=FALSE)
y <- matrix(1:10, nrow=10, ncol=10, byrow=TRUE)
z <- x + 3*y
image(z)

# analyze the pattern of z across space
space <- cbind(as.vector(x), as.vector(y))
z <- as.vector(z)
space.d <- distance(space, "eucl")
z.d <- distance(z, "eucl")
z.mgram <- mgram(z.d, space.d, nperm=0)
plot(z.mgram)

#

data(graze)
space.d <- dist(graze$sitelocation)
forest.d <- dist(graze$forestpct)

grasses <- graze[, colnames(graze) \%in\% c("DAGL", "LOAR10", "LOPE", "POPR")]
legumes <- graze[, colnames(graze) \%in\% c("LOCO6", "TRPR2", "TRRE3")]

grasses.bc <- bcdist(grasses)
legumes.bc <- bcdist(legumes)

# Does the relationship of composition with distance vary for
# grasses and legumes?
par(mfrow=c(2, 1))
plot(mgram(grasses.bc, space.d, nclass=8))
plot(mgram(legumes.bc, space.d, nclass=8))

}
\keyword{ multivariate }
\keyword{ hplot }

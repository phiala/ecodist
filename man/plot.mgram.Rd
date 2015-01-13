\name{plot.mgram}
\alias{plot.mgram}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Plot a Mantel correlogram }
\description{
  Plot a Mantel correlogram from an object of S3 class \code{mgram}, using solid symbols for significant values.
}
\usage{
\method{plot}{mgram}(x, pval = 0.05, xlab = "Distance", ylab = "Mantel r", ...)
}
%- maybe also 'usage' for other objects documented here.
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
\author{ Sarah Goslee, Sarah.Goslee@ars.usda.gov }

\seealso{ \code{\link{mgram}} }
\examples{
\dontrun{
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
}
}
\keyword{ multivariate }
\keyword{ hplot }

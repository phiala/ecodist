\name{residuals.mgram}
\alias{residuals.mgram}
\title{ Residuals of a Mantel correlogram }
\description{
  Extracts residuals from an S3 object of class \code{mgram} (only relevant for objects created by \code{pmgram{}}).
}
\usage{
\method{residuals}{mgram}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
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
\dontrun{
z.pmgram <- pmgram(z.d, space.d, resids=TRUE)
residuals(z.pmgram)
}

# A full example is available in the Mantel correlogram 
# section of the main help file for \link{ecodist}.
}
\keyword{ multivariate }


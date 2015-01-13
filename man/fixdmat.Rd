\name{fixdmat}
\alias{fixdmat}
\title{ Distance matrix conversion }
\description{
  Converts a row-order lower-triangular distance matrix to a full symmetric matrix.
}
\usage{
fixdmat(v)
}
\arguments{
  \item{v}{ lower-triangular distance matrix in row order. }
}
\details{
R and S-Plus distance functions such as dist and bcdist return a lower triangular distance matrix in column order. Some other programs, such as EUCLID, return the lower triangular matrix in row order. To use this matrix in R/S-Plus functions, it must be converted from row order to column order.
}
\value{
  Full symmetric distance matrix.
}
\author{ Sarah Goslee, Sarah.Goslee@ars.usda.gov }

\seealso{ \code{\link{lower}},  \code{\link{full}},  \code{\link{dist}}, \code{\link{bcdist}} }

\keyword{ manip }% at least one, from doc/KEYWORDS

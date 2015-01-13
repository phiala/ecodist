\name{full}
\alias{full}
\title{ Full symmetric matrix }
\description{
  Converts a column order distance matrix to a full symmetric matrix.
}
\usage{
full(v)
}
\arguments{
  \item{v}{ lower-triangular column order distance matrix. }
}
\details{
  Converts lower-triangular distance matrix as written by R/S-Plus functions 
into a symmetric matrix.
  Note that lower() used on a 1x1 matrix will return the single element, which may not be the correct behavior in all cases, while full() used on a single element will return a 2x2 matrix.
}
\value{
 a full symmetric matrix.
}
\author{ Sarah Goslee, Sarah.Goslee@ars.usda.gov }

\seealso{ \code{\link{lower}}  }

\keyword{ manip }% at least one, from doc/KEYWORDS


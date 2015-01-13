\name{lower}
\alias{lower}
\title{ Lower-triangular matrix }
\description{
  Converts a symmetric distance matrix to a column order lower triangular matrix.
}
\usage{
lower(m)
}
\arguments{
  \item{m}{ a symmetric distance matrix. }
}
\details{
  Converts a symmetric matrix, for example a dissimilarity matrix, into a lower-
triangular matrix. This may be useful to format the input for certain clustering and ordination functions.
  Note that lower() used on a 1x1 matrix will return the single element, which may not be the correct behavior in all cases, while full() used on a single element will return a 2x2 matrix.
}
\value{
 a column order lower triangular matrix.
 }
\author{ Sarah Goslee, Sarah.Goslee@ars.usda.gov }

\seealso{ \code{\link{full}}, \code{\link{fixdmat}}, }

\keyword{ manip }% at least one, from doc/KEYWORDS


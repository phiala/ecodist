\name{full}
\alias{full}
\title{ Full symmetric matrix }
\description{
  Convert a column order distance matrix to a full symmetric matrix.
}
\usage{
  full(v)
}
\arguments{
  \item{v}{ lower-triangular column order distance matrix. }
}
\details{
Converts a column order lower-triangular distance matrix as
written by R functions into a symmetric matrix.
Note that \code{lower()} used on a 1x1 matrix will return the single
element, which may not be the correct behavior in all cases,
while \code{full()} used on a single element will return a 2x2 matrix.
}
\value{
  full symmetric matrix.
}
\author{ Sarah Goslee }

\seealso{ \code{\link{lower}}, \code{\link{fixdmat}} }

\examples{
# Given a vector:
x.vec <- seq_len(6)
x.vec

# Make an R-style column order symmetric matrix
full(x.vec)

# Extract the lower triangle from a symmetric matrix
# in column order
lower(full(x.vec))

# Convert to or from a row order symmetric matrix
fixdmat(x.vec)
lower(fixdmat(x.vec))

fixdmat(c(1, 2, 4, 3, 5, 6))
}

\keyword{ manip }
\keyword{ multivariate }

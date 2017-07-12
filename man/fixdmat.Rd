\name{fixdmat}
\alias{fixdmat}
\title{ Distance matrix conversion }
\description{
  Convert a row-order lower-triangular distance matrix to a full symmetric matrix.
}
\usage{
  fixdmat(v)
}
\arguments{
  \item{v}{ lower-triangular distance matrix in row order. }
}
\details{
R distance functions such as dist and bcdist return a lower
triangular distance matrix in column order. Some other programs
return the lower triangular matrix in row order. To use this
matrix in R functions, it must be converted from row order
to column order.
}
\value{
  full symmetric distance matrix.
}
\author{ Sarah Goslee }
\seealso{ \code{\link{lower}},  \code{\link{full}} }
\examples{
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

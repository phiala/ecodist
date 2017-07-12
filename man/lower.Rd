\name{lower}
\alias{lower}
\title{ Lower-triangular matrix }
\description{
  Convert a symmetric distance matrix to a column order lower triangular matrix.
}
\usage{
  lower(m)
}
\arguments{
  \item{m}{ a symmetric distance matrix. }
}
\details{
Converts a symmetric matrix, for example a dissimilarity matrix,
into a column order lower-triangular matrix. This may be useful
to format the input for certain clustering and ordination functions.
Note that \code{lower()} used on a 1x1 matrix will return the single
element, which may not be the correct behavior in all cases,
while \code{full()} used on a single element will return a 2x2 matrix.
}
\value{
  column order lower triangular matrix.
}
\author{ Sarah Goslee }

\seealso{ \code{\link{full}}, \code{\link{fixdmat}} }

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

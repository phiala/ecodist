\name{mstdist}
\alias{mstdist}
\title{
  Minimum spanning tree extension of dissimilarities
}
\description{
	Uses minimum spanning tree path length to estimate the distance between samples with pairwise distances greater than maxv.
}
\usage{
mstdist(v, maxv = 1)
}
\arguments{
  \item{v}{
	lower-triangular distance vector, possibly as produced by dist() or distance().
}
  \item{maxv}{
	cutoff for distances: values greater or equal to this will be estimated from the minimum spanning tree.
}
}
\details{
 	Pairwise samples with no species will have distances greater than a cutoff. A minimum spanning tree connecting these samples by way of intermediate samples with some species in common can be used to interpolate distances by adding up the path length connecting those samples.
}
\value{
	Returns a lower-triangular distance matrix.
}


\author{ Sarah Goslee }

\seealso{
	\code{\link{dist}}, \code{\link{distance}}
}
\examples{

	# samples 1 and 2, and 3 and 4, have no species in common
	x <- matrix(c(	1, 0, 1, 0,
			0, 1, 0, 1,
			1, 0, 0, 0,
			0, 1, 1, 1,
			1, 1, 1, 0,
			1, 0, 1, 1,
			0, 0, 1, 1), ncol = 4, byrow = TRUE)

	# the maximum Jaccard distance is 1
	# regardless of how different the samples are
	x.jd <- dist(x, "binary")

	# using the MST, estimate the true distance between those pairs
	mstdist(x.jd)

}
\keyword{ multivariate }


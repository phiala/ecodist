\name{xdistance}
\alias{xdistance}
\title{
Cross-distance between two datasets.
}
\description{
Pairwise dissimilarity calculation between rows of one dataset and rows of another, for instance across different sampling periods for the same set of sites.
}
\usage{
xdistance(x, y, method = "euclidean")
}

\arguments{
  \item{x}{A site by species or other matrix or data frame.}
  \item{y}{A a second site by species dataset, which must have at least the same columns.}
  \item{method}{This function calls \code{\link{distance}} to do the calculations, and will accept any method used there, currently: "euclidean", "bray-curtis", "manhattan", "mahalanobis" (squared Mahalanobis distance), "jaccard", "difference", "sorensen", "gower", "modgower10" (modified Gower, base 10), "modgower2" (modified Gower, base 2). Partial matching will work for selecting a method.}
}
\details{
This function will calculate rowwise dissimilarities between any pair of matrices or data frames with the same number of columns. If they also have the same number of rows, the resulting square non-symmetric matrices can be used with \code{xmantel} or \code{xmgram}.
}
\value{
A non-symmetric and possibly not square matrix of dissimilarities, where \code{result <- xdistance(x, y)} produces a matrix with \code{result[a, b]} containing the dissimilarity between \code{x[a, ]} and \code{y[b, ]}.
}
\author{ Sarah Goslee <Sarah.Goslee@ars.usda.gov> }

\seealso{ \code{\link{distance}}, \code{\link{xmantel}}, \code{\link{xmgram}}}
\examples{
# Need to develop a cross-dissimilarity example
}

\keyword{ multivariate }

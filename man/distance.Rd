\name{distance}
\alias{distance}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Calculate dissimilarity/distance metrics }
\description{
  This function calculates a variety of dissimilarity or distance metrics. Although it duplicates the functionality of dist() and bcdist(), it is written in such a way that new metrics can easily be added.
  distance() was written for extensibility and understandability, and is not an efficient choice for use with large matrices.
}
\usage{
distance(x, method = "euclidean", sprange=NULL, spweight=NULL, icov)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ matrix or data frame with rows as samples and columns as variables (such as 
species). Distances will be calculated for each pair of rows. }
  \item{method}{ Currently 7 dissimilarity metrics can be calculated: "euclidean", "bray-curtis", "manhattan", "mahalanobis" (squared Mahalanobis distance), "jaccard", "difference", "sorensen", "gower", "modgower10" (modified Gower, base 10), "modgower2" (modified Gower, base 2). Partial matching will work for selecting a method.}
  \item{sprange}{Gower dissimilarities offer the option of dividing by the species range. If sprange=NULL no range is used. If sprange is a vector of length nrow(x) it is used for standardizing the dissimilarities.}
  \item{spweight}{Euclidean, Manhattan, and Gower dissimilarities allow weighting. If spweight=NULL, no weighting is used. If spweight="absence", then W=0 if both species are absent and 1 otherwise, thus deleting joint absences.}
  \item{icov}{Optional; only used if method="mahalanobis" since Mahalanobis distance requires calculating the variance-covariance matrix for the entire dataset. Providing icov directly makes it possible to calculate distances for a subset of the full dataset.}
}

\value{
  Returns a lower-triangular distance matrix as an object of class "dist".
  
}
\author{ Sarah Goslee <Sarah.Goslee@ars.usda.gov> }

\seealso{ \code{\link{dist}} }
\examples{
data(iris)
iris.md <- distance(iris[,1:4], "mahal")

}
\keyword{ multivariate }

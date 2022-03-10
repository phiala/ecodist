\name{distance}
\alias{distance}
\title{ Calculate dissimilarity/distance metrics }
\description{
  This function calculates a variety of dissimilarity or distance metrics. Although it duplicates the functionality of dist() and bcdist(), it is written in such a way that new metrics can easily be added.
  distance() was written for extensibility and understandability, and is not necessarily an efficient choice for use with large matrices.
}
\usage{
distance(x, method = "euclidean", sprange=NULL, spweight=NULL, icov, inverted = FALSE)
}

\arguments{
  \item{x}{ matrix or data frame with rows as samples and columns as variables (such as
species). Distances will be calculated for each pair of rows. }
  \item{method}{ Currently 7 dissimilarity metrics can be calculated: "euclidean", "bray-curtis", "manhattan", "mahalanobis" (squared Mahalanobis distance), "jaccard", "difference", "sorensen", "gower", "modgower10" (modified Gower, base 10), "modgower2" (modified Gower, base 2). Partial matching will work for selecting a method.}
  \item{sprange}{Gower dissimilarities offer the option of dividing by the species range. If sprange=NULL no range is used. If sprange is a vector of length nrow(x) it is used for standardizing the dissimilarities.}
  \item{spweight}{Euclidean, Manhattan, and Gower dissimilarities allow weighting. If spweight=NULL, no weighting is used. If spweight="absence", then W=0 if both species are absent and 1 otherwise, thus deleting joint absences.}
  \item{icov}{Optional covariance matrix; only used if method="mahalanobis" since Mahalanobis distance requires calculating the variance-covariance matrix for the entire dataset. Providing icov directly makes it possible to calculate distances for a subset of the full dataset.}
  \item{inverted}{If TRUE, the optional covariance matrix for method="mahalanobis" is not inverted before solving. Providing an inverted matrix may speed up calculations.}
}

\value{
  Returns a lower-triangular distance matrix as an object of class "dist".

}
\author{ Sarah Goslee }

\seealso{ \code{\link{dist}},  \code{\link{bcdist}} }
\examples{
data(iris)
iris.bc <- distance(iris[, 1:4], "bray-curtis")

# The effect of specifying icov:

# calculate Mahalanobis distance for the full iris dataset
iris.md <- full(distance(iris[, 1:4], "mahal"))
iris.md[1, 2] # Mahalanobis distance between samples 1 and 2 

# calculate Mahalanobis for just one species
setosa.md <- full(distance(iris[iris$Species == "setosa", 1:4], "mahal"))
setosa.md[1, 2] # Mahalanobis distance between samples 1 and 2 

# use the covariance matrix for the full dataset to scale for one species
setosa.scaled.md <- full(distance(iris[iris$Species == "setosa", 1:4],
  "mahal", icov=var(iris[,1:4])))
setosa.scaled.md[1, 2] # Mahalanobis distance between samples 1 and 2 

}
\keyword{ multivariate }

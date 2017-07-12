\name{mgroup}
\alias{mgroup}
\title{ Mantel test for groups }
\description{
  Simple and partial Mantel tests, with options for ranked data, 
  permutation tests, and bootstrapped confidence limits.
}
\usage{
mgroup(edist, groups, nperm=1000)
}

\arguments{
  \item{edist}{ a dist object or lower triangular distance vector. }
  \item{groups}{ a vector of group memberships (numeric, character, or 
  factor), or a matrix or data frame with columns describing multiple 
  sets of groups. }
  \item{nperm}{ number of permutations to use. If set to 0, the 
  permutation test will be omitted. }
}
\details{
 	\code{mgroup} returns the Mantel correlations for group contrast 
	matrices computed from cluster groups across a range of clustering 
	levels.
}
\value{
  \item{nclust }{Number of groups tested.}
  \item{mantelr }{Mantel coefficient.}
  \item{pval1 }{one-tailed p-value (null hypothesis: r <= 0).}
}
\references{ Legendre, P. and M. Fortin. 1989. Spatial pattern and ecological 
analysis. Vegetatio 80:107-138. }

\author{ Sarah Goslee }

\seealso{ \code{\link{mantel}} }
\examples{

# Using a model matrix to test group membership

data(iris)
iris.d <- dist(iris[,1:4])
mgroup(iris.d, iris[,5])

# clustering-based example

data(graze)
graze.d <- dist(graze[, -c(1:2)])
graze.hclust <- hclust(graze.d)

clust.groups <- data.frame(
	k2 = cutree(graze.hclust, k = 2),
	k4 = cutree(graze.hclust, k = 4),
	k6 = cutree(graze.hclust, k = 6),
	k8 = cutree(graze.hclust, k = 8))

mgroup(graze.d, clust.groups, nperm=1000)

}

\keyword{ multivariate }

\name{nmds}
\alias{nmds}
\title{ Non-metric multidimensional scaling }
\description{
  Non-metric multidimensional scaling.
}
\usage{
nmds(dmat, mindim = 1, maxdim = 2, nits = 10, iconf = 0, epsilon = 1e-12, 
    maxit = 500, trace = FALSE, stresscalc="default")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dmat}{ lower-triangular dissimilarity matrix. }
  \item{mindim}{ optional, minimum number of dimensions to use. }
  \item{maxdim}{ optional, maximum number of dimensions to use. }
  \item{nits}{ optional, number of separate ordinations to use. }
  \item{iconf}{ optional, initial configuration. If not specified, then a random configuration is used. } 
  \item{epsilon}{ optional, acceptable difference in stress. }
  \item{maxit}{ optional, maximum number of iterations. }
  \item{trace}{ if TRUE, will write progress indicator to the screen. }
  \item{stresscalc}{ if "default" uses modified stress calculation; if "kruskal" uses unmodified Kruskal stress calculation. }
}
\details{
  The goal of NMDS is to find a configuration in a given number of dimensions which preserves rank-order dissimilarities as closely as possible. The number of dimensions must be specified in advance. Because NMDS is prone to finding local minima, several random starts must be used. Stress is used as the measure of goodness of fit. A lower stress indicates a better match between dissimilarity and ordination.
}
\value{
  \item{conf }{list of configurations.}
  \item{stress }{list of final stress values.}
  \item{r2 }{total variance explained by each configuration.}

  The first results are for the lowest number of dimensions (total number is (mindim - maxdim + 1) * nits).

}
\references{ Kruskal, J.B. 1964. Multidimensional scaling by optimizing goodness of fit to a 
nonmetric hypothesis. Psychometrika 29:1-27.

Minchin, P.R. 1987. An evaluation of the relative robustness of techniques for ecological ordination. Vegetatio 96:89-108. }
\author{ Sarah Goslee, Sarah.Goslee@ars.usda.gov }

\seealso{ \code{\link{nmds.min}}, \code{\link{vf}} }
\examples{

\dontrun{
# Example of multivariate analysis using built-in iris dataset
data(iris)
iris <- iris[seq(1, 150, by=3),]
iris.md <- distance(iris[,1:4], "mahal")

# Minimum-stress 2-dimensional nonmetric multidimensional scaling configuration
# Uses small number of separate ordinations (5) to increase speed of example.
# Use more for final analysis.
iris.nmds <- nmds(iris.md, mindim=2, maxdim=2, nits=3)
iris.nmin <- nmds.min(iris.nmds)

# Plot NMDS result with symbols denoting species
plot(iris.nmin, pch=as.numeric(iris[,5]))

# Fit vectors for the main variables to the NMDS configuration
iris.vf <- vf(iris.nmin, iris[,1:4], nperm=10)
plot(iris.vf, col="blue")
}
}
\keyword{ multivariate }


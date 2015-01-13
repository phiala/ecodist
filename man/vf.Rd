\name{vf}
\alias{vf}
\title{ Vector fitting }
\description{
  Fits environmental variables to an ordination configuration.
}
\usage{
vf(ord, vars, nperm = 100)
}
\arguments{
  \item{ord}{ matrix containing a 2-dimensional ordination result with axes as columns. }
  \item{vars}{ matrix with ancillary variables as columns. }
  \item{nperm}{ number of permutation for the significance test. If nperm = 0, the test will be omitted. }
}
\details{
 Vector fitting finds the maximum correlation of the individual variables with a configuration of samples in ordination space. 
}
\value{
  an object of class vf containing matrix with the first 2 columns containing the scores for every variable in each of the 2 dimensions of the ordination space. r is the maximum correlation of the variable with the ordination space, and pval is the result of the permutation test.
  }
\references{ Jongman, R.H.G., C.J.F. ter Braak and O.F.R. van Tongeren. 1995. Data analysis in community and landscape ecology. Cambridge University Press, New York. }
\author{ Sarah Goslee, Sarah.Goslee@ars.usda.gov }

\seealso{ \code{\link{plot.vf}} }
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

\name{vf}
\alias{vf}
\title{ Vector fitting }
\description{
  Fits ancillary variables to an ordination configuration.
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

# Example of multivariate analysis using built-in iris dataset
data(iris)
iris.d <- dist(iris[,1:4])

### nmds() is timeconsuming, so this was generated
### in advance and saved.
### set.seed(1234)
### iris.nmds <- nmds(iris.d, nits=20, mindim=1, maxdim=4)
### save(iris.nmds, file="ecodist/data/iris.nmds.rda")
data(iris.nmds)

# examine fit by number of dimensions
plot(iris.nmds)

# choose the best two-dimensional solution to work with
iris.nmin <- min(iris.nmds, dims=2)

# fit the data to the ordination as vectors
### vf() is timeconsuming, so this was generated
### in advance and saved.
### set.seed(1234)
### iris.vf <- vf(iris.nmin, iris[,1:4], nperm=1000)
### save(iris.vf, file="ecodist/data/iris.vf.rda")
data(iris.vf)
plot(iris.nmin, col=as.numeric(iris$Species), pch=as.numeric(iris$Species), main="NMDS")
plot(iris.vf)

}
\keyword{ multivariate }

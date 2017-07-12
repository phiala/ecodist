\name{nmds}
\alias{nmds}
\title{ Non-metric multidimensional scaling }
\description{
  Non-metric multidimensional scaling.
}
\usage{
nmds(dmat, mindim = 1, maxdim = 2, nits = 10, iconf, epsilon = 1e-12,
    maxit = 500, trace = FALSE)
}

\arguments{
  \item{dmat}{ lower-triangular dissimilarity matrix. }
  \item{mindim}{ optional, minimum number of dimensions to use. }
  \item{maxdim}{ optional, maximum number of dimensions to use. }
  \item{nits}{ optional, number of separate ordinations to use. }
  \item{iconf}{ optional, initial configuration. If not specified, then a random configuration is used. }
  \item{epsilon}{ optional, acceptable difference in stress. }
  \item{maxit}{ optional, maximum number of iterations. }
  \item{trace}{ if TRUE, will write progress indicator to the screen. }
}
\details{
  The goal of NMDS is to find a configuration in a given number of dimensions which preserves rank-order dissimilarities as closely as possible. The number of dimensions must be specified in advance. Because NMDS is prone to finding local minima, several random starts must be used. Stress is used as the measure of goodness of fit. A lower stress indicates a better match between dissimilarity and ordination.
  As of ecodist 1.9, the stress calculation used is the same as in \code{MASS:isoMDS}. In previous versions it was monotonically related, so the same configurations were produced, but the absolute value was different.
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

\seealso{ \code{\link{plot.nmds}}, \code{\link{min.nmds}}, \code{\link{vf}}, \code{\link{addord}} }
\examples{
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

# rotate the configuration to maximize variance
iris.rot <- princomp(iris.nmin)$scores

# rotation preserves distance apart in ordination space
cor(dist(iris.nmin), dist(iris.rot))

# fit the data to the ordination as vectors
### vf() is timeconsuming, so this was generated
### in advance and saved.
### set.seed(1234)
### iris.vf <- vf(iris.nmin, iris[,1:4], nperm=1000)
### save(iris.vf, file="ecodist/data/iris.vf.rda")
data(iris.vf)

# repeat for the rotated ordination
### vf() is timeconsuming, so this was generated
### in advance and saved.
### set.seed(1234)
### iris.vfrot <- vf(iris.rot, iris[,1:4], nperm=1000)
### save(iris.vfrot, file="ecodist/data/iris.vfrot.rda")
data(iris.vfrot)

par(mfrow=c(1,2))
plot(iris.nmin, col=as.numeric(iris$Species), pch=as.numeric(iris$Species), main="NMDS")
plot(iris.vf)
plot(iris.rot, col=as.numeric(iris$Species), pch=as.numeric(iris$Species),
    main="Rotated NMDS")
plot(iris.vfrot)


# generate new data points to add to the ordination
# this might be new samples, or a second dataset

iris.new <- structure(list(Sepal.Length = c(4.6, 4.9, 5.4, 5.2, 6, 6.5, 6, 
6.8, 7.3), Sepal.Width = c(3.2, 3.5, 3.6, 2.3, 2.8, 3, 2.7, 3.1, 
3.2), Petal.Length = c(1.2, 1.5, 1.5, 3.5, 4.1, 4.2, 4.8, 5, 
5.7), Petal.Width = c(0.26, 0.26, 0.26, 1.2, 1.3, 1.4, 1.8, 2, 
2), Species = structure(c(1L, 1L, 1L, 2L, 2L, 2L, 3L, 3L, 3L), .Label = c("setosa", 
"versicolor", "virginica"), class = "factor")), .Names = c("Sepal.Length", 
"Sepal.Width", "Petal.Length", "Petal.Width", "Species"), class = "data.frame",
row.names = c(NA, -9L))

# provide a dist object containing original and new data
# provide a logical vector indicating which samples were used to
# construct the original configuration

iris.full <- rbind(iris, iris.new)
all.d <- dist(iris.full[,1:4])
is.orig <- c(rep(TRUE, nrow(iris)), rep(FALSE, nrow(iris.new)))

### addord() is timeconsuming, so this was generated
### in advance and saved.
### set.seed(1234)
### iris.fit <- addord(iris.nmin, iris.full[,1:4], all.d, is.orig, maxit=100)
### save(iris.fit, file="ecodist/data/iris.fit.rda")
data(iris.fit)

plot(iris.fit$conf, col=iris.full$Species, pch=c(18, 4)[is.orig + 1],
    xlab="NMDS 1", ylab="NMDS 2")
title("Demo: adding points to an ordination")
legend("bottomleft", c("Training set", "Added point"), pch=c(4, 18))
legend("topright", levels(iris$Species), fill=1:3)
}
\keyword{ multivariate }

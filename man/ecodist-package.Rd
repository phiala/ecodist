\name{ecodist-package}
\alias{ecodist-package}
\alias{ecodist}
\docType{package}
\title{
Dissimilarity-based functions for ecological analysis
}
\description{
Dissimilarity-based analysis functions including
ordination and Mantel test functions with multiple partials, intended for use with
spatial and community data.
}
\details{
\tabular{ll}{
Package: \tab ecodist\cr
Version: \tab 1.2.2\cr
Date: \tab 2011-04-02\cr
Depends: \tab stats\cr
License: \tab GPL version 2 or newer\cr
}

Index:
\preformatted{
bcdist                  Bray-Curtis distance
cor2m                   Generates a correlation table between the
                        variables of 2 matrices
corgen                  Generate correlated data
crosstab                Data formatting
distance                Calculate dissimilarity/distance metrics
fixdmat                 Distance matrix conversion
full                    Full symmetric matrix
lower                   Lower-triangular matrix
mantel                  Mantel test
mgram                   Mantel correlogram
nmds                    Non-metric multidimensional scaling
nmds.min                Find minimum stress configuration
pco                     Principal coordinates analysis
plot.mgram              Plot a Mantel correlogram
plot.vf                 Plot fitted vectors onto an ordination diagram
pmgram                  Partial Mantel correlogram
vf                      Vector fitting
}
}
\author{
Sarah Goslee and Dean Urban

Maintainer: Sarah Goslee <Sarah.Goslee@ars.usda.gov>
}
\references{
Goslee, S.C. and Urban, D.L. 2007. The ecodist package for dissimilarity-based analysis of ecological data. Journal of Statistical Software 22(7):1-19.
}
\keyword{ multivariate }
\examples{

## Distance example using a subset of the built-in iris dataset
data(iris)
iris <- iris[seq(1, 150, by=3),]
iris.md <- distance(iris[,1:4], "mahal")
iris.bc <- bcdist(iris[,1:4])

# compare Mahalanobis and Bray-Curtis dissimilarities
plot(iris.md, iris.bc, xlab="Mahalanobis", ylab="Bray-Curtis")


## NMDS example
# Example of multivariate analysis using built-in iris dataset

# Minimum-stress 2-dimensional nonmetric multidimensional scaling configuration
# Uses only 3 starting configuration to increase speed of example.
# Use more for final analysis to make it more likely that you find the
# global minimum configuration.
iris.nmds <- nmds(iris.md, mindim=2, maxdim=2, nits=3)
iris.nmin <- nmds.min(iris.nmds)

# Plot NMDS result with symbols denoting species
plot(iris.nmin, pch=as.numeric(iris[,5]))

# Fit vectors for the main variables to the NMDS configuration
# Use more permutations for analysis and publication.
iris.vf <- vf(iris.nmin, iris[,1:4], nperm=50)
plot(iris.vf, col="blue")


## PCO example

iris.pco <- pco(iris.md)
# scatterplot of the first two dimensions
plot(iris.pco$vectors[,1], iris.pco$vectors[,2], pch=as.numeric(iris[,5]))


## Mantel test example

# Example of multivariate analysis using built-in iris dataset

# Create a model matrix for testing species differences
iris.model <- distance(as.numeric(iris[,5]), "eucl")
iris.model[iris.model > 0] <- 1

# Test whether samples within the same species are more similar than those not
# Use very low numbers of permutations to increase speed of example.
# Use many permutations for final analyses!
mantel(iris.md ~ iris.model, nperm=50, nboot=0)


## Mantel correlogram example

# generate a simple surface
x <- matrix(1:10, nrow=10, ncol=10, byrow=FALSE)
y <- matrix(1:10, nrow=10, ncol=10, byrow=TRUE)
z <- x + 3*y
image(z)

# analyze the pattern of z across space
space <- cbind(as.vector(x), as.vector(y))
z <- as.vector(z)
space.d <- distance(space, "eucl")
z.d <- distance(z, "eucl")
z.mgram <- mgram(z.d, space.d, nperm=0)
plot(z.mgram)


## Partial Mantel correlogram example
# generate a simple surface
x <- matrix(1:10, nrow=10, ncol=10, byrow=FALSE)
y <- matrix(1:10, nrow=10, ncol=10, byrow=TRUE)
z1 <- x + 3*y
z2 <- 2*x - y

# look at patterns
par(mfrow=c(1,2))
image(z1)
image(z2)

# analyze the pattern of z across space
z1 <- as.vector(z1)
z2 <- as.vector(z2)
z1.d <- distance(z1, "eucl")
z2.d <- distance(z2, "eucl")

space <- cbind(as.vector(x), as.vector(y))
space.d <- distance(space, "eucl")

# take partial correlogram of z2 on the residuals of z1 ~ space.d
z.pmgram <- pmgram(z1.d, space.d, z2.d, nperm=0)
par(mfrow=c(1,1))
plot(z.pmgram, pval=0.1)

}


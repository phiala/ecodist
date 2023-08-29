\name{mgram}
\alias{mgram}
\title{ Mantel correlogram }
\description{
	Calculates simple Mantel correlograms.
}
\usage{
mgram(species.d, space.d, breaks, nclass, stepsize, equiprobable = FALSE, nperm = 1000,
    mrank = FALSE, nboot = 500, pboot = 0.9, cboot = 0.95,
    alternative = "two.sided", trace = FALSE)
}
\arguments{
  \item{species.d}{ lower-triangular dissimilarity matrix. }
  \item{space.d}{ lower-triangular matrix of geographic distances. }
  \item{breaks}{ locations of class breaks. If specified, overrides nclass and stepsize. }
  \item{nclass}{ number of distance classes. If not specified, Sturge's rule will be used
to determine an appropriate number of classes. }
  \item{stepsize}{ width of each distance class. If not specified, nclass and the range of space.d will be used to calculate an appropriate default. }
  \item{equiprobable}{ if TRUE, create nclass classes of equal number of distances; if FALSE, create nclass classes of equal width }  
  \item{nperm}{ number of permutations to use. If set to 0, the permutation test will be omitted. }
  \item{mrank}{ if this is set to FALSE (the default option), Pearson correlations will be used. If
 set to TRUE, the Spearman correlation (correlation ranked distances) will be used. }
  \item{nboot}{ number of iterations to use for the bootstrapped confidence limits. If set to 0,
 the bootstrapping will be omitted. }
  \item{pboot}{ the level at which to resample the data for the bootstrapping procedure. }
  \item{cboot}{ the level of the confidence limits to estimate. }
  \item{alternative}{ default is "two.sided", and returns p-values for H0: rM = 0. The alternative is "one.sided", which returns p-values for H0: rM <= 0. }
  \item{trace}{ if TRUE, returns progress indicators. }
}
\details{
This function calculates Mantel correlograms, and tests the hypothesis that the mean compositional dissimilarity within a distance class differs from the mean of all the other distance classes combined. The Mantel correlogram is essentially a multivariate autocorrelation function. The Mantel r represents the dissimilarity in variable composition  (often species composition) at a particular lag distance, and significance is tested in reference to all distance classes.
}
\value{
	Returns an object of class mgram, which is a list with two elements.
    mgram is a matrix with one row for each distance class and 6 columns:
  \item{lag }{midpoint of the distance class.}
  \item{ngroup }{number of distances in that class.}
  \item{mantelr }{Mantel r value.}
  \item{pval }{p-value for the test chosen.}
  \item{llim }{lower bound of confidence limit for mantelr.}
  \item{ulim }{upper bound of confidence limit for mantelr.}
  resids is NA for objects calculated by mgram().
}
\references{ Legendre, P. and M. Fortin. 1989. Spatial pattern and ecological analysis.
Vegetatio 80:107-138. }
\author{ Sarah Goslee }

\seealso{ \code{\link{mantel}}, \code{\link{plot.mgram}}, \code{\link{pmgram}} }
\examples{

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

#

data(graze)
space.d <- dist(graze$sitelocation)
forest.d <- dist(graze$forestpct)

grasses <- graze[, colnames(graze) \%in\% c("DAGL", "LOAR10", "LOPE", "POPR")]
legumes <- graze[, colnames(graze) \%in\% c("LOCO6", "TRPR2", "TRRE3")]

grasses.bc <- bcdist(grasses)
legumes.bc <- bcdist(legumes)

# Does the relationship of composition with distance vary for
# grasses and legumes?
par(mfrow=c(2, 1))
plot(mgram(grasses.bc, space.d, nclass=8))
plot(mgram(legumes.bc, space.d, nclass=8))

}
\keyword{ multivariate }

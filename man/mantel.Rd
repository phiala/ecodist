\name{mantel}
\alias{mantel}
\title{ Mantel test }
\description{
  Simple and partial Mantel tests, with options for ranked data, permutation tests, and bootstrapped confidence limits.
}
\usage{
mantel(formula = formula(data), data, nperm = 1000,
    mrank = FALSE, nboot = 500, pboot = 0.9, cboot = 0.95)
}

\arguments{
  \item{formula}{ formula describing the test to be conducted. For this test, y ~ x will perform a simple Mantel test, while y ~ x + z1 + z2 + z3 will do a partial Mantel test of the relationship between x and y given z1, z2, z3. All variables can be either a distance matrix of class dist or vectors of dissimilarities. }
  \item{data}{ an optional dataframe containing the variables in the model as columns of dissimilarities. By default the variables are taken from the current environment. }
  \item{nperm}{ number of permutations to use. If set to 0, the permutation test will be omitted. }
  \item{mrank}{ if this is set to FALSE (the default option), Pearson correlations will be used. If
 set to TRUE, the Spearman correlation (correlation ranked distances) will be used. }
  \item{nboot}{ number of iterations to use for the bootstrapped confidence limits. If set to 0,
 the bootstrapping will be omitted. }
  \item{pboot}{ the level at which to resample the data for the bootstrapping procedure. }
  \item{cboot}{ the level of the confidence limits to estimate. }
}
\details{
  If only one independent variable is given, the simple Mantel r (r12) is calculated. If 
  more than one independent variable is given, the partial Mantel r  (ryx|x1 ...) is 
  calculated by permuting one of the original dissimilarity matrices.
The bootstrapping is actually resampling without replacement, because duplication of 
samples is not useful in a dissimilarity context (the dissimilarity of a sample with 
itself is zero). Resampling within dissimilarity values is inappropriate, just as for 
permutation.
}
\value{
  \item{mantelr }{Mantel coefficient.}
  \item{pval1 }{one-tailed p-value (null hypothesis: r <= 0).}
  \item{pval2 }{one-tailed p-value (null hypothesis: r >= 0).}
  \item{pval3 }{two-tailed p-value (null hypothesis: r = 0).}
  \item{llim }{lower confidence limit.}
  \item{ulim }{upper confidence limit.}
}
\references{ Mantel, N. 1967. The detection of disease clustering and a generalized
regression approach. Cancer Research 27:209-220.

Smouse, P.E., J.C. Long and R.R. Sokal. 1986. Multiple regression and correlation
extensions of the Mantel test of matrix correspondence. Systematic Zoology 35:62
7-632.

Goslee, S.C. and Urban, D.L. 2007. The ecodist package for dissimilarity-based analysis
of ecological data. Journal of Statistical Software 22(7):1-19.

Goslee, S.C. 2010. Correlation analysis of dissimilarity matrices. Plant
Ecology 206(2):279-286.

}
\author{ Sarah Goslee }

\seealso{ \code{\link{mgram}}, \code{\link{mgroup}} }
\examples{

data(graze)

grasses <- graze[, colnames(graze) \%in\% c("DAGL", "LOAR10", "LOPE", "POPR")]
legumes <- graze[, colnames(graze) \%in\% c("LOCO6", "TRPR2", "TRRE3")]

grasses.bc <- bcdist(grasses)
legumes.bc <- bcdist(legumes)

space.d <- dist(graze$sitelocation)
forest.d <- dist(graze$forestpct)

# Mantel test: is the difference in forest cover between sites
# related to the difference in grass composition between sites?
mantel(grasses.bc ~ forest.d)

# Mantel test: is the geographic distance between sites
# related to the difference in grass composition between sites?
mantel(grasses.bc ~ space.d)

# Partial Mantel test: is the difference in forest cover between sites
# related to the difference in grass composition once the
# linear effects of geographic distance are removed?
mantel(grasses.bc ~ forest.d + space.d)


# Mantel test: is the difference in forest cover between sites
# related to the difference in legume composition between sites?
mantel(legumes.bc ~ forest.d)

# Mantel test: is the geographic distance between sites
# related to the difference in legume composition between sites?
mantel(legumes.bc ~ space.d)

# Partial Mantel test: is the difference in forest cover between sites
# related to the difference in legume composition once the
# linear effects of geographic distance are removed?
mantel(legumes.bc ~ forest.d + space.d)


# Is there nonlinear pattern in the relationship with geographic distance?
par(mfrow=c(2, 1))
plot(mgram(grasses.bc, space.d, nclass=8))
plot(mgram(legumes.bc, space.d, nclass=8))
}

\keyword{ multivariate }

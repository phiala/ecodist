\name{mantel}
\alias{mantel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Mantel test }
\description{
  Simple and partial Mantel tests, with options for ranked data, permutation tests, and bootstrapped confidence limits.
}
\usage{
mantel(formula = formula(data), data = sys.parent(), nperm = 1000,
    mrank = FALSE, nboot = 500, pboot = 0.9, cboot = 0.95)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{ formula in R/S-Plus format describing the test to be conducted. For this test, y ~ x will perform a simple Mantel test, while y ~ x + z1 + z2 + z3 will do a partial Mantel test of the relationship between x and y given z1, z2, z3. All variables can be either a distance matrix of class dist or vectors of dissimilarities. }
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
  If only one independent variable is given, the simple Mantel r (r12) is calculated. If more than one independent variable is given, the partial Mantel r  (ryx|x1 ...) is calculated by permuting one of the original dissimilarity matrices.
The bootstrapping is actually resampling without replacement, because duplication of samples is not useful in a dissimilarity context (the dissimilarity of a sample with itself is zero). Resampling within dissimilarity values is inappropriate, just as for permutation.
}
\value{
  \item{mantelr }{Mantel coefficient.}
  \item{pval1 }{one-tailed p-value (null hypothesis: r <= 0).}
  \item{pval2 }{one-tailed p-value (null hypothesis: r >= 0).}
  \item{pval3 }{two-tailed p-value (null hypothesis: r = 0).}
  \item{llim }{lower confidence limit.}
  \item{ulim }{upper confidence limit.}
}
\references{ Mantel, N. 1967. The detection of disease clustering and a generalized regressio
n 
approach. Cancer Research 27:209-220.

Smouse, P.E., J.C. Long and R.R. Sokal. 1986. Multiple regression and correlatio
n 
extensions of the Mantel test of matrix correspondence. Systematic Zoology 35:62
7-632. }
\author{ Sarah Goslee, Sarah.Goslee@ars.usda.gov }

\seealso{ \code{\link{mgram}} }
\examples{
\dontrun{
# Example of multivariate analysis using built-in iris dataset
data(iris)
iris.md <- distance(iris[,1:4], "mahal")

# Create a model matrix for testing species differences
iris.model <- distance(as.numeric(iris[,5]), "eucl")
iris.model[iris.model > 0] <- 1

# Test whether samples within the same species are more similar than those not
mantel(iris.md ~ iris.model, nperm=10000)
}
}

\keyword{ multivariate }

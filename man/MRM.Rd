\name{MRM}
\alias{MRM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Multiple Regression on distance Matrices }
\description{
  Multiple regression on distance matrices (MRM) using permutation tests of significance for regression coefficients and R-squared.
}

\usage{
MRM(formula = formula(data), data = sys.parent(), nperm = 1000, method = "linear",
mrank = FALSE)
}

\arguments{
  \item{formula}{ formula in R/S-Plus format describing the test to be conducted. }
  \item{data}{ an optional dataframe containing the variables in the model as columns of dissimilarities. By default the variables are taken from the current environment. }
  \item{nperm}{ number of permutations to use. If set to 0, the permutation test will be omitted. }
  \item{mrank}{ if this is set to FALSE (the default option), Pearson correlations will be used. If
 set to TRUE, the Spearman correlation (correlation ranked distances) will be used.  }
  \item{method}{ if "linear", the default, uses multiple regression analysis. If "logistic", performs logistic regression with appropriate permutation testing. Note that this may be substantially slower.}
}
\details{
  Performs multiple regression on distance matrices following the methods outlined in Legendre et al. 1994.
}
\value{
  \item{coef }{A matrix with regression coefficients and associated p-values from the permutation test (using the pseudo-t of Legendre et al. 1994).}
  \item{r.squared }{Regression R-squared and associated p-value from the permutation test (linear only). } 
  \item{F.test }{F-statistic and p-value for overall F-test for lack of fit (linear only).}
  \item{dev }{Residual deviance, degrees of freetom, and associated p-value (logistic only).}
}

\references{ 
Lichstein, J. 2007. Multiple regression on distance matrices: A multivariate spatial analysis tool. Plant Ecology 188: 117-131.

Legendre, P.; Lapointe, F. and Casgrain, P. 1994. Modeling brain evolution from behavior: A permutational regression approach. Evolution 48: 1487-1499.
 }

\author{ Sarah Goslee, Sarah.Goslee@ars.usda.gov }
\seealso{ \code{\link{mantel}} }

\examples{
data(graze)
LOAR10.mrm <- MRM(dist(LOAR10) ~ dist(sitelocation) + dist(forestpct), data=graze, nperm=100)
}

\keyword{ multivariate }


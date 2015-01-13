\name{cor2m}
\alias{cor2m}
\title{ Generates a correlation table between the variables of 2 matrices }
\description{
  cor2m generates a correlation table between the variables of 2 matrices--by design, species and environment. It stores these in a table with species as columns and envnvironmental variables as rows, so it's easy to scan. It also removes correlations less than a user-specified alpha (0.05 by default) by setting them to NA.
}
\usage{
cor2m(x, y, trim = TRUE, alpha = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ A matrix of species (or other) variables }
  \item{y}{ A matrix of environmental (or other) variables matching the sites of x }
  \item{trim}{ If trim is TRUE, set rho<critical value(alpha) to 0 }
  \item{alpha}{ alpha p-value to use with trim }
}
\details{
Correlate variables from 1 matrix with variables from another matrix. If TRIM, set rho<critical value(alpha) to NA. Computes this critical value as a t-test with n-2 df.

}
\value{
  Returns a correlation table between the variables of 2 matrices.
  }
\author{ Dean Urban }
\examples{
# toy example
speciesdata <- matrix(runif(50), nrow=10)
colnames(speciesdata) <- c("spA", "spB", "spC", "spD", "spE")
envdata <- matrix(runif(30), nrow=10)
colnames(envdata) <- c("var1", "var2", "var3")
sppenv.cor <- cor2m(speciesdata, envdata)
print(sppenv.cor, na.print="")
}
\keyword{ multivariate }

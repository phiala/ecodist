\name{cor2m}
\alias{cor2m}
\title{ Two-matrix correlation table }
\description{
  Generate a correlation table between the variables of
  two data sets, originally for comparing species abundances and
  environmental variables.
}
\usage{
cor2m(x, y, trim = TRUE, alpha = 0.05)
}

\arguments{
  \item{x}{ A matrix or data frame of environmental (or other)
    variables matching the sites of x }
  \item{y}{ A matrix or data frame of species (or other) variables }
  \item{trim}{ If trim is TRUE, set rho<critical value(alpha) to 0 }
  \item{alpha}{ alpha p-value to use with trim, by default 0.05 }
}
\details{
  cor2m generates a correlation table between the variables of
  two matrices. The original use case is to compare species
  abundances and environmental variables. It results in a data frame
  with species (or the first matrix) as columns and environmental
  variables (or the second matrix) as rows, so it's easy to scan.
  Correlations less than a user-specified alpha (0.05 by default)
  can be set to NA.  cor2m generates a correlation table between
  the variables of two matrices. The original use case is to compare
  species abundances and environmental variables. The result has
  species (or the first matrix) as columns and environmental
  variables (or the second matrix) as rows, so it's easy to scan.
  Correlations less than a user-specified alpha can be set to NA.
  If trim, correlations less than the critical value for the
  provided alpha are set to to NA. The critical value is computed
  as a t-test with n-2 df.
  cor2m(x, y, trim=FALSE) is equivalent to cor(x, y)
}
\value{
  Returns a data frame of correlations between the variables
    of 2 data frames.
  }
\author{ Dean Urban }
\examples{
data(graze)
speciesdata <- graze[, 3:7]
envdata <- graze[, 1:2]
sppenv.cor <- cor2m(envdata, speciesdata)
print(sppenv.cor, na.print="")
}
\keyword{ multivariate }

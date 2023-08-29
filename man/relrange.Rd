\name{relrange}
\alias{relrange}
\title{
	Relativize a compositional data matrix.
}
\description{
   Relativizes the range of each column of a data frame or matrix x to 0-1. If globalmin and/or globalmax are provided, those are used to scale the columns, for instance to scale a subset to match a larger sample. If they are NA, the minimum and maximum values for each column are used.
}
\usage{
relrange(x, globalmin = NA, globalmax = NA)
}
\arguments{
  \item{x}{
	The data frame or matrix to be relativized.
}
  \item{globalmin}{
	A value other than the population minimum to be used. Should be the same length as the number of columns of x.
}
  \item{globalmax}{
	A value other than the population maximum to be used. Should be the same length as the number of columns of x.

}
}
\details{
	Relativizes the data using the minimum and maximum values. If globalmin and global max are not used, the range will be 0-1 for each variable. This can be useful for putting disparate variables to the same magnitude while keeping all non-negative values.
}
\value{

	Returns an object of the same class as x (matrix or data frame) with the columns rescaled.
}

\author{ Sarah Goslee }

\seealso{
	\code{\link{scale}}
}
\examples{

	x <- matrix(1:15, ncol = 3)

	# uses min and max of the data
	relrange(x)

	# uses min and max determined by other knowledge of the variables
	relrange(x, globalmin = c(0, 0, 0), globalmax = c(6, 10, 20))
}

\keyword{ multivariate }


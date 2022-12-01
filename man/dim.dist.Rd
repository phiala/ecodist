\name{dim.dist}
\alias{dim.dist}
\title{ Dimension of a distance object }
\description{
  Returns NULL for the dimensions of a distance object.
}
\usage{
  \method{dim}{dist}(x)
}
\arguments{
  \item{x}{ object of class \code{dist} }
}
\details{
The spdep package overwrites the base R behavior of dim.dist() to return c(n, n) where n is the size of the full matrix. The base R behavior returns NULL. This function restores base R behavior within ecodist, because otherwise spdep being loaded breaks ecodist functionality.
}
\value{
  NULL
}
\author{ Sarah Goslee }


\examples{

	data(graze)
	dim(dist(graze))

}

\keyword{ manip }
\keyword{ multivariate }

\name{addord}
\alias{addord}
\title{
Fit new points to an existing NMDS configuration.
}
\description{
Uses a brute force algorithm to find the location for each new point that minimizes overall stress.
}
\usage{
addord(origconf, fulldat, fulldist, isTrain, bfstep = 10, maxit = 50, epsilon = 1e-12)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{origconf}{
    The original ordination configuration.
}
  \item{fulldat}{
    The dataset containing original and new points.
}
  \item{fulldist}{
    A dissimilarity matrix calculated on \code{fulldat}.
}
  \item{isTrain}{
    A boolean vector of length \code{nrow(fulldat)} indicating which rows were training data used in determining \code{origconf} (TRUE), or are new points (FALSE).
}
  \item{bfstep}{
    A tuning parameter for the brute force algorithm describing the size of grid to use.
}
  \item{maxit}{
    The maximum number of iterations to use.	
}
  \item{epsilon}{
	Tolerance value for convergence.
}
}
\details{
   A region comprising the original ordination configuration plus one standard deviation is divided into a grid of \code{bfstep} rows and columns. For a new point, the grid cell with the lowest stress is identified. That cell is divided into a finer grid, and the lowest-stress cell identified. This process is repeated up to \code{maxit} times, or until stress changes less than \code{epsilon}.
}
\value{
  \item{fullfitconf }{The new ordination configuration containing training and new points.}
  \item{stress }{The stress value for each point.}
  \item{isTrain }{The boolean vector indicating training set membership, for reference.}
}
\author{
  Sarah Goslee <Sarah.Goslee@ars.usda.gov>
}

\examples{
\dontrun{
data(iris)
set.seed(123)
iris.sub <- sample(1:150, 15)

iris.istrain <- rep(TRUE, nrow(iris))
iris.istrain[iris.sub] <- FALSE

iris.full <- iris[,1:4]
iris.full.d <- distance(iris.full, "mahal")

iris.train.d <- distance(iris.full[iris.istrain, ], "mahal")
set.seed(456)
# use larger nits in practice
iris.train.nmds <- nmds.min(nmds(iris.train.d, mindim=2, maxdim=2, nits=3)) 

# use larger maxit in practice
iris.fit <- addord(iris.train.nmds, iris.full, iris.full.d, iris.istrain, maxit=2)

plot(iris.fit$conf, col=iris$Species, pch=c(18, 4)[iris.istrain + 1], xlab="NMDS 1", ylab="NMDS 2")
title("Demo: adding points to an ordination")
legend("bottomleft", c("Training set", "Added point"), pch=c(4, 18))
legend("topleft", levels(iris$Species), fill=1:3)

}
}

\keyword{ multivariate }


\name{xdistance}
\alias{xdistance}
\title{
Cross-distance between two datasets.
}
\description{
Pairwise dissimilarity calculation between rows of one dataset and rows of another, for instance across different sampling periods for the same set of sites.
}
\usage{
xdistance(x, y, method = "euclidean")
}

\arguments{
  \item{x}{A site by species or other matrix or data frame.}
  \item{y}{A a second site by species dataset, which must have at least the same columns.}
  \item{method}{This function calls \code{\link{distance}} to do the calculations, and will accept any method used there, currently: "euclidean", "bray-curtis", "manhattan", "mahalanobis" (squared Mahalanobis distance), "jaccard", "difference", "sorensen", "gower", "modgower10" (modified Gower, base 10), "modgower2" (modified Gower, base 2). Partial matching will work for selecting a method.}
}
\details{
This function will calculate rowwise dissimilarities between any pair of matrices or data frames with the same number of columns.
  Note that the cross-dissimilarity functions are for research purposes, and are not well-tested.
}
\value{
A non-symmetric and possibly not square matrix of dissimilarities of class \code{xdist}, where \code{result <- xdistance(x, y)} produces a matrix with \code{result[a, b]} containing the dissimilarity between \code{x[a, ]} and \code{y[b, ]}.
}
\author{ Sarah Goslee }

\seealso{ \code{\link{distance}}, \code{\link{xmantel}}, \code{\link{xmgram}}}
\examples{
data(graze)

### EXAMPLE 1: Square matrices

# take two subsets of sites with different dominant grass abundances
# use cut-offs that produce equal numbers of sites
dom1 <- subset(graze, POPR > 50 & DAGL < 20) #  8 sites
dom2 <- subset(graze, POPR < 50 & DAGL > 20) #  8 sites

# first two columns are site info
dom.xd <- xdistance(dom1[, -c(1,2)], dom2[, -c(1,2)], "bray")

# environmental and spatial distances; preserve rownames
forest.xd <- xdistance(dom1[, "forestpct", drop=FALSE], 
    dom2[, "forestpct", drop=FALSE])
sitelocation.xd <- xdistance(dom1[, "sitelocation", drop=FALSE], 
    dom2[, "sitelocation", drop=FALSE])

# permutes rows and columns of full nonsymmetric matrix
xmantel(dom.xd ~ forest.xd)
xmantel(dom.xd ~ forest.xd + sitelocation.xd)

plot(xmgram(dom.xd, sitelocation.xd))


### EXAMPLE 2: Non-square matrices

# take two subsets of sites with different dominant grass abundances
# this produces a non-square matrix

dom1 <- subset(graze, POPR > 45 & DAGL < 20) # 13 sites
dom2 <- subset(graze, POPR < 45 & DAGL > 20) #  8 sites

# first two columns are site info
dom.xd <- xdistance(dom1[, -c(1,2)], dom2[, -c(1,2)], "bray")

# environmental and spatial distances; preserve rownames
forest.xd <- xdistance(dom1[, "forestpct", drop=FALSE], 
    dom2[, "forestpct", drop=FALSE])
sitelocation.xd <- xdistance(dom1[, "sitelocation", drop=FALSE], 
    dom2[, "sitelocation", drop=FALSE])

# permutes rows and columns of full nonsymmetric matrix
xmantel(dom.xd ~ forest.xd, dims=c(13, 8))
xmantel(dom.xd ~ forest.xd + sitelocation.xd, dims=c(13, 8))

plot(xmgram(dom.xd, sitelocation.xd))
}

\keyword{ multivariate }

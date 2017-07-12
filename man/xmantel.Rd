\name{xmantel}
\alias{xmantel}
\title{ Cross-Mantel test }
\description{
  Simple and partial cross-Mantel tests, with options for ranked data and permutation tests.
}
\usage{
xmantel(formula = formula(data), data = sys.parent(), dims = NA,
   nperm = 1000, mrank = FALSE)
}
\arguments{
  \item{formula}{ formula describing the test to be conducted. For this test, y ~ x will perform a simple Mantel test, while y ~ x + z1 + z2 + z3 will do a partial Mantel test of the relationship between x and y given z1, z2, z3. All variables should be either non-symmetric square cross-dissimilary matrices of class xdist, or vector forms thereof. }
  \item{data}{ an optional dataframe containing the variables in the model as columns of dissimilarities. By default the variables are taken from the current environment. }
  \item{dims}{ if the dissimilarity matrices are not square, the dimensions must be provided as \code{c(nrow, ncol)}}
  \item{nperm}{ number of permutations to use. If set to 0, the permutation test will be omitted. }
  \item{mrank}{ if this is set to FALSE (the default option), Pearson correlations will be used. If
 set to TRUE, the Spearman correlation (correlation ranked distances) will be used. }
}
\details{
  If only one independent variable is given, the simple Mantel r (r12) is calculated. If more than one independent variable is given, the partial Mantel r  (ryx|x1 ...) is calculated by permuting one of the original dissimilarity matrices.
  Note that the cross-dissimilarity functions are for research purposes, and are not well-tested.

}
\value{
  \item{mantelr }{Mantel coefficient.}
  \item{pval1 }{one-tailed p-value (null hypothesis: r <= 0).}
  \item{pval2 }{one-tailed p-value (null hypothesis: r >= 0).}
  \item{pval3 }{two-tailed p-value (null hypothesis: r = 0).}
}

\author{ Sarah Goslee, Sarah.Goslee@ars.usda.gov }

\seealso{  \code{\link{xdistance}}, \code{\link{xmgram}} }
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

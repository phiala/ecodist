\name{xmgram}
\alias{xmgram}
\title{Cross-Mantel correlogram }
\description{
	Calculates simple Mantel correlograms from cross-distance matrices.
}
\usage{
xmgram(species.xd, space.xd, breaks, nclass, stepsize, equiprobable = FALSE, nperm = 1000,
    mrank = FALSE, alternative = "two.sided", trace = FALSE)
}
\arguments{
  \item{species.xd}{ non-symmetric square cross-distance matrix. }
  \item{space.xd}{ non-symmetric square matrix of geographic distances. }
  \item{breaks}{ locations of class breaks. If specified, overrides nclass and stepsize. }
  \item{nclass}{ number of distance classes. If not specified, Sturge's rule will be used
to determine an appropriate number of classes. }
  \item{stepsize}{ width of each distance class. If not specified, nclass and the range of space.d will be used to calculate an appropriate default. }
  \item{equiprobable}{ if TRUE, create nclass classes of equal number of distances; if FALSE, create nclass classes of equal width }  
  \item{nperm}{ number of permutations to use. If set to 0, the permutation test will be omitted. }

  \item{mrank}{ if this is set to FALSE (the default option), Pearson correlations will be used. If
 set to TRUE, the Spearman correlation (correlation ranked distances) will be used. }
  \item{alternative}{ default is "two.sided", and returns p-values for H0: rM = 0. The alternative is "one.sided", which returns p-values for H0: rM <= 0. }
  \item{trace}{ if TRUE, returns progress indicators. }
}
\details{
This function calculates cross-Mantel correlograms. The Mantel correlogram is essentially a multivariate autocorrelation function.
The Mantel r represents the dissimilarity in variable composition  (often
species composition) at a particular lag distance.
  Note that the cross-dissimilarity functions are for research purposes, and are not well-tested.
}
\value{
	Returns an object of class mgram, which is a list with two elements.
    mgram is a matrix with one row for each distance class and 6 columns:
  \item{lag }{midpoint of the distance class.}
  \item{ngroup }{number of distances in that class.}
  \item{mantelr }{Mantel r value.}
  \item{pval }{p-value for the test chosen.}
  resids is NA for objects calculated by mgram().
}
\references{ Legendre, P. and M. Fortin. 1989. Spatial pattern and ecological analysis.
Vegetatio 80:107-138. }
\author{ Sarah Goslee }

\seealso{\code{\link{xdistance}} \code{\link{xmantel}}, \code{\link{plot.mgram}} }
\examples{

# Need to develop a cross-dissimilarity example
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

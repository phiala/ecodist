\name{permresults}
\alias{permresults}
\docType{data}
\title{ Mantel test results for number of permutations, correlation and dataset size.}
\description{
  These results explore the variability in randomization tests of significance and in bootstrap confidence limits for different sample sizes, strengths of relationship and numbers of permutations.
}
\usage{data(permresults)}
\format{
  A data frame with 15000 observations on the following 17 variables.
  \describe{
    \item{\code{n}}{The number of samples.}
    \item{\code{r}}{The Mantel r of the dissimilarity vectors.}
    \item{\code{p500}}{The p-values generated using 500 permutations.}
    \item{\code{l100}}{The lower confidence limit generated using 100 subsamples.}
    \item{\code{u100}}{The upper confidence limit generated using 100 subsamples.}
    \item{\code{p1000}}{p-values with 1,000 permutations}
    \item{\code{l500}}{lower CI with 500 permutations}
    \item{\code{u500}}{upper CI with 500 permutations}
    \item{\code{p10000}}{p-values with 10,000 permutations}
    \item{\code{l1000}}{lower CI with 1,000 permutations}
    \item{\code{u1000}}{upper CI with 1,000 permutations}
    \item{\code{p100000}}{p-values with 100,000 permutations}
    \item{\code{l10000}}{lower CI with 10,000 permutations}
    \item{\code{u10000}}{upper CI with 10,000 permutations}
    \item{\code{p1000000}}{p-values with 1,000,000 permutations}
    \item{\code{l100000}}{lower CI with 100,000 permutations}
    \item{\code{u100000}}{upper CI with 100,000 permutations}
  }
}
\details{
These simulation results are included here because the entire simulation takes an extremely long time to run. The tests explored how the p-value and confidence limits varied with:
\preformatted{
- number of samples (10 (45 dissimilarities), 25 (300), 50 (1225), 100 (4950))
- number of permutations: 500, 1000, 10000, 100000, 1000000
- number of bootstrap samples: 100, 500, 1000, 10000, 100000
- correlation of the original data: 0.01, 0.05, 0.10, 0.25, 0.50
}
The larger correlation values were omitted for the larger sample sizes because they were always highly signficant, and an additional low value of 0.025 was added for the largest sample size.

The method used to generate the data is shown in the example section.


}
\examples{
data(permresults)

\dontrun{
# Simulation Method
# example for sample size 10, correlation of 0.01

permresults <- data.frame(matrix(0, nrow=1000, ncol=17))
colnames(permresults) <- c("n", "r", "p500", "l100", "u100", "p1000", 
    "l500", "u500", "p10000", "l1000", "u1000", "p100000", "l10000", "u10000", 
    "p1000000", "l100000", "u100000")

# generate correlated data
corxy.len010.r010 <- corgen(len= 10, r= 0.01, epsilon=0.0001)

all.nperm <- c(500, 1000, 10000, 100000, 1000000)
all.nboot <- c(100, 500, 1000, 10000, 100000)

for(i in 1:1000) {
    permresults[i, 1] <- 10
    permresults[i, 2] <- 0.01
    for(j in 1:5) {
# record p-value, lower and upper confidence limits
        thismantel <- mantel(corxy.len010.r010$y ~ corxy.len010.r010$x, 
		  nperm=all.nperm[j], nboot=all.nboot[j])
        permresults[i, (3*j):(3*j + 2)] <- c(thismantel[c(2,5,6)])
    }
}
}

}
\keyword{datasets}

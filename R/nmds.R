# Non-metric multidimensional scaling function
# using the majorization algorithm from
# Borg & Groenen 1997, Modern Multidimensional Scaling.
#
# Sarah Goslee
# 20 December 1999
#
# dmat is a lower-triangular distance matrix.
# mindim is the minimum number of dimensions to calculate.
# maxdim is the maximum number of dimensions to calculate.
# nits is the number of repetitions.
# iconf is the initial configuration to use.
# If iconf is not specified, then a random configuration is used.
# epsilon and maxit specify stopping points.
# Returns a list of configurations (conf)
# and a vector of final stress values (stress),
# along with the cumulative and incremental r^2 for each axis.
# The first nits elements are for the lowest number of dimensions.
# mindim, maxdim, nits are saved as part of the returned list
#
# stresscalc has been updated 2016-12-27 to be compatible with vegan::metaMDS
# and MASS::isoMDS.
# The method of finding the optimum is unchanged.

sstress <- function(dmat, newconf)
{
    # Calculates the stress-1 function for the original and
    # new NMDS configurations (Kruskal 1964).

    # Calculate stress based on isotonic regression
    # of distance matrices
    # uses stats::isoreg instead of MASS::Shepard

    dmat <- as.vector(dmat)
    dord <- order(dmat)
    cmat <- as.vector(dist(newconf, "minkowski", 2))

    dmat <- dmat[dord]
    cmat <- cmat[dord]

    ir <- isoreg(dmat, cmat)
    dmat <- ir$y
    cmat <- ir$yf

    sstresscalc <- (dmat - cmat) ^ 2
    sstresscalc <- sum(sstresscalc) / sum(dmat ^ 2)

    sqrt(sstresscalc)
}

nmdscalc <- function(dmat, ndim, iconf, epsilon, maxit, trace) {

    # This is the optimization routine for NMDS ordination.
    # Use front-end nmds.
    # not exported

    n <- (1 + sqrt(1 + 8 * length(dmat))) / 2

    if(missing(iconf)) {
       iconf <- matrix(runif(n * ndim), nrow = n, ncol = ndim)
    }

    k <- 0
    conf <- iconf
    stress2 <- sstress(dmat, iconf)
    stress1 <- stress2 + 1 + epsilon

    while(k < maxit && abs(stress1 - stress2) > epsilon) {
       stress1 <- stress2

       dmat.full <- full(dmat)
       confd.full <- full(dist(conf))

       confd2.full <- confd.full
       confd2.full[confd.full == 0] <- 1

       b <- dmat.full / confd2.full
       b[confd.full == 0] <- 0
       bsum <- apply(b, 2, sum)
       b <- -1 * b
       diag(b) <- bsum

       conf <- (1 / n) * b %*% conf

       stress2 <- sstress(dmat, conf)

       if(trace) cat(k, ",\t", stress1, "\n")

       k <- k + 1
    }

    list(conf = conf, stress = stress1)

}

nmds <- function(dmat, mindim = 1, maxdim = 2, nits = 10, iconf, epsilon = 1e-12, maxit = 500, trace=FALSE) {

    conf <- list(1:((maxdim - mindim + 1) * nits))
    stress <- list(1:((maxdim - mindim + 1) * nits))
    r2 <- list(1:((maxdim - mindim + 1) * nits))

    k <- 1

    for(i in mindim:maxdim) {
       if(trace) cat("Number of dimensions: ", i, "\n")
       for(j in seq_len(nits)) {
          if(trace) cat("Iteration number: ", j, "\n")
            if(!missing(iconf) && ncol(iconf) == i) {
                nmdsr <- nmdscalc(dmat, ndim = i, iconf=iconf, epsilon=epsilon, maxit=maxit, trace=trace)
            } else {
                nmdsr <- nmdscalc(dmat, ndim = i, epsilon=epsilon, maxit=maxit, trace=trace)
            }

          conf[[k]] <- nmdsr$conf
          stress[[k]] <- nmdsr$stress
          r2[[k]] <- cor(dmat, dist(nmdsr$conf)) ^ 2
          k <- k + 1
       }
    }

    results <- list(conf = conf, stress = unlist(stress), r2 = unlist(r2), mindim=mindim, maxdim=maxdim, nits=nits)

    class(results) <- "nmds"
    results
}


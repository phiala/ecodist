xmgram <- function(species.xd, space.xd, breaks, nclass, stepsize, nperm = 1000, mrank = FALSE, alternative = "two.sided", trace = FALSE)


# Cross-Mantel correlogram developed from mgram()
# Sarah Goslee 2017-02-17
#
# This function calculates a mantel correlogram for a full cross-distance (non-symmetric
# but square) based on the geographic distances given in space.d, also a full cross-distance
# (nonsymmetric but square) distance matrix. 
# nclass: number of distance classes 
# stepsize: width of distance classes
# nperm: number of permutations for mantel test
#
# Default is two-tailed test (H0: rM = 0; alternative = "two.sided")
# May also use one-sided test (H0: rM <= 0; alternative = "one.sided")

{

    dims <- dim(species.xd)

    space.xd <- as.vector(space.xd)

# use breaks if it exists.
# If nclass or stepsize aren't specified, use Sturge's rule to calculate nclass
# classes are shifted so that they don't have to start with zero
    if(missing(breaks)) {
        if (missing(nclass)) {
            if (missing(stepsize)) {
                nclass <- round(1 + 3.3 * log10(length(space.xd)))
                stepsize <- (max(space.xd) - min(space.xd))/nclass
            }
            else {
                nclass <- round((max(space.xd) - min(space.xd))/stepsize)
            }
        }
        else {
            if (missing(stepsize)) {
                stepsize <- (max(space.xd) - min(space.xd))/nclass
            }
        }
        breaks <- seq(0, stepsize * nclass, stepsize)
    }
    else {
        nclass <- length(breaks) - 1
    }

    answer.m <- matrix(0, ncol=4, nrow=nclass)
    dimnames(answer.m) <- list(NULL, c("lag", "ngroup", "mantelr", "pval"))
    answer.m[,4] <- rep(1, nrow(answer.m))

    for(i in seq_len(nclass)) {
      dmin <- breaks[i]
      dmax <- breaks[i + 1]
        answer.m[i,1] <- (dmin + dmax) / 2

        space.dclass <- rep(0, length(space.xd))
        space.dclass[space.xd <= dmin] <- 1
        space.dclass[space.xd > dmax] <- 1

        ngroup <- length(space.dclass) - sum(space.dclass)
        answer.m[i,2] <- ngroup

        if(ngroup > 0) {
            space.dclass <- matrix(space.dclass, nrow=dims[1], ncol=dims[2])
            mant <- xmantel(species.xd ~ space.dclass, dims=dims, nperm=nperm, mrank=mrank)
            answer.m[i,3] <- mant[1]
            if(alternative == "two.sided")
                answer.m[i,4] <- mant[4]
            else
                answer.m[i,4] <- mant[2]
        }            
            
        if(trace) cat(i, "\t", answer.m[i,2], "\t", answer.m[i, 3], "\n")    

    }

   results <- list(mgram = answer.m, resids = NA)
   class(results) <- "mgram"
   results
}


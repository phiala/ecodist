pmgram <- function(data, space, partial, breaks, nclass, stepsize, equiprobable = FALSE, resids = FALSE, nperm = 1000)

# Piecewise Mantel correlogram
# 
# results are:
#    1 dmax
#    2 ngroup
#    3 piecewise mgram
#    4 two-sided p-value

# use breaks if it exists.
# If nclass or stepsize aren't specified, use Sturge's rule to calculate nclass
# classes are shifted so that they don't have to start with zero

# 2023-08-29 added equiprobable option for distance classes of equal number
# rather than equal width

{

    epsilon <- 0.0000001

    if(inherits(data, "dist")) {
        data <- as.matrix(as.vector(data))
    } else {
        data <- as.matrix(data)
    }
    space <- as.vector(space)
    n <- ncol(full(space))

# 2023-08-18: changed from using round to ceiling in Sturges' rule 
# calculation for compatibility with nclass.Sturges
    if(missing(breaks)) {
      if(missing(nclass)) {
         if(missing(stepsize)) {
            nclass <- ceiling(1 + log2(length(space)))
            stepsize <- (max(space) - min(space)) / nclass
         } else {
            nclass <- round((max(space) - min(space))/stepsize)
         }
      } else {
         if(missing(stepsize)) {
            stepsize <- (max(space) - min(space)) / nclass
         }
      }
           if(equiprobable) {
            breaks <- quantile(space, seq(0, 1, length.out = nclass + 1))
        } else {
            # equal width breaks
            breaks <- seq(0, stepsize * nclass, stepsize)
        }
    } else {
        nclass <- length(breaks) - 1
    }

    answer.m <- matrix(NA, nrow=nclass, ncol=4)
    dimnames(answer.m) <- list(NULL, c("lag", "ngroup", "pieceR", "pval"))
    answer.m[,4] <- rep(NA, nrow(answer.m))

# standardize so mean = 0, variance = 1
    for(i in seq_len(ncol(data))) {
        thiscol <- data[,i]
        ydif <- thiscol - mean(thiscol)
        yvar <- sqrt(sum(ydif^2)/length(ydif))
        thiscol <- ydif / yvar
        data[,i] <- thiscol
    }

    if(!missing(partial)) {
        partial <- as.matrix(as.vector(partial))
        for(i in seq_len(ncol(partial))) {
            thiscol <- partial[,i]
            ydif <- thiscol - mean(thiscol)
            yvar <- sqrt(sum(ydif^2)/length(ydif))
            thiscol <- ydif / yvar
            partial[,i] <- thiscol
        }
    }

   if(resids) {
      mgresids <- rep(0, length(space))
   }
   else {
      mgresids <- NA
   }

    if(missing(partial)) {
        if(ncol(data) == 1) {
            colnames(answer.m)[3] <- "wtI"
            for(i in seq_len(nclass)) {
            dmin <- breaks[i]
            dmax <- breaks[i + 1]
         
                answer.m[i,1] <- (dmin + dmax) / 2
    
                space.dclass <- rep(0, length(space))
                space.dclass[space <= dmin] <- 1
                space.dclass[space > dmax] <- 1

                if(sum(space.dclass== 0) > 0) {
                    ngroup <- length(space.dclass) - sum(space.dclass)
                    answer.m[i,2] <- ngroup

                    answer.m[i,3] <- (-1/ngroup) * sum(data[space.dclass == 0]) # similar to Moran's I
                    if(resids == TRUE) {
                        mgresids[space.dclass == 0] <- residuals(lm(data[space.dclass == 0] ~ space[space.dclass == 0]))
                    }
                    if(nperm > 0) {
                        zstats <- rep(0, nperm)
                        tmat <- matrix(0, n, n)
                        rarray <- rep(0, n)
        
    #                    cat(i, "\t", answer.m[i,1], "\t", answer.m[i,2], "\t", answer.m[i,3], "\n")
                        cresults <- .C("newpermone",
                            as.double(data),
                            as.integer(space.dclass),
                            as.integer(n),
                            as.integer(length(data)),
                            as.integer(nperm),
                            zstats = as.double(zstats),
                            as.double(as.vector(tmat)),
                            as.integer(rarray),
                            PACKAGE = "ecodist")
            
                        zstats <- cresults$zstats
                        answer.m[i,4] <- length(zstats[abs(zstats) >= abs(zstats[1])])/nperm
                    }
                }
            }
        }
        else {
            for(i in seq_len(nclass)) {
            dmin <- breaks[i]
            dmax <- breaks[i + 1]
                answer.m[i,1] <- dmax
    
                space.dclass <- rep(0, length(space))
                space.dclass[space <= dmin] <- 1
                space.dclass[space > dmax] <- 1

                if(sum(space.dclass== 0) > 0) {
                    ngroup <- length(space.dclass) - sum(space.dclass)
                    answer.m[i,2] <- ngroup

                    if(ncol(data) == 2) {
                        answer.m[i, 3] <- cor(data[space.dclass == 0, 1], data[space.dclass == 0, 2])
                        if(resids == TRUE) {
                            mgresids[space.dclass == 0] <- residuals(lm(data[space.dclass == 0,1] ~ data[space.dclass == 0,2]))
                        }
                        if(nperm > 0) {
                            if(is.na(answer.m[i, 3])) 
                                answer.m[i, 4] <- 1
                            else {
                                zstats <- rep(0, nperm)
                                tmat <- matrix(0, n, n)
                                rarray <- rep(0, n)
                                xmat <- data[,1]
                                xmat[space.dclass == 1] <- -9999

                            #    cat(i, "\t", answer.m[i,1], "\t", answer.m[i,2], "\t", answer.m[i,3], "\n")
                                cresults <- .C("newpermtwo",
                                    as.double(xmat),
                                    as.double(data[,2]),
                                    as.integer(n),
                                    as.integer(length(xmat)),
                                    as.integer(nperm),
                                    zstats = as.double(zstats),
                                    as.double(as.vector(tmat)),
                                    as.integer(rarray),
                                    PACKAGE = "ecodist")
            
                                zstats <- cresults$zstats
                                answer.m[i,4] <- length(zstats[abs(zstats) >= abs(zstats[1])])/nperm
                            }
                        }
                    }
                }    
            }
        }
    }
    else {
        if(ncol(data) == 1) {
            colnames(answer.m)[3] <- "wtI"
            for(i in seq_len(nclass)) {
            dmin <- breaks[i]
            dmax <- breaks[i + 1]
                answer.m[i,1] <- dmax
    
                space.dclass <- rep(0, length(space))
                space.dclass[space <= dmin] <- 1
                space.dclass[space > dmax] <- 1

                ngroup <- sum(space.dclass== 0)
                answer.m[i,2] <- ngroup

                if(ngroup > 0) {

                    ngroup <- sum(space.dclass== 0)

                    data.lm <- residuals(lm(data ~ partial))    
                    answer.m[i,3] <- (-1/ngroup) * sum(data.lm[space.dclass == 0]) # similar to Moran's I
                    if(resids == TRUE) {
                        mgresids[space.dclass == 0] <- residuals(lm(data[space.dclass == 0] ~ partial[space.dclass == 0]))
                    }
                    if(nperm > 0) {
                        zstats <- rep(0, nperm)
                        zstats[1] <- answer.m[i,3]
                        for(j in 2:nperm) {
                            xmat <- full(data)
                            xsamp <- sample(ncol(xmat))
                            xmat <- xmat[xsamp, xsamp]
                            xmat <- lower(xmat)
                            data.lm <- residuals(lm(xmat ~ partial))
                            zstats[j] <- (-1/ngroup) * sum(data.lm[space.dclass == 0])
                        }
                        answer.m[i,4] <- length(zstats[abs(zstats) >= abs(zstats[1])])/nperm
                    }
                }
            }
        }
        else {
            for(i in seq_len(nclass)) {
            dmin <- breaks[i]
            dmax <- breaks[i + 1]
                answer.m[i,1] <- (dmin + dmax) / 2
    
                space.dclass <- rep(0, length(space))
                space.dclass[space <= dmin] <- 1
                space.dclass[space > dmax] <- 1
                
                if(sum(space.dclass== 0) > 0) {

                    ngroup <- length(space.dclass) - sum(space.dclass)
                    answer.m[i,2] <- ngroup

                    if(ncol(data) == 2) {
                        data1.lm <- residuals(lm(data[space.dclass == 0,1] ~ partial[space.dclass == 0,]))
                        data2.lm <- residuals(lm(data[space.dclass == 0,2] ~ partial[space.dclass == 0,]))
                        answer.m[i, 3] <- cor(data1.lm, data2.lm)
                        if(resids == TRUE) {
                            mgresids[space.dclass == 0] <- residuals(lm(data1.lm ~ data2.lm))
                        }
                        if(nperm > 0) {
                            if(is.na(answer.m[i, 3])) 
                                answer.m[i, 4] <- 1
                            else {
                                zstats <- rep(0, nperm)
                                zstats[1] <- answer.m[i,3]
                                for(j in 2:nperm) {
                                    xmat <- data[,1]
                                    xmat[space.dclass == 1] <- -9999
                                    xmat <- full(xmat)
                                    xsamp <- sample(ncol(xmat))
                                    xmat <- xmat[xsamp, xsamp]
                                    xmat <- lower(xmat)
                                    xmat <- xmat[xmat != -9999]
                                    data1.lm <- residuals(lm(xmat ~ partial[space.dclass == 0,]))
                                    zstats[j] <- cor(data1.lm, data2.lm)
                                }
                                answer.m[i,4] <- length(zstats[abs(zstats) >= abs(zstats[1])])/nperm
                            }
                        }
                    }
                }    
            }
        }
    }

   results <- list(mgram = data.frame(answer.m), resids = mgresids)
   class(results) <- "mgram"
   results
}


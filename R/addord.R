addord <- function(origconf, fulldat, fulldist, isTrain, bfstep=10, maxit = 50, epsilon = 1e-12) {

    ## add new points to an ordination configuration
    ## by minimizing stress of new point location
    ## brute force, similar to PC-ORD method

    setgrid <- function(minvec, maxvec, nstep) {
        # find midpoints for each
        midlist <- vector(mode="list", length=length(minvec))
        midlist <- lapply(seq_along(minvec), function(i){
            thisstep <- (maxvec[i] - minvec[i])/nstep
            seq(minvec[i]+thisstep/2, maxvec[i]-thisstep/2, length=nstep)
        })

        as.matrix(expand.grid(midlist))
    }

    ##########
    if(missing(isTrain)) {
        isTrain <- c(rep(TRUE, nrow(origconf)), rep(FALSE, nrow(fulldat) - nrow(origconf)))
    }

    # set up extent to sample: +/- 1 sd
    osd <- apply(origconf, 2, sd)
    omin <- apply(origconf, 2, min) - osd
    omax <- apply(origconf, 2, max) + osd

    startgrid <- setgrid(omin, omax, bfstep)
    colnames(startgrid) <- colnames(origconf)

    # set up output objects
    fullfitconf <- data.frame(matrix(NA, nrow=nrow(fulldat), ncol=ncol(origconf)))
    colnames(fullfitconf) <- colnames(origconf)
    fullfitconf[isTrain, ] <- origconf

    stress.fullfit <- rep(NA, nrow(fullfitconf)) 

    for(thispoint in seq_along(isTrain)) {
        if(!isTrain[thispoint]) {

            # fit a new point to the ordination configuration by brute force
            # 1. make a distance matrix that includes selected new point as the last row/col

            usethis <- isTrain + 0
            usethis[thispoint] <- 2
            # orig points are 1, focus point is 2, points to not use are 0

            pointdist <- full(fulldist)[usethis > 0, usethis > 0]

            useord <- order(usethis[usethis > 0])

            pointdist <- pointdist[useord, useord]
            pointdist <- lower(pointdist)

            # 2. which of the startgrid gives the lowest stress?

            pointstress <- apply(startgrid, 1, function(x) {
                pointconf <- rbind(origconf, x)
                sstress(pointdist, pointconf)
            })

            conf <- startgrid[which.min(pointstress), ]
            stress2 <- min(pointstress)
            stress1 <- stress2 + 10 * epsilon # for stress, decreases
            thisstep <- (omax - omin) / (bfstep * 2)
            k <-  0

            while(k < maxit && abs(stress1 - stress2) > epsilon) {
                # go finer and finer into the ordination space

                newmin <- as.vector(conf) - thisstep
                newmax <- as.vector(conf) + thisstep

                newgrid <- setgrid(newmin, newmax, bfstep)
                colnames(newgrid) <- colnames(origconf)

                # find stress for newgrid

                pointstress <- apply(newgrid, 1, function(x) {
                    pointconf <- rbind(origconf, x)
                    sstress(pointdist, pointconf)
                })

                conf <- newgrid[which.min(pointstress), ]
                stress1 <- stress2
                stress2 <- min(pointstress)
                thisstep <- (newmax - newmin) / (bfstep * 2)
                k <- k + 1
            }
            fullfitconf[thispoint, ] <- conf
            stress.fullfit[thispoint] <- stress2
        }
    }

    list(conf=fullfitconf, stress=stress.fullfit, isTrain=isTrain)
}


xmantel <- function(formula = formula(data), data = sys.parent(), nperm = 1000, mrank = FALSE, nboot = 500, pboot = 0.90, cboot = 0.95)
{
# Cross-mantel test 
# Written by Sarah C. Goslee
# 01/01/01
# Updated 6 April 2001
#
# formula is y ~ x + n1 + n2 + n3 + ...
# NOT y ~ x | n1 + n2 + n3 + ... 
# The | means something else in S-Plus formulas.
#
# Uses C for permutation and bootstrap routines.
#
# Cross-mantel is a variation of the standard Mantel test that
# takes nonsymmetric full matrices such as those returned by
# xdist and xbcdist.
#
# This version calculates partial coefficients by permuting the y matrix.
#
# Will calculate the simple correlation or n-th order partial correlation  
# between two distance matrices in either of two ways: Pearson (mrank=F) 
# or Spearman (mrank=T)
#
# A permutation test is used to calculate the significance of r.
# The permutation test was designed to be relatively fast, but because of the
# way this was done, there is a possibility of repeating permutations of
# 1/n! where the distance matrix is n by n. In particular, for small matrices 
# n < 8 or so, it may be better to enumerate the permutations. 
#
#
# As an added bonus, this function offers the option of calculating 
# bootstrapped confidence limits for the correlation coefficient.
# nboot is the number of iterations.
# pboot is the level to resample at.
# cboot is the desired confidence limit.
#
# xmantel returns a five-element list:
# mantelr is the correlation.
# pval1 is the one-sided p-value (null hypothesis r <= 0) (0 if nperm == 0).
# pval2 is the one-sided p-value (null hypothesis r >= 0) (0 if nperm == 0).
# pval3 is the two-sided p-value (null hypothesis r = 0) (0 if nperm == 0).
# llim is the lower confidence limit.
# ulim is the upper confidence limit.


# Stuff R needs to be able to use a formula
        call <- match.call()
        m <- match.call(expand.dots = FALSE)
        m$family <- m$model <- m$... <- NULL
        m[[1]] <- as.name("model.frame")
	m <- m[1:2]
	m <- eval(m, sys.parent())
# End of S-Plus stuff. m is now the data for the Mantel test as
# columns y, x, n1, n2, n3, ...

# Determine the size of the matrices & do some error checking.

	n <- sqrt(dim(m)[[1]])
    if(abs(n - round(n)) > 0.0000001)
		stop("Matrix not square.\n")
	n <- round(n)
	if(dim(m)[[2]] < 2)
		stop("Not enough data. \n")

# If there are only x and y, then use the data as is.

	if(dim(m)[[2]] == 2) {
		ymat <- m[,1]
		xmat <- m[,2]
		if(mrank) {
			ymat <- rank(ymat)
			xmat <- rank(xmat)
		}
		ycor <- ymat
		xcor <- xmat
	} else {

# If this is a partial Mantel test, get the regression residuals
# for y ~ n1 + n2 + n3 + ... and x ~ n1 + n2 + n3 + ...
	    ymat <- as.vector(m[, 1])
        omat <- m[, -1]
            if(mrank) {
                    ymat <- rank(ymat)
                    omat <- apply(omat, 2, rank)
            }
            omat <- cbind(rep(1, length(ymat)), omat)
            xmat <- as.vector(omat[, 2])
            omat <- omat[, -2]
    omat <- as.matrix(omat)
            ycor <- lm.fit(omat, ymat)$residuals
            xcor <- lm.fit(omat, xmat)$residuals

	}

# Calculate the Mantel r

	mantelr <- cor(xcor, ycor)

# Standardize the columns of the matrices so
# that z = r and we can do 2-tailed tests.
                ncor <- length(xmat)

		w1 <- sum(xmat) / ncor
		w2 <- sum(xmat ^ 2)
		w2 <- sqrt(w2 / ncor - w1 ^ 2)
		xmat <- (xmat - w1) / w2

		w1 <- sum(ymat) / ncor
		w2 <- sum(ymat ^ 2)
		w2 <- sqrt(w2 / ncor - w1 ^ 2)
		ymat <- (ymat - w1) / w2

		if(dim(m)[[2]] > 2) {
			for(i in 2:dim(omat)[[2]]){
				curcoll <- omat[,i]
				w1 <- sum(curcoll) / ncor
				w2 <- sum(curcoll ^ 2)
				w2 <- sqrt(w2 / ncor - w1 ^ 2)
				curcoll <- (curcoll - w1) / w2
				omat[,i] <- curcoll
			}
		}

# If using a permutation test, start here:
	if(nperm > 0) {

# Set up the arrays needed.
		zstats <- numeric(nperm)
		tmat <- xmat
		rarray <- rep(0, n)      


		if(dim(m)[[2]] == 2) {
             cresults <- .C("xpermute",
                    as.double(xmat),
                    as.double(ymat),
                    as.integer(n),
                    as.integer(length(xmat)),
                    as.integer(nperm),
                    zstats = as.double(zstats),
                    as.double(tmat),
                    as.integer(rarray),
                    PACKAGE = "ecodist")

		} else {
			hmat <- solve((t(omat) %*% omat))
			hmat <- omat %*% hmat %*% t(omat)
			hmat <- diag(dim(hmat)[[1]]) - hmat
			xcor <- as.vector(lm.fit(omat, xmat)$residuals)
			ycor <- rep(0, length(xcor))

			cresults <- .C("xpermpart",
				as.double(as.vector(hmat)),
				as.double(as.vector(ymat)),
				as.double(xcor),
				as.double(ycor),
				as.integer(n),
				as.integer(length(xmat)),
				as.integer(nperm),
				zstats = as.double(zstats),
				as.double(as.vector(tmat)),
				as.integer(rarray),
                PACKAGE = "ecodist")
		}

                zstats <- cresults$zstats   

# Calculate the p-values.
		pval1 <- length(zstats[zstats >= zstats[1]])/nperm
		pval2 <- length(zstats[zstats <= zstats[1]])/nperm
		pval3 <- length(zstats[abs(zstats) >= abs(zstats[1])])/nperm

	}		
# If not using a permutation test, return 0 for the p-values.
	else {
		pval1 <- 0
		pval2 <- 0
		pval3 <- 0
	}


	if(nboot > 0) {

		if(dim(m)[[2]] == 2) {
			ycor <- ymat
			xcor <- xmat
		} else {
			xcor <- as.vector(lm.fit(omat, xmat)$residuals)
			ycor <- as.vector(lm.fit(omat, ymat)$residuals)
		}
		bootcor <- numeric(nboot)
		rarray <- numeric(n)
		rmat <- rep(1, length(xcor))
		xdif <- numeric(length(xcor))
		ydif <- numeric(length(xcor))
		cresults <- .C("xbootstrap", 
			as.double(xcor),
			as.double(ycor),
			as.integer(n),
			as.integer(length(xcor)),
			as.integer(nboot),
			as.double(pboot),
			bootcor=as.double(bootcor),
			as.integer(rarray),
			as.integer(rmat),
			as.double(xdif),
			as.double(ydif),
            PACKAGE = "ecodist")

		bootcor <- cresults$bootcor
		bootcor <- sort(bootcor)

		pval <- (1-cboot)/2

		llim <- quantile(bootcor, pval)
		ulim <- quantile(bootcor, 1-pval)
	} else {
		llim <- 0
		ulim <- 0
	}

# Return the Mantel r and the p-value.

	unlist(list(mantelr = mantelr, pval1 = pval1, pval2 = pval2, pval3 = pval3, llim = llim, ulim = ulim))

}


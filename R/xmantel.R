xmantel <- function(formula = formula(data), data = sys.parent(), dims = NA, nperm = 1000, mrank = FALSE)
{
# Cross-mantel test 
# Written by Sarah C. Goslee
# 01/01/01
# Updated 6 April 2001
# added to ecodist 10 July 2017
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
    m <- match.call(expand.dots = FALSE)
    m2 <- match(c("formula", "data"), names(m), nomatch=0)
    m <- m[c(1, m2)]
    m[[1]] <- as.name("model.frame")
    m <- eval(m, parent.frame())
    m <- as.matrix(m)
# m is now the data for the Mantel test as
# cbinded matrices

# Determine the size of the matrices & do some error checking.

    # if dims is NA, matrices are assumed to be square
    # otherwise dims must be specified
    nr <- nrow(m)
    if(is.na(dims[1])) {
        nc <- nr
    } else {
        if(nr != dims[1]) stop("dims doesn't match data \n")
        nc <- dims[2]
    }
    nmat <- ncol(m) / nc

	if(round(nmat) != nmat)
		stop("Matrices don't match dims. \n")

	if(nmat < 2)
		stop("Not enough data. \n")

# If there are only x and y, then use the data as is.

	if(nmat == 2) {
		ymat <- m[,1:nc]
		xmat <- m[,(nc+1):ncol(m)]
		if(mrank) {
			ymat <- rank(ymat)
			xmat <- rank(xmat)
		}
		ycor <- as.vector(ymat)
		xcor <- as.vector(xmat)
	} else {

# If this is a partial Mantel test, get the regression residuals
# for y ~ n1 + n2 + n3 + ... and x ~ n1 + n2 + n3 + ...
		ymat <- as.vector(m[,1:nc])
        omat <- vector(nmat - 1, mode="list")
        for(i in seq(1, nmat - 1)) {
            omat[[i]] <- as.vector(m[, seq(i * nc + 1, (i+1) * nc)])
        }
        omat <- do.call("cbind", omat)

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

		if(nmat > 2) {
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
		rarray <- rep(0, nr)      
		carray <- rep(0, nc)      


		if(nmat == 2) {
             cresults <- .C("xpermute",
                    as.double(xmat),
                    as.double(ymat),
                    as.integer(nr),
                    as.integer(nc),
                    as.integer(length(xmat)),
                    as.integer(nperm),
                    zstats = as.double(zstats),
                    as.double(xmat),
                    as.integer(rarray),
                    as.integer(carray),
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
				as.integer(nr),
				as.integer(nc),
				as.integer(length(xmat)),
				as.integer(nperm),
				zstats = as.double(zstats),
                as.double(as.vector(ymat)),
				as.integer(rarray),
				as.integer(carray),
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



# Return the Mantel r and the p-value.

	unlist(list(mantelr = mantelr, pval1 = pval1, pval2 = pval2, pval3 = pval3))

}


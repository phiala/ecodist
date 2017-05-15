# ClusterLevel returns the Mantel correlations for group contrast
# matrices computed from cluster groups across a range of 
# clustering levels. 
# The inputs are an ecological distance matrix
# and a vector or matrix of cluster levels in any order. 

clusterlevel <- function(edist, groups, nperm=10000, plot=TRUE, ...)  
{
	nl <- ncol(groups)
    if(is.null(nl)) {
        groups <- matrix(groups, ncol=1)
        nl <- 1
    }

	nclust <- rep(NA, nl)
	mantelr <- rep(NA, nl)
	pval <- rep(NA, nl)

	outtable <- data.frame(nclust, mantelr, pval)

    for (i in seq_len(nl)) {
		# number of groups at this level:
		cl <- length(unique(groups[,i]))
		# create group contrast:
		gdist <- dist(groups[,i])
		gdist[gdist!=0] <- 1
		# run mantel:
		m <- mantel(edist ~ gdist, nperm=nperm, nboot=0) 
		outtable[["nclust"]][i] <- cl
		outtable[["mantelr"]][i] <- m[1]
		outtable[["pval"]][i] <- m[2]
    }

	if (plot) {
		plot(outtable[,1], outtable[,2], xlab="Number of Groups", ylab="Mantel r", type="b", ...)
	}   
	return(outtable)
}


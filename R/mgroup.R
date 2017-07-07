# mgroup returns the Mantel correlations for group contrast
# matrices computed from cluster groups across a range of 
# clustering levels. 
# The inputs are an ecological distance matrix
# and a vector or matrix of cluster levels in any order. 
#
# dlu
#

mgroup <- function(edist, groups, nperm=10000)  
{
	nl <- ncol(groups)
    if(is.null(nl)) {
        groups <- matrix(groups, ncol=1)
        nl <- 1
    }

	outtable <- data.frame(nclust = rep(NA, nl), mantelr = rep(NA, nl), pval = rep(NA, nl))

    for (i in seq_len(nl)) {
		# number of groups at this level:
        thisgroups <- groups[,i]
        if(is.factor(thisgroups)) {
            thisgroups <- as.numeric(thisgroups)
        }
        if(is.character(thisgroups)) {
            thisgroups <- as.numeric(factor(thisgroups))
        }

		cl <- length(unique(groups[,i]))
		# create group contrast:
		gdist <- dist(thisgroups)

		gdist[gdist > 0] <- 1

		# run mantel:
		m <- mantel(edist ~ gdist, nperm=nperm, nboot=0) 
		outtable[["nclust"]][i] <- cl
		outtable[["mantelr"]][i] <- m[1]
		outtable[["pval"]][i] <- m[2]
    }

	outtable
}


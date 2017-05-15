plot.nmds <- function(x, plot=TRUE) {
	dims <- seq(x$mindim, x$maxdim, by=1)

	# summarize and plot the nmds object x: stress
	x.stress <- data.frame(matrix(x$stress, nrow=length(dims), byrow=TRUE))
	colnames(x.stress) <- paste0("iter", sprintf("%02d", seq_len(x$nits)))
	x.stress$ndim <- dims
	x.stress$min <- apply(x.stress[, seq_len(x$nits)], 1, min, na.rm=TRUE)
	x.stress$mean <- apply(x.stress[, seq_len(x$nits)], 1, mean, na.rm=TRUE)
	x.stress$max <- apply(x.stress[, seq_len(x$nits)], 1, max, na.rm=TRUE)

	# summarize and plot the nmds object x: r2
	x.r2 <- data.frame(matrix(x$r2, nrow=length(dims), byrow=TRUE))
	colnames(x.r2) <- paste0("iter", sprintf("%02d", seq_len(x$nits)))
	x.r2$ndim <- dims
	x.r2$min <- apply(x.r2[, seq_len(x$nits)], 1, min, na.rm=TRUE)
	x.r2$mean <- apply(x.r2[, seq_len(x$nits)], 1, mean, na.rm=TRUE)
	x.r2$max <- apply(x.r2[, seq_len(x$nits)], 1, max, na.rm=TRUE)

	if(plot) {
		par(mfrow=c(1,2))
		
		plot(dims, x.stress$mean, type="b", pch=19, lwd=2, xlab="Dimensions", ylab="Stress")
		lines(dims, x.stress$min, type="l", lty=2)
		lines(dims, x.stress$max, type="l", lty=2)
				
		plot(dims, x.r2$mean, type="b", pch=19, lwd=2, xlab="Dimensions", ylab="r2")
		lines(dims, x.r2$min, type="l", lty=2)
		lines(dims, x.r2$max, type="l", lty=2)
	}

	invisible(list(stress=x.stress, r2=x.r2))
}


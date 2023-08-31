plot.mgram <- function(x, pval = 0.05, xlab = "Distance", ylab = 
                NULL, ...)
{
# x is the output from mgram
# pval is the p-value to be considered signficant
# ... are additional graphics parameters

        x <- x$mgram

        if(is.null(ylab)) {
            if(colnames(x)[3] == "wtI") ylab <- "Piecewise autocorrelation (I)"
            if(colnames(x)[3] == "mantelr") ylab <- "Mantel autocorrelation (r)"
        }
        pval.v <- x[, 4]
        pval.v[is.na(pval.v)] <- 1
        plot(x[, 1], x[, 3], type = "l", xlab = xlab, ylab = 
                ylab, ...)
        points(x[pval.v <= pval, 1], x[pval.v <= pval, 3], pch = 16, cex=2)
        points(x[pval.v > pval, 1], x[pval.v > pval, 3], pch = 1, cex=2)
        abline(h=0, lty=2, col="gray")
        invisible()
}

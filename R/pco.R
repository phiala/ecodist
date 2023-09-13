pco <- function(x, negvals = "zero", dround=0)

{
# Principal coordinates analysis (classical scaling)
# Sarah Goslee
# x is a lower-triangular dissimilarity matrix.
# negvals: if "zero" then negative eigenvalues are set to zero.
# if "rm" then correction method 1 from Legendre & Anderson 1999
# dround is an attempt to correct for round-off error.

n <- (1 + sqrt(1 + 8 * length(x)))/2
if(abs(n - round(n)) > 0.0000001)
    stop("Matrix not square.\n")
n <- round(n)

xlabels <- labels(x)

x <- full(x)

# note cmdscale does this in two steps
# so dmat does not match cmdscale()$x
dmat <- -0.5 * x * x

# Double-center the dissimilarity matrix
# Subtract row and column means and add the grand mean.

dr <- matrix(apply(dmat, 1, mean), nrow=n, ncol=n)
dmat <- dmat - (dr + t(dr)) + mean(dmat)

deigen <- eigen(dmat)

if(negvals == "rm") {
   c1 <- min(deigen$values)
   if(c1 < 0) {
      c1 <- abs(c1)
      x2 <- sqrt(x^2 + 2 * c1)
      diag(x2) <- 0
      dmat <- -0.5 * x2 * x2
      dr <- matrix(apply(dmat, 1, mean), nrow=n, ncol=n)
      dmat <- dmat - (dr + t(dr)) + mean(dmat)

      deigen <- eigen(dmat)
   }
} else {
   deigen$values[deigen$values < 0] <- 0
}

if(dround > 0) {
   deigen$values <- round(deigen$values, dround)
   deigen$vectors <- round(deigen$vectors, dround)
}

eigenscale <- deigen$values
eigenscale[eigenscale > 0.000000001] <- sqrt(eigenscale[eigenscale > 0.000000001])
deigen$vectors <- sweep(deigen$vectors, 2, eigenscale, "*")

deigen$vectors <- data.frame(deigen$vectors)
rownames(deigen$vectors) <- xlabels
colnames(deigen$vectors) <- paste0("X", seq_len(ncol(deigen$vectors)))

deigen

}

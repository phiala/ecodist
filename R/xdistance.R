xdistance <- function(x, y, method = "euclidean") {
    # calculate cross-dissimilarities between
    # rows of x and rows of y
    # returns a nonsymmetric matrix where
    # d[a, b] is the dissimilarity between
    # x[a, ] and y[b, ]

    # Sarah Goslee 2017-02-17, modified from legacy Splus code dated 01/01/01

    if(is.null(ncol(x))) {
        x <- matrix(x, ncol=1)
        rownames(x) <- seq_len(nrow(x))
    }
    if(is.null(ncol(y))) {
        y <- matrix(y, ncol=1)
        rownames(y) <- seq_len(nrow(y))
    }


    if(!(ncol(x) == ncol(y)))
        stop("Matrices must have the same number of columns\n")

    x.names <- paste0("x", row.names(x))
    y.names <- paste0("y", row.names(y))
    x <- as.matrix(x)
    y <- as.matrix(y)

    d <- rbind(x, y)

    d <- full(distance(d, method=method))


    d <- d[seq(1, nrow(x)), seq(nrow(x) + 1, nrow(x) + nrow(y)), drop=FALSE]
    rownames(d) <- x.names
    colnames(d) <- y.names

    class(d) <- "xdist"
    d
}


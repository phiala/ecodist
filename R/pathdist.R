pathdist <- function(v, maxv = 1) {

	# v is a lower-triangular distance matrix
	# maxv is the cutoff for distances: values greater or equal to this will be estimated 
	# from the shortest path along the distance-weighted graph connecting the samples
    # note that this will not work with completely disconnected subsets

	n <- (1 + sqrt(1 + 8 * length(v)))/2
	if (abs(n - round(n)) > 1e-07)
		stop("Matrix not square.")
	n <- round(n)

	###


    v.ind <- data.frame(v1 = lower(row(full(v))), v2 = lower(col(full(v))))

    v.ind$d <- v
    v.ind <- subset(v.ind, v.ind$d < maxv)

    v.graph <- add_edges(make_empty_graph(n, directed = FALSE), edges = t(v.ind[, 1:2]), weight = v.ind$d)

	v.dist <- distances(v.graph)
	v.dist <- lower(v.dist)

    ## give the result appropriate attributes based on the original object
    attributes(v.dist) <- attributes(v)

    if(!is.null(attr(v.dist, "method"))) {
        attr(v.dist, "method") <- paste(attr(v.dist, "method"), "mst", sep="-")
    } else {
        attr(v.dist, "method") <- "mst"
    }

    v.dist
}


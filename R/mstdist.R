mstdist <- function(v, maxv = 1) {

	# v is a lower-triangular distance matrix
	# maxv is the cutoff for distances: values greater or equal to this will be estimated 
	# from the minimum spanning tree

	n <- (1 + sqrt(1 + 8 * length(v)))/2
	if (abs(n - round(n)) > 1e-07)
		stop("Matrix not square.")
	n <- round(n)

	###


	v.ind <- expand.grid(v1 = seq_len(n), v2 = seq_len(n))

	# lower triangular matrix: row > column
	v.ind <- subset(v.ind, v1 > v2)

	v.graph <- add_edges(make_empty_graph(n), edges = t(v.ind), weight = ifelse(v >= maxv, sum(v) * 10, v))

	v.dist <- distances(v.graph)

	lower(v.dist)

}


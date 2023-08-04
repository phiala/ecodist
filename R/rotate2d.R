rotate2d <- function(ord, x) {

	# rotates a two-dimensional ordination configuration to match a vector x
    # so that x is along the horizontal axis to the right
	# facilitates display

    # get theta from coordinates of x
    # assumes first two values of x are the first two coordinates
    # this is correct for vf() output
    # if ord is a vector, also takes the first two values
    #
    # returns a matrix with ncol = 2 regardless of ord input type

    # extract configuration from pco() ordination output
    # does NOT preserve original, just rotated coordinates
    if(any(names(ord) == "vectors")) {
        ord <- ord$vectors
    }

    # extract configuration from vf() output
    # DOES preserve original, so vf.plot can still be used
    isvf <- FALSE
    if(inherits(ord, "vf")) {
        vf.orig <- ord
        isvf <- TRUE
    }


    # works for a vector of length 2 or a two-dimensional configuration
    if(!is.null(dim(ord))) {
        ord <- as.matrix(ord[, 1:2, drop = FALSE])
    } else {
        # assumes first two values of x are the coordinates, as from vf()
        ord <- matrix(ord[1:2], ncol = 2)
    }

    # assumes first two values of x are the coordinates, as from vf()
    x <- x[1:2]

    ###

    theta <- - atan2(x[2], x[1])

    ord.rot <- data.frame(x = ord[, 1] * cos(theta) - ord[, 2] * sin(theta), 
                      y = ord[, 1] * sin(theta) + ord[, 2] * cos(theta))

    if(isvf) {
        vf.orig[, 1:2] <- as.matrix(ord.rot)
        ord.rot <- vf.orig
    }

    ord.rot

}


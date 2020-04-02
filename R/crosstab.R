crosstab <- function(rowlab, collab, values, type="sum", data, allrows, allcols, na.as.0 = TRUE, check.names=TRUE, ...)

{

# Converts field data in the form:
# site, species, observation
# into a site by species table
#
# By default, takes the sum of the data
# Can also use "mean", "max" or "min"
# Sarah C. Goslee
# 24 Sept 2003
#
###
#
# if allrows or allcols exists, will expand the matrix to
# include all of those elements into the rows or colums of
# the results matrix
# scg 20 Apr 2006
#
# added data argument to simplify using data frames
# added count option to type (length of matching data)
# scg 2 Nov 2012
#
# fixed bug with expanding single-row or single-column values using allrows or allcols
# scg 17 Feb 2017
# added na.as.0 argument; scg 18 Apr 2017

# added check.names and data frame return 21 Jun 2017

    if(!missing(data)) {
        if(mode(substitute(rowlab)) == "name") rowlab <- data[, deparse(substitute(rowlab))]
        if(mode(substitute(collab)) == "name") collab <- data[, deparse(substitute(collab))]
        if(!missing(values) & mode(substitute(values)) == "name") values <- data[, deparse(substitute(values))]
    }

        rowlab <- as.vector(rowlab)
        collab <- as.vector(collab)

        # if values are not provided, a count of combinations of rowlab and collab is returned
        # equivalent to table(rowlab, collab)
        if(missing(values)) {
            values <- rep(1, length(rowlab))
            type <- "sum" 
        }

        # if type is count and values are provided, combinations of unique values are counted
        if(type == "count") {
            values <- paste(rowlab, collab, values)
            values <- as.numeric(!duplicated(values))
            type <- "sum" 
        }

        values <- as.vector(values)

    results <- switch(type,
        mean = tapply(values, list(rowlab, collab), mean, ...),
        max = tapply(values, list(rowlab, collab), max, ...),
        min = tapply(values, list(rowlab, collab), min, ...),
        sum = tapply(values, list(rowlab, collab), sum, ...)
    )

    if(!missing(allrows)) {
        allrows <- as.vector(allrows)
        allrows <- c(rownames(results), allrows)
        allrows <- sort(unique(allrows))
        newrows <- allrows[!(allrows %in% rownames(results))]
        temp <- matrix(NA, ncol=ncol(results), nrow=length(newrows))
        colnames(temp) <- colnames(results)
        rownames(temp) <- newrows
        results <- rbind(results, temp)
        if(is.numeric(rowlab)) {
            results <- results[order(as.numeric(rownames(results))), , drop=FALSE]
        } else {
            results <- results[order(rownames(results)), , drop = FALSE]
        }
    }

    if(!missing(allcols)) {
        allcols <- as.vector(allcols)
        allcols <- c(colnames(results), allcols)
        allcols <- sort(unique(allcols))
        newcols <- allcols[!(allcols %in% colnames(results))]
        temp <- matrix(NA, nrow=nrow(results), ncol=length(newcols))
        rownames(temp) <- rownames(results)
        colnames(temp) <- newcols
        results <- cbind(results, temp)
        if(is.numeric(collab)) {
            results <- results[, order(as.numeric(colnames(results))), drop=FALSE]
        } else {
            results <- results[, order(colnames(results)), drop = FALSE]
        }
    }

    if(na.as.0) results[is.na(results)] <- 0

    data.frame(results, check.names=check.names)

}


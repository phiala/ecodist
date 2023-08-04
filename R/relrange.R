# modified 2017-11-20

relrange <- function(x, globalmin=NA, globalmax=NA) {
    # relativize the range of each column to 0-1
    # if globalmin, globalmax are provided, uses those, eg to
    # scale a subset to match a larger sample
    # otherwise uses population min and max
    if(is.na(globalmin[1])) globalmin <- apply(x, 2, min)
    if(is.na(globalmax[1])) globalmax <- apply(x, 2, max)
    globalmax <- globalmax - globalmin
    x <- sweep(x, 2, globalmin, "-")
    x <- sweep(x, 2, globalmax, "/")
    x
}


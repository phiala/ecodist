\name{crosstab}
\alias{crosstab}
\title{ Data formatting }
\description{
  Converts field data of the form site, species, observation into a site
  by species data frame.
}
\usage{
crosstab(rowlab, collab, values, type = "sum", data, allrows, allcols,
na.as.0 = TRUE, check.names = TRUE, \dots)
}
\arguments{
  \item{rowlab}{ row labels, e.g. site names. }
  \item{collab}{ column labels, e.g. species names. }
  \item{values}{ data values. }
  \item{data}{ optional data frame from which to take rowlab,
    collab and/or values. }
  \item{type}{ function to use to combine data, one of "sum"
    (default), "min", "max", "mean", "count". }
  \item{allrows}{ optional, list of all desired row names that
    may not appear in rowlab. }
  \item{allcols}{ optional, list of all desired column names
    that may not appear in collab. }
  \item{na.as.0}{ if TRUE, all NA values are replaced with 0. }
  \item{check.names}{ if FALSE, data frame names are not checked
    for syntactic validity, so that they match the input categories.
    Otherwise make.names() is used to adjust them. }
  \item{\dots}{ optional arguments to the function specified in type,
    such as na.rm=TRUE}
}
\details{
Field data are often recorded as a separate row for each
site-species combination. This function reformats such data into
a data frame for further analysis based on unique row and column labels.
The three vectors should all be the same length (including duplicates).
The three vectors may also be provided as names of columns in the data
frame specified by the data argument.

If allrows or allcols exists, rows and/or columns of zeros are inserted
for any elements of allrows/allcols not present in rowlab/collab.

If values is missing the number of occurrences of combinations of rowlab
and collab will be returned. Thus, crosstab(rowlab, collab) is equivalent
to table(rowlab, collab).

If type is "count", the unique combinations of rowlab, collab and values
will be returned.
}
\value{
  data frame with rowlab as row headings, collab as columns,
  and values as the data.
}
\author{ Sarah Goslee }
\examples{
# Make a random example
plotnames <- rep(1:5, each = 6)
speciesnames <- rep(c("A", "B", "C"), 10)
freqdata <- runif(30)

# number of samples of each species and plot
crosstab(plotnames, speciesnames)

# can use the data argument
speciesdata <- data.frame(plots = plotnames, species = speciesnames,
  freq = freqdata, stringsAsFactors=FALSE)

# mean frequency by species and plot
crosstab(plots, species, freq, data=speciesdata, type="mean")

# can specify additional possible row or column levels
crosstab(plots, species, freq, data=speciesdata, type="mean", allcols=LETTERS[1:5])
}
\keyword{ manip }
\keyword{ multivariate }

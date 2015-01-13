\name{crosstab}
\alias{crosstab}
\title{ Data formatting }
\description{
  Converts field data of the form site, species, observation into a matrix.
}
\usage{
crosstab(rowlab, collab, values, type = "sum", data, allrows, allcols, ...)
}
\arguments{
  \item{rowlab}{ row labels, e.g. site names. }
  \item{collab}{ column labels, e.g. species names. }
  \item{values}{ data values. }
  \item{data}{ optional data argument to take rowlab, collab and/or values from. }
  \item{type}{ function to use to combine data, one of "sum" (default), "min", "max", "mean", "count". }
  \item{allrows}{ optional, list of all desired row names that may not appear in rowlab. }
  \item{allcols}{ optional, list of all desired column names that may not appear in collab. }
  \item{...}{ optional arguments to the function specified in type, such as na.rm=TRUE}
}

\details{
Field data are often recorded as a separate row for each site-species
combination. This function turns them into a matrix for further analysis based on unique row and column labels.

The three vectors should all be the same length (including duplicates). The three vectors may also be provided as names of columns in the data frame specified by the data argument.

If allrows or allcols exists, rows and/or columns of zeros are inserted for any elements of allrows/allcols not present in rowlab/collab.

If values is missing the number of occurrences of combinations of rowlab and collab will be returned. Thus, crosstab(rowlab, collab) is equivalent to table(rowlab, collab).

If type is "count", the unique combinations of rowlab, collab and values will be returned.

}
\value{
matrix with rowlab as row headings, collab as columns, and
values as the data.
}
\author{ Sarah Goslee, Sarah.Goslee@ars.usda.gov }

\examples{
# Make a random example
plotnames <- sort(rep(1:5, 6))
speciesnames <- rep(c("A", "B", "C"), 10)
freqdata <- runif(30)

# number of samples of each species and plot
crosstab(plotnames, speciesnames)

# mean frequency by species and plot
crosstab(plotnames, speciesnames, freqdata, type="mean")


}
\keyword{ manip }


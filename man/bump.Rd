\name{bump}
\alias{bump}
\docType{data}
\title{Nine-bump spatial pattern}
\description{
  A two-dimensional artificial "landscape" illustrating the kind of spatial pattern that might be seen across mountain peaks.
}
\usage{data(bump)}
\format{
  The format is:
 int [1:25, 1:25] 2 2 2 2 2 2 2 2 2 2 ...
 - attr(*, "dimnames")=List of 2
  ..$ : chr [1:25] "1" "3" "5" "7" ...
  ..$ : chr [1:25] "V1" "V3" "V5" "V7" ...
}
\author{ Sarah Goslee }
\seealso{ \code{\link{bump.pmgram}}, \code{\link{pmgram}} }

\examples{
data(bump)
image(bump)
}
\keyword{datasets}

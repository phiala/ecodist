\name{rotate2d}
\alias{rotate2d}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Rotate a 2D ordination.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
rotate2d(ord, x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ord}{
%%     ~~Describe \code{ord} here~~
}
  \item{x}{
%%     ~~Describe \code{x} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (ord, x) 
{
    if (any(names(ord) == "vectors")) {
        ord <- ord$vectors
    }
    isvf <- FALSE
    if (inherits(ord, "vf")) {
        vf.orig <- ord
        isvf <- TRUE
    }
    if (!is.null(dim(ord))) {
        ord <- as.matrix(ord[, 1:2, drop = FALSE])
    }
    else {
        ord <- matrix(ord[1:2], ncol = 2)
    }
    x <- x[1:2]
    theta <- -atan2(x[2], x[1])
    ord.rot <- data.frame(x = ord[, 1] * cos(theta) - ord[, 2] * 
        sin(theta), y = ord[, 1] * sin(theta) + ord[, 2] * cos(theta))
    if (isvf) {
        vf.orig[, 1:2] <- as.matrix(ord.rot)
        ord.rot <- vf.orig
    }
    ord.rot
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.

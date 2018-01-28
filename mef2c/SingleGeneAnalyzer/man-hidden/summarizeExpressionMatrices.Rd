\name{summarizeExpressionMatrices}
\alias{summarizeExpressionMatrices}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
summarizeExpressionMatrices()
}
%- maybe also 'usage' for other objects documented here.
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
function () 
{
    matrix.count <- length(expression.matrices)
    empty.vec <- rep(0, matrix.count)
    tbl.mtx <- data.frame(name = names(expression.matrices), 
        nrow = empty.vec, ncol = empty.vec, min = empty.vec, 
        q1 = empty.vec, median = empty.vec, q3 = empty.vec, max = empty.vec, 
        stringsAsFactors = FALSE)
    rownames(tbl.mtx) <- tbl.mtx$name
    tbl.mtx <- tbl.mtx[, -1]
    for (name in rownames(tbl.mtx)) {
        mtx <- expression.matrices[[name]]
        summary.stats <- fivenum(mtx)
        dimensions <- dim(mtx)
        tbl.mtx[name, c("min", "q1", "median", "q3", "max")] <- summary.stats
        tbl.mtx[name, c("nrow", "ncol")] <- dimensions
    }
    dataFrameToPandasFriendlyList(tbl.mtx)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

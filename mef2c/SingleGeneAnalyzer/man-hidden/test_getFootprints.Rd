\name{test_getFootprints}
\alias{test_getFootprints}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
test_getFootprints()
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
    printf("--- test_getFootprints")
    test.roi <- fp.roi
    test.roi$start <- test.roi$start + 10000
    test.roi$end <- test.roi$start + 20000
    test.roi.string <- with(test.roi, sprintf("\%s:\%d-\%d", chrom, 
        start, end))
    fp.roi.string <- with(fp.roi, sprintf("\%s:\%d-\%d", chrom, 
        start, end))
    tbl.test <- getFootprints(test.roi.string)
    checkEquals(nrow(tbl.test), 203)
    checkEquals(ncol(tbl.test), 12)
    checkEquals(colnames(tbl.test)[1:3], c("chrom", "start", 
        "end"))
    outOfRange.roi <- sprintf("\%s:\%d-\%d", fp.roi$chrom, fp.roi$end + 
        10, fp.roi$end + 1000)
    tbl.outOfRange <- getFootprints(outOfRange.roi)
    checkTrue(nrow(tbl.outOfRange) > 20)
    checkEquals(ncol(tbl.outOfRange), 12)
    checkEquals(colnames(tbl.outOfRange)[1:3], c("chrom", "start", 
        "end"))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

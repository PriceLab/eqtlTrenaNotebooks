\name{test_getDHSMotifs}
\alias{test_getDHSMotifs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
test_getDHSMotifs()
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
    printf("--- test_getDHSMotifs")
    test.roi <- fp.roi
    test.roi$start <- test.roi$start + 10000
    test.roi$end <- test.roi$start + 20000
    test.roi.string <- with(test.roi, sprintf("\%s:\%d-\%d", chrom, 
        start, end))
    fp.roi.string <- with(fp.roi, sprintf("\%s:\%d-\%d", chrom, 
        start, end))
    tbl.dhs.new <- getDHSMotifs(fp.roi.string)
    checkEquals(dim(tbl.dhsMotifs), dim(tbl.dhs.new))
    checkEquals(colnames(tbl.dhs.new)[1:3], c("chrom", "start", 
        "end"))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

\name{test_assay}
\alias{test_assay}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
test_assay(tbl.snp.tf, tf, tbl.geno, mtx)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tbl.snp.tf}{
%%     ~~Describe \code{tbl.snp.tf} here~~
}
  \item{tf}{
%%     ~~Describe \code{tf} here~~
}
  \item{tbl.geno}{
%%     ~~Describe \code{tbl.geno} here~~
}
  \item{mtx}{
%%     ~~Describe \code{mtx} here~~
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
function (tbl.snp.tf, tf, tbl.geno, mtx) 
{
    assay(tbl.20, "EMX1", tbl.geno, tbl.pheno, mtx.tcx.normalized)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

\name{getFootprints}
\alias{getFootprints}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getFootprints(roiString)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{roiString}{
%%     ~~Describe \code{roiString} here~~
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
function (roiString) 
{
    roi <- parseChromLocString(roiString)
    if (roi$chrom == fp.roi$chrom & roi$start >= fp.roi$start & 
        roi$end <= fp.roi$end) {
        tbl.roi <- subset(tbl.fp, chrom == roi$chrom & start >= 
            roi$start & end <= roi$end)
        return(tbl.roi)
    }
    trena <- Trena("hg38")
    source.1 <- "postgres://bddsrds.globusgenomics.org/brain_wellington_16"
    source.2 <- "postgres://bddsrds.globusgenomics.org/brain_wellington_20"
    source.3 <- "postgres://bddsrds.globusgenomics.org/brain_hint_16"
    source.4 <- "postgres://bddsrds.globusgenomics.org/brain_hint_20"
    sources <- c(source.1, source.2, source.3, source.4)
    names(sources) <- c("well_16", "well_20", "hint_16", "hint_20")
    x <- getRegulatoryChromosomalRegions(trena, roi$chrom, roi$start, 
        roi$end, sources, targetGene, targetGene.tss)
    print(1)
    names(x) <- names(sources)
    print(2)
    x2 <- lapply(names(x), function(name) {
        tbl <- x[[name]]
        if (nrow(tbl) > 0) 
            tbl$db <- name
        return(tbl)
    })
    print(3)
    tbl.reg <- do.call(rbind, x2)
    print(4)
    rownames(tbl.reg) <- NULL
    print(5)
    tbl.reg <- unique(tbl.reg[grep("Hsapiens-jaspar2016", tbl.reg$motifName, 
        ignore.case = TRUE), ])
    tbl.reg <- tbl.reg[order(tbl.reg$motifStart), ]
    printf("---- getFootprints, before associateTranscriptionFactors")
    print(head(tbl.reg))
    tbl.reg <- associateTranscriptionFactors(MotifDb, tbl.reg, 
        source = "MotifDb", expand.rows = TRUE)
    colnames(tbl.reg)[2:3] <- c("start", "end")
    print(6)
    tbl.reg
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

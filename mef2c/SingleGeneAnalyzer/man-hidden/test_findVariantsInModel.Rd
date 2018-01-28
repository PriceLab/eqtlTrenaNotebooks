\name{test_findVariantsInModel}
\alias{test_findVariantsInModel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
test_findVariantsInModel()
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
    printf("--- test_findVariantsInModel")
    tbl.00.fp <- findVariantsInModel(modelName = "mef2c.tcx", 
        shoulder = 0, motifSource = c("footprints"))
    checkEquals(dim(tbl.00.fp), c(3, 24))
    tbl.00.dhs <- findVariantsInModel(modelName = "mef2c.tcx", 
        shoulder = 0, motifSource = c("dhs"))
    checkEquals(dim(tbl.00.dhs), c(2, 24))
    tbl.00.allDNA <- findVariantsInModel(modelName = "mef2c.tcx", 
        shoulder = 0, motifSource = c("allDNA"))
    checkEquals(dim(tbl.00.allDNA), c(14, 24))
    tbl.15 <- findVariantsInModel(modelName = "mef2c.tcx", shoulder = 15, 
        motifSource = c("footprints", "dhs"))
    checkEquals(dim(tbl.15), c(15, 24))
    tbl.xtab <- as.data.frame(table(tbl.00.fp$source), stringsAsFactors = FALSE)
    checkEquals(tbl.xtab$Var1, c("hint+wellington, 16+20 footprints"))
    checkEquals(tbl.xtab$Freq, c(3))
    tbl.xtab <- as.data.frame(table(tbl.00.fp$geneSymbol), stringsAsFactors = FALSE)
    checkEquals(tbl.xtab$Var1, c("FOXP1", "MEF2A"))
    checkEquals(tbl.xtab$Freq, c(2, 1))
    coi <- c(1, 2, 3, 4, 6, 7, 17, 19, 23, 24, 20, 21, 22)
    tbl.05 <- findVariantsInModel(modelName = "mef2c.tcx", motifSources = c("footprints", 
        "dhs", "allDNA"), shoulder = 5)
    checkEquals(dim(tbl.05), c(34, 24))
    tbl.xtab <- as.data.frame(table(tbl.05$source), stringsAsFactors = FALSE)
    checkEquals(tbl.xtab$Var1, c("all DNA motifs", "encode DHS motifs", 
        "hint+wellington, 16+20 footprints"))
    checkEquals(tbl.xtab$Freq, c(25, 3, 6))
    tbl.xtab <- as.data.frame(table(tbl.05$geneSymbol), stringsAsFactors = FALSE)
    checkEquals(tbl.xtab$Var1, c("BHLHE22", "EMX1", "FOXP1", 
        "FOXP2", "HLF", "LBX2", "MEF2A", "NFE2L2", "SP3"))
    checkEquals(tbl.xtab$Freq, c(1, 10, 2, 3, 2, 10, 3, 2, 1))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

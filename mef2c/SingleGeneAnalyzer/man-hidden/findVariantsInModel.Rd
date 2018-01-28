\name{findVariantsInModel}
\alias{findVariantsInModel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
findVariantsInModel(modelName, shoulder, motifSources)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{modelName}{
%%     ~~Describe \code{modelName} here~~
}
  \item{shoulder}{
%%     ~~Describe \code{shoulder} here~~
}
  \item{motifSources}{
%%     ~~Describe \code{motifSources} here~~
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
function (modelName, shoulder, motifSources) 
{
    stopifnot(all(motifSources \%in\% c("footprints", "dhs", "allDNA")))
    tbl.snp.clean <- tbl.snp[, c("chrom", "pos", "pos")]
    colnames(tbl.snp.clean) <- c("chrom", "start", "end")
    gr.snp <- GRanges(tbl.snp.clean)
    tbl.fpHits <- data.frame()
    tbl.dhsMotifHits <- data.frame()
    tbl.allMotifHits <- data.frame()
    if ("footprints" \%in\% motifSources) {
        gr.fp <- GRanges(tbl.fp)
        tbl.fp.ov <- as.data.frame(findOverlaps(gr.snp, gr.fp, 
            maxgap = shoulder))
        if (nrow(tbl.fp.ov) > 0) {
            colnames(tbl.fp.ov) <- c("snp", "fp")
            tbl.fpHits <- cbind(tbl.snp[tbl.fp.ov$snp, ], tbl.fp[tbl.fp.ov$fp, 
                "geneSymbol", drop = FALSE])
            tbl.fpHits <- unique(subset(tbl.fpHits, geneSymbol \%in\% 
                gene.models[[modelName]]$gene))
            if (nrow(tbl.fpHits) > 0) {
                tbl.fpHits$modelName <- modelName
                tbl.fpHits$shoulder <- shoulder
                tbl.fpHits$source <- "hint+wellington, 16+20 footprints"
            }
        }
    }
    if ("dhs" \%in\% motifSources) {
        gr.dhsMotifs <- GRanges(tbl.dhsMotifs)
        tbl.dhs.ov <- as.data.frame(findOverlaps(gr.snp, gr.dhsMotifs, 
            maxgap = shoulder))
        if (nrow(tbl.dhs.ov) > 0) {
            colnames(tbl.dhs.ov) <- c("snp", "dhsMotif")
            tbl.dhsMotifHits <- cbind(tbl.snp[tbl.dhs.ov$snp, 
                ], tbl.dhsMotifs[tbl.dhs.ov$dhsMotif, c("geneSymbol"), 
                drop = FALSE])
            tbl.dhsMotifHits <- unique(subset(tbl.dhsMotifHits, 
                geneSymbol \%in\% gene.models[[modelName]]$gene))
            if (nrow(tbl.dhsMotifHits) > 0) {
                tbl.dhsMotifHits$modelName <- modelName
                tbl.dhsMotifHits$shoulder <- shoulder
                tbl.dhsMotifHits$source <- "encode DHS motifs"
            }
        }
    }
    if ("allDNA" \%in\% motifSources) {
        gr.allDNAMotifs <- GRanges(tbl.allDNAMotifs)
        tbl.all.ov <- as.data.frame(findOverlaps(gr.snp, gr.allDNAMotifs, 
            maxgap = shoulder))
        if (nrow(tbl.all.ov) > 0) {
            colnames(tbl.all.ov) <- c("snp", "allMotif")
            tbl.allMotifHits <- cbind(tbl.snp[tbl.all.ov$snp, 
                ], tbl.allDNAMotifs[tbl.all.ov$allMotif, c("geneSymbol"), 
                drop = FALSE])
            tbl.allMotifHits <- unique(subset(tbl.allMotifHits, 
                geneSymbol \%in\% gene.models[[modelName]]$gene))
            if (nrow(tbl.allMotifHits) > 0) {
                tbl.allMotifHits$modelName <- modelName
                tbl.allMotifHits$shoulder <- shoulder
                tbl.allMotifHits$source <- "all DNA motifs"
            }
        }
    }
    tbl.out <- rbind(tbl.fpHits, tbl.dhsMotifHits)
    tbl.out <- rbind(tbl.out, tbl.allMotifHits)
    genes.in.model <- gene.models[[modelName]]$gene
    motif.tfs <- tbl.out$geneSymbol
    match.order <- match(tbl.out$geneSymbol, gene.models[[modelName]]$gene)
    tbl.out$tfRank <- match.order
    tbl.out$pearson <- gene.models[[modelName]]$pearson.coef[match.order]
    rownames(tbl.out) <- NULL
    tbl.out
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

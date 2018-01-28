\name{processWellStructuredMessage}
\alias{processWellStructuredMessage}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
processWellStructuredMessage(msg)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{msg}{
%%     ~~Describe \code{msg} here~~
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
function (msg) 
{
    if (msg$cmd == "summarizeExpressionMatrices") {
        tbl.summary.as.list <- summarizeExpressionMatrices()
        response <- list(cmd = msg$callback, status = "success", 
            callback = "", payload = tbl.summary.as.list)
    }
    else if (msg$cmd == "getSessionInfo") {
        info <- as.character(sessionInfo())
        response <- list(cmd = msg$callback, status = "success", 
            callback = "", payload = info)
    }
    else if (msg$cmd == "getExpressionMatrixNames") {
        info <- sort(names(expression.matrix.files))
        response <- list(cmd = msg$callback, status = "success", 
            callback = "", payload = info)
    }
    else if (msg$cmd == "getModelNames") {
        info <- sort(mef2c.model.names)
        response <- list(cmd = msg$callback, status = "success", 
            callback = "", payload = info)
    }
    else if (msg$cmd == "getModel") {
        modelName <- msg$payload
        key <- as.character(as.numeric(Sys.time()) * 1e+05)
        tbl <- gene.models[[modelName]]
        cache[[key]] <- tbl
        tbl.fp.as.list <- dataFrameToPandasFriendlyList(tbl)
        payload <- list(tbl = tbl.fp.as.list, key = key)
        response <- list(cmd = msg$callback, status = "success", 
            callback = "", payload = payload)
    }
    else if (msg$cmd == "getVariants") {
        thresholdScore <- msg$payload$minScore
        tbl.snp$filteringScore <- -log10(tbl.snp$CER_P)
        tbl.var <- subset(tbl.snp, filteringScore >= thresholdScore)[, 
            c("chrom", "pos", "pos", "rsid", "filteringScore")]
        colnames(tbl.var) <- c("chrom", "start", "end", "id", 
            "score")
        rownames(tbl.var) <- NULL
        key <- as.character(as.numeric(Sys.time()) * 1e+05)
        tbl.fp.as.list <- dataFrameToPandasFriendlyList(tbl.var)
        payload <- list(tbl = tbl.fp.as.list, key = key)
        response <- list(cmd = msg$callback, status = "success", 
            callback = "", payload = payload)
    }
    else if (msg$cmd == "getFootprintsInRegion") {
        stopifnot("roi" \%in\% names(msg$payload))
        roi <- msg$payload$roi
        tbl.reg <- getFootprints(roi)
        key <- as.character(as.numeric(Sys.time()) * 1e+05)
        cache[[key]] <- tbl.reg
        tbl.fp.as.list <- dataFrameToPandasFriendlyList(tbl.reg)
        payload <- list(tbl = tbl.fp.as.list, key = key)
        response <- list(cmd = msg$callback, status = "success", 
            callback = "", payload = payload)
    }
    else if (msg$cmd == "getDHSRegionsInRegion") {
        stopifnot("roi" \%in\% names(msg$payload))
        roi = msg$payload$roi
        tbl.dhs.roi <- getDHSRegions(roi)
        if (nrow(tbl.dhs.roi) == 0) {
            response <- list(cmd = msg$callback, status = "failure", 
                callback = "", payload = "no overlapping regions")
        }
        else {
            key <- as.character(as.numeric(Sys.time()) * 1e+05)
            cache[[key]] <- tbl.dhs.roi
            tbl.fp.as.list <- dataFrameToPandasFriendlyList(tbl.dhs.roi)
            payload <- list(tbl = tbl.fp.as.list, key = key)
            response <- list(cmd = msg$callback, status = "success", 
                callback = "", payload = payload)
        }
    }
    else if (msg$cmd == "getDHSMotifsInRegion") {
        stopifnot("roi" \%in\% names(msg$payload))
        roi <- msg$payload$roi
        tbl.reg <- getDHSMotifs(roi)
        if (nrow(tbl.reg) == 0) {
            response <- list(cmd = msg$callback, status = "failure", 
                callback = "", payload = "no motifs in region")
        }
        else {
            key <- as.character(as.numeric(Sys.time()) * 1e+05)
            cache[[key]] <- tbl.reg
            tbl.fp.as.list <- dataFrameToPandasFriendlyList(tbl.reg)
            payload <- list(tbl = tbl.fp.as.list, key = key)
            response <- list(cmd = msg$callback, status = "success", 
                callback = "", payload = payload)
        }
    }
    else if (msg$cmd == "findVariantsInModel") {
        stopifnot("modelName" \%in\% names(msg$payload))
        stopifnot("shoulder" \%in\% names(msg$payload))
        modelName <- msg$payload$modelName
        shoulder <- msg$payload$shoulder
        tbl.var <- findVariantsInModel(modelName, shoulder, "footprints")
        key <- as.character(as.numeric(Sys.time()) * 1e+05)
        cache[[key]] <- tbl.var
        tbl.fp.as.list <- dataFrameToPandasFriendlyList(tbl.var)
        payload <- list(tbl = tbl.fp.as.list, key = key)
        response <- list(cmd = msg$callback, status = "success", 
            callback = "", payload = payload)
    }
    else if (msg$cmd == "listSharedData") {
        filenames <- list.files("/home/trena/sharedData")
        response <- list(cmd = msg$callback, status = "success", 
            callback = "", payload = filenames)
    }
    else if (msg$cmd == "createTaggedDataFrame") {
        rows <- msg$payload$rows
        cols <- msg$payload$cols
        tbl <- as.data.frame(matrix(runif(32), nrow = rows, ncol = cols, 
            dimnames = list(LETTERS[1:rows], letters[1:cols])))
        key <- as.character(as.numeric(Sys.time()) * 1e+05)
        cache[[key]] <- tbl
        tbl.for.pandas <- dataFrameToPandasFriendlyList(tbl)
        payload <- list(tbl = tbl.for.pandas, key = key)
        response <- list(cmd = msg$callback, status = "success", 
            callback = "", payload = payload)
    }
    else if (msg$cmd == "identifyTaggedDataFrame") {
        key <- msg$payload
        found <- key \%in\% names(cache)
        tbl.retrieved <- cache[[key]]
        printf("retrieved tbl with dimensions \%d, \%d", nrow(tbl.retrieved), 
            ncol(tbl.retrieved))
        print(tbl.retrieved)
        sum <- sum(as.matrix(tbl.retrieved))
        payload <- list(found = found, sum = sum)
        response <- list(cmd = msg$callback, status = "success", 
            callback = "", payload = payload)
    }
    else if (msg$cmd == "createGeneModel") {
        trena <- Trena("hg38")
        payload <- msg$payload
        solver.names <- payload$solverNames
        printf("---- createGeneModel about to extract tbl.motifs from cache")
        key <- payload$tblRegulatoryRegionsCacheKey
        printf("    key: \%s", key)
        found <- key \%in\% names(cache)
        printf("    found? \%s", found)
        tbl.motifs <- cache[[key]]
        printf("    tbl.motifs: \%d, \%d", nrow(tbl.motifs), ncol(tbl.motifs))
        tfMap <- payload$tfMap
        stopifnot(all(c("targetGene", "matrixName") \%in\% names(msg$payload)))
        targetGene <- toupper(msg$payload$targetGene)
        matrixName <- msg$payload$matrixName
        tbl.motifs <- tbl.motifs[, 1:5]
        colnames(tbl.motifs) <- c("chrom", "start", "end", "motifName")
        print(head(tbl.motifs))
        tbl.motifs$motifName <- sub("_well_16", "", tbl.motifs$motifName)
        tbl.motifs$motifName <- sub("_well_20", "", tbl.motifs$motifName)
        tbl.motifs$motifName <- sub("_hint_16", "", tbl.motifs$motifName)
        tbl.motifs$motifName <- sub("_hint_20", "", tbl.motifs$motifName)
        print(head(tbl.motifs))
        tbl.motifs.tfs <- associateTranscriptionFactors(MotifDb, 
            tbl.motifs, source = tfMap, expand.rows = FALSE)
        mtx <- expression.matrices[[matrixName]]
        stopifnot(targetGene \%in\% rownames(mtx))
        print(dim(mtx))
        print(mean(mtx[targetGene, ]))
        tbl.geneModel <- createGeneModel(trena, targetGene, solver.names, 
            tbl.motifs.tfs, mtx)
        tbl.geneModel <- tbl.geneModel[order(tbl.geneModel$rfScore, 
            decreasing = TRUE), ]
        key <- as.character(as.numeric(Sys.time()) * 1e+05)
        cache[[key]] <- tbl.geneModel
        print(tbl.geneModel)
        payload = list(tbl = dataFrameToPandasFriendlyList(tbl.geneModel), 
            key = key)
        response <- list(cmd = msg$callback, status = "success", 
            callback = "", payload = payload)
    }
    else if (msg$cmd == "buildMultiModelGraph") {
        printf("----- buildMultiModelGraph")
        payload = msg$payload
        print(payload)
        targetGene <- payload$targetGene
        models <- payload$models
        models.decached <- list()
        for (name in names(models)) {
            model <- models[[name]]
            tbl.geneModel <- cache[[model$model]]
            tbl.regions <- cache[[model$regions]]
            tbl.regions <- subset(tbl.regions, geneSymbol \%in\% 
                tbl.geneModel$gene)
            models.decached[[name]]$model <- tbl.geneModel
            models.decached[[name]]$regions <- tbl.regions
        }
        printf("---- after decaching")
        print(models.decached)
        g <- trenaViz::buildMultiModelGraph(targetGene, models.decached)
        xCoordinate.span <- 1500
        g.lo <- trenaViz::addGeneModelLayout(g, xPos.span = xCoordinate.span)
        printf("---- g.lo")
        print(g.lo)
        g.json <- trenaViz:::.graphToJSON(g.lo)
        response <- list(cmd = msg$callback, status = "success", 
            callback = "", payload = g.json)
    }
    else {
        response <- list(cmd = msg$callback, status = "success", 
            callback = "", payload = "well-structured (unparsed) message")
    }
    response
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

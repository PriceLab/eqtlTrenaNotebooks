.SingleGeneAnalyzer <- setClass ("SingleGeneAnalyzer",
                            representation = representation(
                               genomeName="character",
                               targetGene="character",
                               targetGene.TSS="numeric",
                               singleGeneData="SingleGeneData",
                               trena="Trena",
                               quiet="logical"
                               )
                            )

PORT <- 5548
#------------------------------------------------------------------------------------------------------------------------
setGeneric('summarizeExpressionMatrices', signature='obj', function(obj) standardGeneric ('summarizeExpressionMatrices'))
setGeneric('getFootprintsForRegion', signature='obj', function(obj, roiString, score.threshold=NA)
              standardGeneric ('getFootprintsForRegion'))
setGeneric('getVariantsForRegion', signature='obj', function(obj, roiString, score.threshold=NA)
              standardGeneric ('getVariantsForRegion'))
setGeneric('getDHSForRegion', signature='obj', function(obj, roiString, score.threshold=NA) standardGeneric ('getDHSForRegion'))
setGeneric('getEnhancersForRegion',
           signature='obj', function(obj, roiString, score.threshold=NA) standardGeneric ('getEnhancersForRegion'))
#------------------------------------------------------------------------------------------------------------------------
SingleGeneAnalyzer = function(genomeName, targetGene, targetGene.TSS, singleGeneData, quiet=TRUE)
{
   trena <- Trena(genomeName)

   obj <- .SingleGeneAnalyzer(genomeName=genomeName,
                              targetGene=targetGene,
                              targetGene.TSS=targetGene.TSS,
                              singleGeneData=singleGeneData,
                              trena=trena,
                              quiet=quiet)
   obj

} # constructor
#----------------------------------------------------------------------------------------------------
setMethod('summarizeExpressionMatrices', 'SingleGeneAnalyzer',

    function(obj){
       matrix.list <- getExpressionMatrices(obj@singleGeneData)
       matrix.count <- length(matrix.list)
       matrix.names <- names(matrix.list)
       empty.vec <- rep(0, matrix.count)
       tbl.mtx <- data.frame(name=matrix.names,
                             nrow=empty.vec,
                             ncol=empty.vec,
                             min=empty.vec,
                             q1=empty.vec,
                             median=empty.vec,
                             q3=empty.vec,
                             max=empty.vec,
                             stringsAsFactors=FALSE);
       rownames(tbl.mtx) <- tbl.mtx$name
       tbl.mtx <- tbl.mtx[, -1]

       for(name in rownames(tbl.mtx)){
          mtx <- matrix.list[[name]]
          summary.stats <- fivenum(mtx)
          dimensions <- dim(mtx)
          tbl.mtx[name, c("min", "q1", "median", "q3", "max")] <- summary.stats
          tbl.mtx[name, c("nrow", "ncol")] <- dimensions
          } # for mtx.name

       tbl.mtx
       })

#----------------------------------------------------------------------------------------------------
setMethod('getFootprintsForRegion', 'SingleGeneAnalyzer',

    function(obj, roiString, score.threshold=NA){
       roi <- trena::parseChromLocString(roiString)
       tbl.fp <- getFootprints(obj@singleGeneData, roi)

       if(!is.na(score.threshold))
         tbl.fp <- subset(tbl.fp, score >= score.threshold)

       invisible(tbl.fp)
       })

#----------------------------------------------------------------------------------------------------
setMethod('getVariantsForRegion', 'SingleGeneAnalyzer',

    function(obj, roiString, score.threshold=NA){
       roi <- trena::parseChromLocString(roiString)
       tbl.all <- obj@singleGeneData@misc.data[["eqtl.snps"]]
       tbl.sub <- subset(tbl.all, chrom==roi$chrom & pos >= roi$start & pos <= roi$end)
       if(!is.na(score.threshold) & nrow(tbl.sub) > 0)
          tbl.sub <- subset(tbl.sub, -log10(CER_P) >= score.threshold)
       return(tbl.sub)
       })

#----------------------------------------------------------------------------------------------------
setMethod('getDHSForRegion', 'SingleGeneAnalyzer',

    function(obj, roiString, score.threshold=NA){
       roi <- trena::parseChromLocString(roiString)
       tbl.all <- obj@singleGeneData@misc.data[["tbl.dhs"]]
       tbl.sub <- subset(tbl.all, chrom==roi$chrom & start >= roi$start & end <= roi$end)
       if(!is.na(score.threshold))
          tbl.sub <- subset(tbl.sub,  score >= score.threshold)
       return(tbl.sub)
       })

#----------------------------------------------------------------------------------------------------
setMethod('getEnhancersForRegion', 'SingleGeneAnalyzer',

    function(obj, roiString, score.threshold=NA){
       roi <- trena::parseChromLocString(roiString)
       tbl.all <- obj@singleGeneData@misc.data[["enhancer.locs"]]
       tbl.sub <- subset(tbl.all, chrom==roi$chrom & start >= roi$start & end <= roi$end)
       return(tbl.sub)
       })

#----------------------------------------------------------------------------------------------------

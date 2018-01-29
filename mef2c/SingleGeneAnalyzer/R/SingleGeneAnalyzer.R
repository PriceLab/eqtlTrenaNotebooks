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

         # this result is destined for JSON and a python pandas dataframe
         # structure the data.frame as a 3-part list for easy uptake with pandas

       dataFrameToPandasFriendlyList(tbl.mtx)
       })

#----------------------------------------------------------------------------------------------------

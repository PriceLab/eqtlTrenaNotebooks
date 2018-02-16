.SingleGeneData <- setClass ("SingleGeneData",
                             representation = representation(
                                  # chrom:start-end -> region for which footprints were precalculated
                                  # these bounds should expand to inclusively include regions of all data
                                  # variants, enhancers, dhs regions, ...
                                chrom="character",
                                start="numeric",
                                end="numeric",
                                tbl.fp="data.frame",
                                expression.matrices="list",
                                models="list",
                                misc.data="environment",    # for anything not completely standardized
                                quiet="logical"
                                )
                             )

#------------------------------------------------------------------------------------------------------------------------
setGeneric('getGenomicBounds', signature='obj', function(obj, asString=FALSE) standardGeneric('getGenomicBounds'))
setGeneric('getExpressionMatrices', signature='obj', function(obj) standardGeneric ('getExpressionMatrices'))
setGeneric('getFootprints', signature='obj', function(obj, roi) standardGeneric('getFootprints'))
setGeneric('getEnhancers', signature='obj', function(obj, roi) standardGeneric('getEnhancers'))
setGeneric('getModels', signature='obj', function(obj) standardGeneric('getModels'))
setGeneric('getVariants', signature='obj', function(obj, source.name, roi, score.1.threshold=NA_real_,
                                    score.2.threshold=NA_real, score.3.threshold=NA_real_) standardGeneric ('getVariants'))
setGeneric('getMotifs', signature='obj', function(obj, source.name, roi, score.threshold=NA_real_)
                                    standardGeneric ('getMotifs'))


#------------------------------------------------------------------------------------------------------------------------
SingleGeneData <- function(chrom, start, end, tbl.fp, expression.matrices, models, misc.data)
{
   .SingleGeneData(chrom=chrom, start=start, end=end,
                   tbl.fp=tbl.fp,
                   expression.matrices=expression.matrices,
                   models=models,
                   misc.data=misc.data)

} # ctor
#------------------------------------------------------------------------------------------------------------------------
setMethod('getGenomicBounds', 'SingleGeneData',

    function(obj, asString=FALSE){
       result <- list(chrom=obj@chrom, start=obj@start, end=obj@end)
       if(asString)
          result <- with(result, sprintf("%s:%d-%d", chrom, start, end))
       result
       })

#----------------------------------------------------------------------------------------------------
setMethod('getExpressionMatrices', 'SingleGeneData',

    function(obj){
       invisible(obj@expression.matrices)
       })

#----------------------------------------------------------------------------------------------------
setMethod('getFootprints', 'SingleGeneData',

    function(obj, roi=NA){
       if(all(is.na(roi)))
          invisible(obj@tbl.fp)
       tbl.sub <- with(roi, subset(tbl.fp, chrom==chrom & start >= start & end <= end))
       invisible(tbl.sub)
       })

#----------------------------------------------------------------------------------------------------
setMethod('getModels', 'SingleGeneData',

    function(obj){
       obj@models
       })

#----------------------------------------------------------------------------------------------------

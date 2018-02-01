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
setGeneric('getFootprintsForRegion', signature='obj', function(obj, roi.string, score.threshold=NA)
              standardGeneric ('getFootprintsForRegion'))
setGeneric('getVariantsForRegion', signature='obj', function(obj, roi.string, score.threshold=NA)
              standardGeneric ('getVariantsForRegion'))
setGeneric('getDHSForRegion', signature='obj', function(obj, roi.string, score.threshold=NA) standardGeneric ('getDHSForRegion'))
setGeneric('getEnhancersForRegion',
           signature='obj', function(obj, roi.string, score.threshold=NA) standardGeneric ('getEnhancersForRegion'))
setGeneric('findVariantsInModelForRegion', signature='obj',
           function(obj, roi.string, model.name, shoulder, tf.count=NA) standardGeneric ('findVariantsInModelForRegion'))
setGeneric('findMotifsInRegion', signature='obj',
           function(obj, roi.string, motifs, pwmMatchPercentage, variants=NA_character)
              standardGeneric ('findMotifsInRegion'))

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

    function(obj, roi.string, score.threshold=NA){
       roi <- trena::parseChromLocString(roi.string)
       tbl.fp <- getFootprints(obj@singleGeneData, roi)

       if(!is.na(score.threshold))
         tbl.fp <- subset(tbl.fp, score >= score.threshold)

       invisible(tbl.fp)
       })

#----------------------------------------------------------------------------------------------------
setMethod('getVariantsForRegion', 'SingleGeneAnalyzer',

    function(obj, roi.string, score.threshold=NA){
       roi <- trena::parseChromLocString(roi.string)
       tbl.all <- obj@singleGeneData@misc.data[["eqtl.snps"]]
       tbl.sub <- subset(tbl.all, chrom==roi$chrom & pos >= roi$start & pos <= roi$end)
       if(!is.na(score.threshold) & nrow(tbl.sub) > 0)
          tbl.sub <- subset(tbl.sub, -log10(CER_P) >= score.threshold)
       return(tbl.sub)
       })

#----------------------------------------------------------------------------------------------------
setMethod('getDHSForRegion', 'SingleGeneAnalyzer',

    function(obj, roi.string, score.threshold=NA){
       roi <- trena::parseChromLocString(roi.string)
       tbl.all <- obj@singleGeneData@misc.data[["tbl.dhs"]]
       tbl.sub <- subset(tbl.all, chrom==roi$chrom & start >= roi$start & end <= roi$end)
       if(!is.na(score.threshold))
          tbl.sub <- subset(tbl.sub,  score >= score.threshold)
       return(tbl.sub)
       })

#----------------------------------------------------------------------------------------------------
setMethod('getEnhancersForRegion', 'SingleGeneAnalyzer',

    function(obj, roi.string, score.threshold=NA){
       roi <- trena::parseChromLocString(roi.string)
       tbl.all <- obj@singleGeneData@misc.data[["enhancer.locs"]]
       tbl.sub <- subset(tbl.all, chrom==roi$chrom & start >= roi$start & end <= roi$end)
       return(tbl.sub)
       })

#----------------------------------------------------------------------------------------------------
setMethod('findMotifsInRegion', 'SingleGeneAnalyzer',

      function(obj, roi.string, motifs, pwmMatchPercentage, variants=NA_character){
          roi <- parseChromLocString(roi.string)
          tbl.regions <- with(roi, data.frame(chrom=chrom, start=start, end=end, stringsAsFactors=FALSE))
          findMotif <- function(motif){
             mm <- MotifMatcher("hg38", as.list(query(query(MotifDb, motif), "jaspar2018")))
             tbl <- findMatchesByChromosomalRegion(mm, tbl.regions, pwmMatchPercentage)
             if(nrow(tbl) > 0)
                tbl$motif <- motif
             tbl
             } # findMotif
          tbls.all <- lapply(motifs, findMotif)
          tbl.all <- do.call(rbind, tbls.all)
          if(nrow(tbl.all) == 0)
             return(data.frame())
          tbl.all <- tbl.all[order(tbl.all$motifStart, tbl.all$motifEnd, decreasing=FALSE),]
             # put the incomping "motif" in the firest column
          preferred.column.order <- c("motif", "chrom","motifStart","motifEnd","strand","motifName","motifScore",
                                      "motifRelativeScore","match","chromStart","chromEnd","seq","status")
          tbl.all[, preferred.column.order]
          }) # findMotifsInR

#----------------------------------------------------------------------------------------------------
# strategy:
#   identify the tfs in the mode
#   get the motifs associated with each
#
setMethod('findVariantsInModelForRegion', 'SingleGeneAnalyzer',

    function(obj, roi.string, model.name, shoulder, tf.count=NA){
       stopifnot(model.name %in% names(getModels(obj@singleGeneData)))

       roi <- parseChromLocString(roi.string)
       tbl.snps <- obj@singleGeneData@misc.data[["eqtl.snps"]]
       tbl.snps <- subset(tbl.snps, chrom==roi$chrom & pos >= roi$start & pos <= roi$end)
       rsids <- tbl.snps$rsid

       tfs <- tbl.model$gene
       if(!is.na(tf.count)){
          count <- min(length(tfs), tf.count)
          tfs <- head(tfs, n=count)
          }

       lookup <- function(gene){
          tbl.mdb <- geneToMotif(MotifDb, gene, "MotifDb", ignore.case=TRUE)
          tbl.tfc <- geneToMotif(MotifDb, gene, "TFClass", ignore.case=TRUE)
          motifs <- unique(c(tbl.mdb$motif, tbl.tfc$motif))
          motifs <- grep("^MA", motifs, value=TRUE)
          motifs
          }

       tf.motifs <- lapply(tfs, lookup)
       names(tf.motifs) <- tfs
       browser()
       xyz <- 99
       })

#----------------------------------------------------------------------------------------------------


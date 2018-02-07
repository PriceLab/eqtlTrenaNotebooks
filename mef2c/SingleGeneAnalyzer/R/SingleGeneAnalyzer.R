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

setGeneric('getWholeGenomeVariantsForRegion', signature='obj', function(obj, roi.string, altToRefRatio, minAltCount)
              standardGeneric ('getWholeGenomeVariantsForRegion'))

setGeneric('getDHSForRegion', signature='obj', function(obj, roi.string, score.threshold=NA) standardGeneric ('getDHSForRegion'))
setGeneric('getEnhancersForRegion',
           signature='obj', function(obj, roi.string, score.threshold=NA) standardGeneric ('getEnhancersForRegion'))

setGeneric('findVariantsInModelForRegion', signature='obj',
            function(obj, roi.string, motif.track, variants.track, candidate.tfs, tfMotifMappingName, shoulder=0)
               standardGeneric ('findVariantsInModelForRegion'))


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
setMethod('getWholeGenomeVariantsForRegion', 'SingleGeneAnalyzer',

   function(obj, roi.string, altToRefRatio, minAltCount){
      roi <- trena::parseChromLocString(roi.string)
      tbl.variants <- getWholeGenomeVariants(obj@singleGeneData, roi, altToRefRatio, minAltCount)
      invisible(tbl.variants)
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
             pfms <- as.list(query(query(MotifDb, motif), "jaspar2018"))
             printf("findMotifsInRegion.findMotif(%s), %d pfms from jaspar2018", motif, length(pfms))
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
             # put the incomping "motif" in the fifth column; first five columns then
             # follow good-enough bed format:  chrom start end name score strand
          preferred.column.order <- c("chrom","motifStart","motifEnd", "motif", "motifRelativeScore", "strand",
                                      "motifName", "motifScore", "match","chromStart","chromEnd","seq","status")
          tbl.all <- tbl.all[, preferred.column.order]
          colnames(tbl.all)[2:3] <- c("start", "end")
          tbl.all
          }) # findMotifsInR

#----------------------------------------------------------------------------------------------------
# strategy:
#   identify the tfs in the mode
#   get the motifs associated with each
#
setMethod('findVariantsInModelForRegion', 'SingleGeneAnalyzer',

    function(obj, roi.string, motif.track, variants.track, candidate.tfs, tfMotifMappingName, shoulder=0){

       roi <- parseChromLocString(roi.string)
       stopifnot(variants.track %in% c("eqtl.snps", "wgVariants"))
       stopifnot(tfMotifMappingName %in% c("MotifDb", "TFClass"))
       stopifnot(motif.track %in% c("DHS.motifs", "footprints", "enhancer.motifs", "allDNA.motifs"))

       browser()
       tbl.motifs <- switch(motif.track,
           "footprints" = obj@singleGeneData@tbl.fp,
           "enhancer.motifs" =
               {tbl.tmp <- obj@singleGeneData@misc.data$enhancer.motifs.mdb;
                tbl.tmp <- tbl.tmp[, c("chrom", "motifStart", "motifEnd", "motifName", "motifRelativeScore", "strand", "geneSymbol")];
                colnames(tbl.tmp) <- c("chrom", "start", "end", "name", "score", "strand", "geneSymbol");
                rownames(tbl.tmp) <- NULL;
                printf("enhancer.motifs: %d rows, %d geneSymbols", nrow(tbl.tmp), length(unique(tbl.tmp$geneSymbol)))
                tbl.motifs <- tbl.tmp
                },
           "allDNA.motifs" = NA,
           "DHS.motifs" =
              {tbl.tmp <- obj@singleGeneData@misc.data$tbl.dhsMotifs
               shortMotifNames <- unlist(lapply(tbl.tmp$motifName,
                               function(name) {tokens <- strsplit(name, "-")[[1]]; return(tokens[length(tokens)])}))
               tbl.tmp$shortMotif <- shortMotifNames
               tbl.tmp <- associateTranscriptionFactors(MotifDb, tbl.tmp, tfMotifMappingName, expand.rows=TRUE)
               tbl.tmp <- tbl.tmp[, c("chrom", "motifStart", "motifEnd", "motifName", "score", "strand", "geneSymbol")];
               colnames(tbl.tmp) <- c("chrom", "start", "end", "name", "score", "strand", "geneSymbol");
               rownames(tbl.tmp) <- NULL;
               printf("DHS.motifs: %d rows, %d geneSymbols", nrow(tbl.tmp), length(unique(tbl.tmp$geneSymbol)))
               tbl.motifs <- tbl.tmp
               }
            ) # switch on motif.track

       tbl.snps <- switch(variants.track,
          "eqtl.snps" = {
               tbl.tmp <- obj@singleGeneData@misc.data[["eqtl.snps"]][, c("chrom", "pos", "pos", "rsid", "score")];
               colnames(tbl.tmp) <- c("chrom", "start", "end", "rsid", "score")
               tbl.snps <- tbl.tmp
               },
          "wgVariants" = {
               tbl.tmp <- obj@singleGeneData@misc.data[["wgVariants"]]
               tbl.tmp$end <- tbl.tmp$pos;
               colnames.ordered <- c("chrom","pos","end","ref","alt","het.altAD","het.altCTL",
                                     "hom.altAD","hom.altCTL","any.altAD","any.altCTL")
               tbl.tmp <- tbl.tmp[, colnames.ordered]
               colnames(tbl.tmp)[2] <- "start"
               tbl.snps <- tbl.tmp
               }
           ) # switch on variants.track

          # find all motifs which intersect with the roi
       browser()

       gr.roi <- GRanges(as.data.frame(roi, stringsAsFactors=FALSE))
       gr.motifs <- GRanges(tbl.motifs[, 1:3])
       tbl.ov <- as.data.frame(findOverlaps(gr.roi, gr.motifs, type="any"))
       colnames(tbl.ov) <- c("roi", "motif")
       tbl.motifs <- tbl.motifs[tbl.ov$motif,]

       gr.motifs.roi <- GRanges(tbl.motifs)
       gr.snps <- GRanges(tbl.snps)
       tbl.ov <- as.data.frame(findOverlaps(gr.motifs.roi, gr.snps, type="any"))
       if(nrow(tbl.ov) == 0)
          return(data.frame())

       colnames(tbl.ov) <- c("motif", "snp")
       tbl.snpsInMotif.roi <- tbl.snps[tbl.ov$snp,]
       tbl.motifs.roi <- tbl.motifs[tbl.ov$motif,]
       tbl.snpsMotifs <- cbind(tbl.snpsInMotif.roi, tbl.motifs.roi)
       subset(tbl.snpsMotifs, geneSymbol %in% candidate.tfs)
       }) # findVariantsInModelForRegion

#----------------------------------------------------------------------------------------------------


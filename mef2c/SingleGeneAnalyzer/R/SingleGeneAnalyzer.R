.SingleGeneAnalyzer <- setClass ("SingleGeneAnalyzer",
                            representation = representation(
                               genomeName="character",
                               targetGene="character",
                               targetGene.TSS="numeric",
                               singleGeneData="SingleGeneData",
                               trena="Trena",
                               quiet="logical",
                               trackerCache="environment"
                               )
                            )

PORT <- 5548
#------------------------------------------------------------------------------------------------------------------------
setGeneric('summarizeExpressionMatrices', signature='obj', function(obj) standardGeneric ('summarizeExpressionMatrices'))
setGeneric('getRegulatoryModelNames', signature='obj', function(obj) standardGeneric ('getRegulatoryModelNames'))
setGeneric('getRegulatoryModel', signature='obj', function(obj, modelName) standardGeneric ('getRegulatoryModel'))

setGeneric('getFootprintsForRegion', signature='obj', function(obj, roi.string, score.threshold=NA)
              standardGeneric ('getFootprintsForRegion'))

setGeneric('getVariantsForRegion', signature='obj', function(obj, source.name, tracking.name, roi.string,
                  score.1.threshold=NA_real_, score.2.threshold=NA_real, score.3.threshold=NA_real_)
              standardGeneric ('getVariantsForRegion'))

setGeneric('getMotifsForRegion', signature='obj', function(obj, source.name, tracking.name, roi.string, score.threshold=NA_real_)
              standardGeneric('getMotifsForRegion'))

setGeneric('getDHSForRegion', signature='obj', function(obj, roi.string, score.threshold=NA) standardGeneric ('getDHSForRegion'))
setGeneric('getEnhancersForRegion',
           signature='obj', function(obj, roi.string, score.threshold=NA) standardGeneric ('getEnhancersForRegion'))

setGeneric('findVariantsInModelForRegion', signature='obj',
            function(obj, roi.string, motif.track, variants.source, candidate.tfs, tfMotifMappingName, shoulder=0)
               standardGeneric ('findVariantsInModelForRegion'))

setGeneric('findMotifsInRegion', signature='obj',
           function(obj, roi.string, motifs, pwmMatchPercentage, variants=NA_character)
              standardGeneric ('findMotifsInRegion'))

setGeneric('getCacheItemNames', signature='obj', function(obj) standardGeneric('getCacheItemNames'))
setGeneric('getFromCache', signature='obj', function(obj, tracker.name) standardGeneric('getFromCache'))
setGeneric('clearCache', signature='obj', function(obj) standardGeneric('clearCache'))
setGeneric('intersectTracks', signature='obj', function(obj, trackName.1, trackName.2, shoulder) standardGeneric('intersectTracks'))

#------------------------------------------------------------------------------------------------------------------------
SingleGeneAnalyzer = function(genomeName, targetGene, targetGene.TSS, singleGeneData, quiet=TRUE)
{
   trena <- Trena(genomeName)

   obj <- .SingleGeneAnalyzer(genomeName=genomeName,
                              targetGene=targetGene,
                              targetGene.TSS=targetGene.TSS,
                              singleGeneData=singleGeneData,
                              trena=trena,
                              quiet=quiet,
                              trackerCache=new.env(parent=emptyenv()))
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
setMethod('getCacheItemNames', signature='SingleGeneAnalyzer',

       function(obj) {
          return(ls(obj@trackerCache))
          })

#----------------------------------------------------------------------------------------------------
setMethod('getFromCache', signature='SingleGeneAnalyzer',

       function(obj, tracker.name){
          if(!tracker.name %in% ls(obj@trackerCache)){
             printf("nothing named '%s' in cache")
             return(NA)
             }
          obj@trackerCache[[tracker.name]]
          })

#----------------------------------------------------------------------------------------------------
setMethod('clearCache', signature='SingleGeneAnalyzer',

       function(obj){
          rm(list=ls(obj@trackerCache), envir=obj@trackerCache)
          })

#----------------------------------------------------------------------------------------------------
setMethod('getRegulatoryModelNames', 'SingleGeneAnalyzer',
       function(obj){
          return(names(getModels(obj@singleGeneData)))
       })

#----------------------------------------------------------------------------------------------------
setMethod('getRegulatoryModel', 'SingleGeneAnalyzer',
      function(obj, modelName){

         tbl <- data.frame()
         if(modelName %in% getRegulatoryModelNames(obj))
            tbl <- getModels(obj@singleGeneData)[[modelName]]
         return(tbl)
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

   function(obj, source.name, tracking.name, roi.string, score.1.threshold=NA_real_,
            score.2.threshold=NA_real_, score.3.threshold=NA_real_){

      roi <- trena::parseChromLocString(roi.string)
      tbl.variants <- getVariants(obj@singleGeneData, source.name, roi,
                                  score.1.threshold=score.1.threshold,
                                  score.2.threshold=score.2.threshold,
                                  score.3.threshold=score.3.threshold)
      if(!is.na(tracking.name))
         obj@trackerCache[[tracking.name]] <- tbl.variants

      invisible(tbl.variants)
      })

#----------------------------------------------------------------------------------------------------
setMethod('getMotifsForRegion', 'SingleGeneAnalyzer',

   function(obj, source.name, tracking.name, roi.string, score.threshold=NA_real_){

      roi <- trena::parseChromLocString(roi.string)
      tbl.motifs <- getMotifs(obj@singleGeneData, source.name, roi,
                                  score.threshold=score.threshold)
      if(!is.na(tracking.name))
         obj@trackerCache[[tracking.name]] <- tbl.motifs

      invisible(tbl.motifs)
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
setMethod('intersectTracks', signature='SingleGeneAnalyzer',

    function(obj, trackName.1, trackName.2, shoulder=0){
        printf("sga::intersectTracks, %s in %s with shoulder: %d",  trackName.1, trackName.2, shoulder)

        if(!trackName.1 %in% ls(obj@trackerCache)){
           printf("'%s% named variable not in SingleGeneAnalyzer cacher", trackName.1)
           return(data.frame())
           }
        if(!trackName.2 %in% ls(obj@trackerCache)){
           printf("'%s% named variable not in SingleGeneAnalyzer cacher", trackName.2)
           return(data.frame())
           }
        tbl.1 <- obj@trackerCache[[trackName.1]]
        tbl.2 <- obj@trackerCache[[trackName.2]]

        stopifnot(all(c("chrom", "start", "end") %in% colnames(tbl.1)))
        stopifnot(all(c("chrom", "start", "end") %in% colnames(tbl.2)))

        tbl.1$start <- tbl.1$start - shoulder
        tbl.1$end   <- tbl.1$end + shoulder

        tbl.ov <- as.data.frame(findOverlaps(GRanges(tbl.1), GRanges(tbl.2)))
        tbl.out <- data.frame()
        if(nrow(tbl.ov) > 0){
           colnames(tbl.ov) <- c("one", "two")
           tbl.1 <- tbl.1[tbl.ov$one,]
           tbl.2 <- tbl.2[tbl.ov$two,]
           colnames(tbl.2) <- sprintf("%s.B", colnames(tbl.2))
           tbl.out <- cbind(tbl.1, tbl.2)
           }
        return(tbl.out)
        }) # intersectTracks

#----------------------------------------------------------------------------------------------------
# strategy:
#   identify the tfs in the mode
#   get the motifs associated with each
#
setMethod('findVariantsInModelForRegion', 'SingleGeneAnalyzer',

    function(obj, roi.string, motif.track, variants.source, candidate.tfs, tfMotifMappingName, shoulder=0){

       roi <- parseChromLocString(roi.string)
       gr.roi <- GRanges(as.data.frame(roi, stringsAsFactors=FALSE))  # used on each motifs table below

       stopifnot(variants.source %in% c("eqtl.snps", "wgVariants"))
       stopifnot(tfMotifMappingName %in% c("MotifDb", "TFClass"))
       stopifnot(motif.track %in% c("DHS.motifs", "footprints", "enhancer.motifs", "allDNA.motifs"))

       tbl.motifs <- switch(motif.track,
          "footprints" =
              {tbl.tmp <- obj@singleGeneData@tbl.fp;
               gr.motifs <- with(tbl.tmp, GRanges(seqnames=chrom, IRanges(start=start, end=end)))
               tbl.ov <- as.data.frame(findOverlaps(gr.roi, gr.motifs, type="any"))
               colnames(tbl.ov) <- c("roi", "motif")
               tbl.tmp <- tbl.tmp[tbl.ov$motif,]
               colnames(tbl.tmp)[grep("^start$", colnames(tbl.tmp))] <- "motifStart";
               colnames(tbl.tmp)[grep("^end$", colnames(tbl.tmp))] <- "motifEnd";
               colnames(tbl.tmp)[grep("score", colnames(tbl.tmp))] <- "motifScore";
               tbl.motifs <- tbl.tmp;
               },
          "enhancer.motifs" =
              {tbl.tmp <- obj@singleGeneData@misc.data$enhancer.motifs.mdb;
               gr.motifs <- with(tbl.tmp, GRanges(seqnames=chrom, IRanges(start=motifStart, end=motifEnd)))
               tbl.ov <- as.data.frame(findOverlaps(gr.roi, gr.motifs, type="any"))
               colnames(tbl.ov) <- c("roi", "motif")
               tbl.tmp <- tbl.tmp[tbl.ov$motif,]
               tbl.tmp <- tbl.tmp[, c("chrom", "motifStart", "motifEnd", "motifName", "motifRelativeScore", "strand", "geneSymbol")];
               colnames(tbl.tmp) <- c("chrom", "motifStart", "motifEnd", "name", "motifScore", "strand", "geneSymbol");
               rownames(tbl.tmp) <- NULL;
               printf("enhancer.motifs: %d rows, %d geneSymbols", nrow(tbl.tmp), length(unique(tbl.tmp$geneSymbol)))
               tbl.motifs <- tbl.tmp
               },
           "allDNA.motifs" = NA,
           "DHS.motifs" =
              {tbl.tmp <- obj@singleGeneData@misc.data$tbl.dhsMotifs
               gr.motifs <- with(tbl.tmp, GRanges(seqnames=chrom, IRanges(start=motifStart, end=motifEnd)))
               tbl.ov <- as.data.frame(findOverlaps(gr.roi, gr.motifs, type="any"))
               colnames(tbl.ov) <- c("roi", "motif")
               tbl.tmp <- tbl.tmp[tbl.ov$motif,]
               shortMotifNames <- unlist(lapply(tbl.tmp$motifName,
                               function(name) {tokens <- strsplit(name, "-")[[1]]; return(tokens[length(tokens)])}))
               tbl.tmp$shortMotif <- shortMotifNames
               tbl.tmp <- associateTranscriptionFactors(MotifDb, tbl.tmp, tfMotifMappingName, expand.rows=TRUE)
               tbl.tmp <- tbl.tmp[, c("chrom", "motifStart", "motifEnd", "motifName", "score", "strand", "geneSymbol")];
               colnames(tbl.tmp) <- c("chrom", "motifStart", "motifEnd", "name", "motifScore", "strand", "geneSymbol");
               rownames(tbl.tmp) <- NULL;
               printf("DHS.motifs: %d rows, %d geneSymbols", nrow(tbl.tmp), length(unique(tbl.tmp$geneSymbol)))
               tbl.motifs <- tbl.tmp
               }
            ) # switch on motif.track

       #browser()

       tbl.snps <- switch(variants.source,
          "eqtl.snps" = {
               obj@singleGeneData@misc.data[["MAYO.eqtl.snps"]]
               },
          "wgVariants" = {
               obj@singleGeneData@misc.data[["wgVariants"]]
               #browser()
               #colnames.ordered <- c("chrom","start","end","id", "ref","alt","het.altAD","het.altCTL",
               #                      "hom.altAD","hom.altCTL","any.altAD","any.altCTL", "alleles_as_ambig")
               #tbl.tmp <- tbl.tmp[, colnames.ordered]
               #tbl.snps <- tbl.tmp
               }
           ) # switch on variants.source

          # find all motifs which intersect with the roi

       gr.motifs.roi <- GRanges(tbl.motifs)
       gr.snps <- GRanges(seqnames=tbl.snps$chrom, IRanges(start=tbl.snps$start-shoulder, end=tbl.snps$end+shoulder))
       tbl.ov <- as.data.frame(findOverlaps(gr.motifs.roi, gr.snps, type="any"))
       if(nrow(tbl.ov) == 0)
          return(data.frame())

       colnames(tbl.ov) <- c("motif", "snp")
       tbl.snpsInMotif.roi <- tbl.snps[tbl.ov$snp,]
       tbl.motifs.roi <- tbl.motifs[tbl.ov$motif,]
       tbl.snpsMotifs <- cbind(tbl.snpsInMotif.roi, tbl.motifs.roi)
       tbl.out <- subset(tbl.snpsMotifs, geneSymbol %in% candidate.tfs)
       tbl.out <- tbl.out[order(tbl.out$start, tbl.out$end),]
       rownames(tbl.out) <- NULL
       cols.to.remove <- sort(match(c("start.1", "end.1", "chrom.1"), colnames(tbl.out)))
       if(length(cols.to.remove) > 0)
          tbl.out <- tbl.out[, -cols.to.remove]
       preferred.colnames.in.order <- c("chrom", "start", "end", "id", "snpScore", "strand", "ref", "alt",
                                        "het.altAD", "het.altCTL", "hom.altAD", "hom.altCTL", "any.altAD",
                                        "any.altCTL", "seqnames", "alleles_as_ambig", "name",
                                        "motifStart", "motifEnd", "motifName", "motifScore", "geneSymbol")
       relevant.preferred.colnames.in.order <- intersect(preferred.colnames.in.order, colnames(tbl.out))
       tbl.out <- tbl.out[, relevant.preferred.colnames.in.order]
       dup.test.columns <- c("chrom", "start", "end", "id", "strand", "motifName", "motifScore")
       dup.test.columns <- intersect(dup.test.columns, colnames(tbl.out))
       dups <- which(duplicated(tbl.out[, dup.test.columns]))
       if(length(dups) > 0)
          tbl.out <- tbl.out[-dups,]
       tbl.out
       }) # findVariantsInModelForRegion

#----------------------------------------------------------------------------------------------------


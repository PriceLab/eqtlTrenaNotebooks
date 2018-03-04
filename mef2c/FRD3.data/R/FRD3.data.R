.FRD3.data <- setClass ("FRD3.data",
                         contains="SingleGeneData"
                         )

#------------------------------------------------------------------------------------------------------------------------
setGeneric('makeModelForRegion', signature='obj', function(obj, dhs.cutoff, region=NA, trenaViz=NA)
             standardGeneric('makeModelForRegion'))
setGeneric('motifTrackForTF', signature='obj',  function(obj, tbl.motifs, tf, trenaViz=NA) standardGeneric('motifTrackForTF'))
#------------------------------------------------------------------------------------------------------------------------
FRD3.data = function()
{
      # frd3.roi <- "3:2566277-2572151"
      # frd3.extended.roi <- "3:2,563,340-2,575,089"

   misc.data <- new.env(parent=emptyenv())

     #----------------------------------------------------------------------------------------------------
     #  load expression data, reduced as described in my log
     #  "cleaned up: create expression matrix from GEO project (a GSE) (30 dec 2017)"
     #  https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE77286
     #  Public on Jan 28, 2016: expression data from Arabidopsis plants under varying zinc supply.
     #----------------------------------------------------------------------------------------------------

   load(system.file(package="FRD3.data", "extdata", "mtx.zinc.22810x42.RData"))
   expression.matrices <- list()
   expression.matrices[["varying.zinc"]] <- mtx

     #----------------------------------------------------------------------------------------------------
     # load precalculated motifs from (only) the region of interest,
     #----------------------------------------------------------------------------------------------------

   load(system.file(package="FRD3.data", "extdata", "tbl.motifs.jaspar2018.athaliana.Chr3.2566277.2572151.ge85.RData"))
   misc.data[["motifs"]] <- list(jaspar2018.athaliana.Chr3.2566277.2572151.ge85 = tbl.motifs)
     # todo: move these to the base class and first-class slots
   misc.data[["TSS"]] <- list(primary=2569502,  secondary=2572149)
   misc.data[["targetGene"]] <- "AT3G08040"

     #----------------------------------------------------------------------------------------------------
     # load precalculated dhs region scores, for leaves and for buds
     #----------------------------------------------------------------------------------------------------

   load(system.file(package="FRD3.data", "extdata", "tbls.frd3.dhs.budAndLeaf.RData"))
   misc.data[["dhs"]] <- list(buds=tbl.frd3buds, leaves=tbl.frd3leaf)

     #----------------------------------------------------------------------------------------------------
     # MotifDb for athaliania maps motifs preferentially to uniprot ids.  load a variety of
     # cross-referencing tables calculated previously.  tbl.xref, at least, is a boon, mapping
     # uniprot id to athaliana orf
     # the createion of tables is document in my log:
     #     "my arabidopsis motif-to-tf mapping has gaps (31 dec 2017)"
     #----------------------------------------------------------------------------------------------------

   load(system.file(package="FRD3.data", "extdata", "tbl.xref.RData"))
   misc.data[["xref"]] <- list(tbl.upOrf=tbl.upOrf, tbl.orfSym=tbl.orfSym,
                               tbl.upOldgene=tbl.upOldgene, tbl.xref=tbl.xref)

   obj <- .FRD3.data(SingleGeneData(chrom="chr3",
                                    start=2566277,
                                    end=2572151,
                                    tbl.fp=data.frame(),
                                    expression.matrices=expression.matrices,
                                    models=list(),
                                    misc.data=misc.data))


   obj

} # constructor
#----------------------------------------------------------------------------------------------------
# generalize, put in base class?
setMethod('getMotifs', 'FRD3.data',

    function(obj, source.name, roi, score.threshold=NA_real_){

      available.motifs <- names(obj@misc.data$motifs)
      stopifnot(source.name %in% available.motifs)
      tbl.motifs <- obj@misc.data$motifs[[source.name]]
      tbl.out <- subset(tbl.motifs, tolower(chrom)==tolower(roi$chrom) & start >= roi$start & end < roi$end)
      if(!is.na(score.threshold))
        tbl.out <- subset(tbl.out, score >= score.threshold)
      tbl.out
      })

#----------------------------------------------------------------------------------------------------
setMethod('makeModelForRegion', 'FRD3.data',

    function(obj, dhs.cutoff, region=NA, trenaViz=NA){
       tss <- obj@misc.data$TSS
       if(is.na(region)){
          if(is.na(trenaViz)){
             printf("if not supplying explicit genomic region, must supply trenaViz so that can be queried")
             return(NA)
             }
          region <- getGenomicRegion(trenaViz)
          }
       region <- parseChromLocString(region)
       tbl.croi <- with(region, data.frame(chrom=chrom, start=start, end=end, stringsAsFactors=FALSE))
       tbl.croi$chrom <- sub("chr", "Chr", tbl.croi$chrom)

       tbl.frd3buds  <- obj@misc.data$dhs$buds
       minSpan <- 5
       peaks <- findRegionsAboveThreshold(tbl.frd3buds$score, threshold=dhs.cutoff, minSpan)

       peak.position.offset <- tbl.frd3buds$start[1]

       tbl.dhs <- with(peaks, data.frame(chrom="Chr3",
                                         start=starts+peak.position.offset,
                                         end=ends+peak.position.offset,
                                         stringsAsFactors=FALSE))

       tbl.ov <- as.data.frame(findOverlaps(GRanges(tbl.croi), GRanges(tbl.dhs), type="any"))
       colnames(tbl.ov) <- c("croi", "dhs")
       if(nrow(tbl.ov) == 0){
          printf("no DHS regions")
          return(NA)
          }

       tbl.dhs.currentView <- tbl.dhs[tbl.ov$dhs,]

       track.title <- sprintf("bud dhs > %f", dhs.cutoff)
       if(!is.na(trenaViz)){
          raiseTab(trenaViz, "IGV")
          addBedTrackFromDataFrame(trenaViz, track.title, tbl.dhs.currentView, color="darkGreen", trackHeight=50)
          }

       pfms <- as.list(query(query(MotifDb, "jaspar2018"), "athaliana"))
       mm <- MotifMatcher("tair10", pfms)
       tbl.motifs <- findMatchesByChromosomalRegion(mm, tbl.dhs.currentView, pwmMatchMinimumAsPercentage=90)

       dim(tbl.motifs)
       uniprotIDs <- mcols(MotifDb[tbl.motifs$motifName])$proteinId
       tbl.motifs$uniprot <- uniprotIDs

       tbl.xref   <- obj@misc.data$xref$tbl.xref
       tbl.orfSym <- obj@misc.data$xref$tbl.orfSym

       match.indices <- match(tbl.motifs$uniprot, tbl.xref$proteinId)
       tbl.motifs$tf <- tbl.xref$orfs[match.indices]

       mtx <- getExpressionMatrices(obj)$varying.zinc

       target.orf <- obj@misc.data$targetGene
       stopifnot(target.orf %in% rownames(mtx))

       tf.candidates <- tbl.motifs$tf
       colnames(tbl.motifs)[grep("tf", colnames(tbl.motifs))] <- "geneSymbol"
       colnames(tbl.motifs)[grep("motifStart", colnames(tbl.motifs))] <- "start"
       colnames(tbl.motifs)[grep("motifEnd", colnames(tbl.motifs))] <- "end"

       solver.names <- c("lasso", "pearson", "randomForest", "ridge", "spearman")

       if(!exists("trena"))
          trena <- Trena("mm10")   # a dodge, one which works if getProximalPromoter method is avoided
       tbl.model <- createGeneModel(trena, target.orf, solver.names, tbl.motifs, mtx)
       match.indices <- match(tbl.model$gene, tbl.orfSym$TAIR)
       tbl.model$sym <- tbl.orfSym$SYMBOL[match.indices]
       mapping.failures <- which(is.na(tbl.model$sym))

       if(length(mapping.failures) > 0)
          tbl.model$sym[mapping.failures] <- tbl.model$gene[mapping.failures]
       tbl.model <- tbl.model[order(tbl.model$rfScore, decreasing=TRUE),]
       tbl.model.strong <- head(tbl.model, n=10)
       tbl.motifs.strong <- subset(tbl.motifs, geneSymbol %in% tbl.model.strong$gene)

       distance <- tbl.motifs.strong$start - tss$primary
       direction <- rep("upstream", length(distance))
       direction[which(distance < 0)] <- "downstream"
       tbl.motifs.strong$distance.from.tss <- distance
       tbl.motifs.strong$id <- sprintf("%s.dhs.%s.%06d.%s", target.orf, direction, abs(distance), tbl.motifs.strong$motifName)

       return(list(model=tbl.model.strong, regions=tbl.motifs.strong))
       }) # makeModelForREgion

#----------------------------------------------------------------------------------------------------
setMethod('motifTrackForTF', 'FRD3.data',

    function(obj, tbl.motifs, tf, trenaViz=NA){
      tbl.bed <- subset(tbl.motifs, geneSymbol==tf)[, c("chrom", "start", "end", "motifName", "motifRelativeScore")]
      if(nrow(tbl.bed) == 0){
         printf("FRD3.data::motifTrackForTF error: no motifs for tf '%s'", tf)
         return(tbl.bed)
         }
      getLastToken <- function(longMotifName){
         tokens <- strsplit(longMotifName, "-")[[1]];
         return(tokens[length(tokens)])
         }
      names <-unlist(lapply(tbl.bed$motifName, getLastToken))
      tbl.bed$motifName <- names
      colnames(tbl.bed) <- c("chrom", "start", "end", "name", "score")
      rownames(tbl.bed) <- NULL
      if(!is.na(trenaViz)){
         raiseTab(trenaViz, "IGV")
         track.title <- sprintf("%s motifs", tf)
         addBedTrackFromDataFrame(trenaViz, track.title, tbl.bed, color="purple", trackHeight=50)
         }
      invisible(tbl.bed)
      })

#----------------------------------------------------------------------------------------------------
findRegionsAboveThreshold <- function(vec, threshold, minSpan, quiet=TRUE)
{
   vec[vec < threshold] <- NA
   vec[vec >= threshold] <- 1
   vec.rle <- rle(vec)
   rle.segments.above.threshold <- which(vec.rle$lengths >= minSpan)
       # make all runs less than desired.minSpan are NA'd
   all.subMinimumSpans <- setdiff(1:length(vec.rle$values), rle.segments.above.threshold)
   vec.rle$values[all.subMinimumSpans] <- NA    # vec.rle$values[18] <- NA

   rle.region.count <- length(vec.rle$values)
   actual.index <- 1

   peak.starts <- vector("numeric", length(vec))
   peak.ends <- vector("numeric", length(vec))
   peak.count <- 0

   for(i in 1:rle.region.count){
      size <- vec.rle$length[i]
      value <- vec.rle$values[i]
      if(!is.na(value)){
         peak.count <- peak.count + 1
         region.start <- actual.index
	 region.end   <- actual.index + size - 1
	 if(!quiet) printf("peak found:  %d-%d", region.start, region.end)
         peak.starts[peak.count] <- region.start
         peak.ends[peak.count] <- region.end
	 } # !is.na
       actual.index <- actual.index + size
       } # for i

   list(starts=peak.starts[1:peak.count], ends=peak.ends[1:peak.count])

} # findRegionsAboveThreshold
#------------------------------------------------------------------------------------------------------------------------

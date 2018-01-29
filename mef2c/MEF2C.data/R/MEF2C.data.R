.MEF2C.data <- setClass ("MEF2C.data",
                         contains="SingleGeneData"
                         )

PORT <- 5548
#------------------------------------------------------------------------------------------------------------------------
setGeneric('getGenomicBounds', signature='obj', function(obj) standardGeneric ('getGenomicBounds'))
setGeneric('getExpressionMatrices', signature='obj', function(obj) standardGeneric ('getExpressionMatrices'))
setGeneric('getFootprints', signature='obj', function(obj, roiString) standardGeneric ('getFootprints'))
#------------------------------------------------------------------------------------------------------------------------
MEF2C.data = function()
{
    mtx.names <- load(system.file(package="MEF2C.data", "extData", "mtx.withDimers.cer.ros.tcx.RData"))
    expression.matrices <- list()
    for(matrix.name in mtx.names){
       mtx <- eval(parse(text=matrix.name))
       expression.matrices[[matrix.name]] <- mtx
       }
   load(system.file(package="MEF2C.data", "extData", "tbl.fp.chr5.88615025-89052115.4sources.noDups.RData"))

   obj <- .MEF2C.data(SingleGeneData(chrom="chr5", start=88618415, end=89052071, tbl.fp, expression.matrices))
   obj

} # constructor
#----------------------------------------------------------------------------------------------------
setMethod('getFootprints', 'MEF2C.data',

    function(obj, roiString){
      roi <- trena::parseChromLocString(roiString)   # a trena function
      if(roi$chrom == obj@tbl.fp$chrom &
         roi$start >= obj@tbl.fp$start &
         roi$end   <= obj@tbl.fp$end) {
         tbl.roi <- subset(obj@tbl.fp, chrom==roi$chrom & start >= roi$start & end <= roi$end)
         return(tbl.roi)
         }

      printf("MEF2C.data::getFootprints, roi not completely contained in precalculated footprints");

      # source.1 <- "postgres://bddsrds.globusgenomics.org/brain_wellington_16"
      # source.2 <- "postgres://bddsrds.globusgenomics.org/brain_wellington_20"
      # source.3 <- "postgres://bddsrds.globusgenomics.org/brain_hint_16"
      # source.4 <- "postgres://bddsrds.globusgenomics.org/brain_hint_20"
      # sources <- c(source.1, source.2, source.3, source.4)
      # names(sources) <- c("well_16", "well_20", "hint_16", "hint_20")

     #  x <- getRegulatoryChromosomalRegions(trena, roi$chrom, roi$start, roi$end, sources, targetGene, targetGene.tss)
      # print(1)
      # names(x) <- names(sources)
      # print(2)

      #     # append a column to each non-empty table, giving it the source name
      # x2 <- lapply(names(x), function(name) {tbl <-x[[name]]; if(nrow(tbl) >0) tbl$db <- name; return(tbl)})
      # print(3)

      #  tbl.reg <- do.call(rbind, x2)
      # print(4)
      # rownames(tbl.reg) <- NULL
      # print(5)

      #     # be strict for now: just the 2016 jaspar human motifs
      # tbl.reg <- unique(tbl.reg[grep("Hsapiens-jaspar2016", tbl.reg$motifName, ignore.case=TRUE),])
      # tbl.reg <- tbl.reg[order(tbl.reg$motifStart),]
      # printf("---- getFootprints, before associateTranscriptionFactors")
      # print(head(tbl.reg))
      # tbl.reg <- associateTranscriptionFactors(MotifDb, tbl.reg, source="MotifDb", expand.rows=TRUE)

      #  colnames(tbl.reg)[2:3] <- c("start", "end")
      #invisible(tbl.reg)
      return(data.frame())
      }) # getFootprints

#----------------------------------------------------------------------------------------------------

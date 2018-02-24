.FRD3.data <- setClass ("FRD3.data",
                         contains="SingleGeneData"
                         )

PORT <- 5548
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

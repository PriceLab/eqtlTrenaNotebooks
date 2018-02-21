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

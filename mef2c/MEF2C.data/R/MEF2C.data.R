.MEF2C.data <- setClass ("MEF2C.data",
                         contains="SingleGeneData"
                         )

PORT <- 5548
#------------------------------------------------------------------------------------------------------------------------
setGeneric('getGenomicBounds', signature='obj', function(obj, asString=FALSE) standardGeneric ('getGenomicBounds'))
setGeneric('getExpressionMatrices', signature='obj', function(obj) standardGeneric ('getExpressionMatrices'))
setGeneric('getFootprints', signature='obj', function(obj, roi) standardGeneric ('getFootprints'))
setGeneric('getEnhancers', signature='obj', function(obj, roi) standardGeneric ('getEnhancers'))
#------------------------------------------------------------------------------------------------------------------------
MEF2C.data = function()
{
       #--------------------------------------------------------------------------------
       # 3 normalized matrices, each with virtual dimer expression added
       #--------------------------------------------------------------------------------

    mtx.names <- load(system.file(package="MEF2C.data", "extdata", "mtx.withDimers.cer.ros.tcx.RData"))
    expression.matrices <- list()
    for(matrix.name in mtx.names){
       mtx <- eval(parse(text=matrix.name))
       expression.matrices[[matrix.name]] <- mtx
       }

       #--------------------------------------------------------------------------------
       # all hint/wellington 16/20 footprints
       #--------------------------------------------------------------------------------

  # load(system.file(package="MEF2C.data", "extdata", "tbl.fp.chr5.88615025-89052115.4sources.noDups.RData")) # tbl.fp
   load(system.file(package="MEF2C.data", "extdata", "tbl.fp.chr5.88391000-89322000.4sources.noDups.RData")) #

       #--------------------------------------------------------------------------------
       # enhancers and their motifs
       #--------------------------------------------------------------------------------

   load(system.file(package="MEF2C.data", "extdata", "tbl.mef2c.eLocs.RData")) # tbl.mef2c.eLocs, dim 22 3

   load(system.file(package="MEF2C.data", "extdata", "tbl.eMotifs.RData"))
         #  tbl.eMotifs.mdb, [1]  43076    15
         #  tbl.eMotifs.tfc, [1] 891914     15
         # all motifs with motifRelativeScore >= 0.85, quite permissive

   misc.data <- new.env(parent=emptyenv())
   misc.data[["enhancer.locs"]] <- tbl.mef2c.eLocs
   misc.data[["enhancer.motifs.mdb"]] <- tbl.eMotifs.mdb
   misc.data[["enhancer.motifs.tfc"]] <- tbl.eMotifs.tfc

       #--------------------------------------------------------------------------------
       # vcf lifted over to hg38.  see extdata/README.txt, "vcf data reduction"
       # per-sample phenotype data loaded here too
       #--------------------------------------------------------------------------------

   load(system.file(package="MEF2C.data", "extdata", "mtx.tcx.pheno.geno.RData"))
   misc.data[["vcf.geno"]] <- tbl.geno
   misc.data[["pheno"]]    <- tbl.pheno

       #--------------------------------------------------------------------------------
       # eqtl snps, from mayo.AD.collaboration/tbl.snp.hg38.score-ref-alt.RData
       # dim(tbl.snps) 155 18
       #--------------------------------------------------------------------------------

   load(system.file(package="MEF2C.data", "extdata", "tbl.snp.hg38.score-ref-alt.RData"))
   misc.data[["eqtl.snps"]] <- tbl.snp

       #--------------------------------------------------------------------------------
       # dhs regions, then their motifs, from UCSC-hosted wgEncodeRegDnaseClustered
       #--------------------------------------------------------------------------------

   load(system.file(package="MEF2C.data", "extdata", "tbl.dhs.RData"))
   misc.data[["tbl.dhs"]] <- tbl.dhs
   load(system.file(package="MEF2C.data", "extdata", "tbl.dhsMotifs.RData"))
   misc.data[["tbl.dhsMotifs"]] <- tbl.dhsMotifs

       #--------------------------------------------------------------------------------
       # load all precalculated trena models
       #--------------------------------------------------------------------------------

   model.names <- load(system.file(package="MEF2C.data", "extdata", "mef2c.tf.5kb.RData"))
   gene.models <- list()
   for(name in model.names){  # 3 data.frames, 20 x 12
     gene.models[[name]] <- eval(parse(text=name))
     }


   obj <- .MEF2C.data(SingleGeneData(chrom="chr5",
                                     start=88391000,  # see extdata/README.txt, padded enhancer span: "chr5:88391337-89321026"
                                       end=89322000,
                                     tbl.fp=tbl.fp,
                                     expression.matrices=expression.matrices,
                                     models=gene.models,
                                     misc.data=misc.data))


   obj

} # constructor
#----------------------------------------------------------------------------------------------------
# TODO: get intersection (perhaps empty) of roi with tbl.fp, then make live calls
# TODO: to the footprint databases to fill in the missing regions

setMethod('getFootprints', 'MEF2C.data',

    function(obj, roi){
      tbl.roi <- subset(obj@tbl.fp, chrom==roi$chrom & start >= roi$start & end <= roi$end)
      if(nrow(tbl.roi) > 0)
         return(tbl.roi)

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

.MEF2C.data <- setClass ("MEF2C.data",
                         contains="SingleGeneData"
                         )

PORT <- 5548
#------------------------------------------------------------------------------------------------------------------------
setGeneric('getGenomicBounds', signature='obj', function(obj, asString=FALSE) standardGeneric ('getGenomicBounds'))
setGeneric('getExpressionMatrices', signature='obj', function(obj) standardGeneric ('getExpressionMatrices'))
setGeneric('getFootprints', signature='obj', function(obj, roi) standardGeneric ('getFootprints'))
setGeneric('getEnhancers', signature='obj', function(obj, roi) standardGeneric ('getEnhancers'))
setGeneric('makeModelForRegion', signature='obj', function(obj, expression.matrix.name, region.string=NA, trenaViz=NA)
              standardGeneric('makeModelForRegion'))
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
   misc.data[["targetGene"]] <- "MEF2C"
   misc.data[["TSS"]] <- 88884466

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
      # now the vcf data for the bounds() region, reporting relative counts at each variant base.
      # subset(tbl.pos, any.altAD > (2.5*any.altCTL) &  any.altAD > 10)[, c(1,2,5:10)]
      #      chrom      pos het.altAD het.altCTL hom.altAD hom.altCTL any.altAD any.altCTL
      # 108   chr5 88410266        11          3         1          0        12          3
      # 1422  chr5 88599367        15          5         0          0        15          5
      # 2622  chr5 88783554        13          4         0          1        13          5
      # 2745  chr5 88799905        13          4         0          1        13          5
      # 4082  chr5 89004814        13          4         0          0        13          4
      # 5226  chr5 89208216         7          2         5          2        12          4

   load(system.file(package="MEF2C.data", "extdata", "tbl.vcf.chr5.88391000.89322000.79AD.73CTL.RData"))
   misc.data[["wgVariants"]] <- tbl.pos

       #--------------------------------------------------------------------------------
       # eqtl snps, from mayo.AD.collaboration/tbl.snp.hg38.score-ref-alt.RData
       # dim(tbl.snps) 155 18
       #--------------------------------------------------------------------------------

   load(system.file(package="MEF2C.data", "extdata", "tbl.snp.hg38.score-ref-alt.bedFormat.RData"))
   misc.data[["MAYO.eqtl.snps"]] <- tbl.snp

       #--------------------------------------------------------------------------------
       # dhs regions, then their motifs, from UCSC-hosted wgEncodeRegDnaseClustered
       #--------------------------------------------------------------------------------

   load(system.file(package="MEF2C.data", "extdata", "tbl.dhs.RData"))
   misc.data[["tbl.dhs"]] <- tbl.dhs
   load(system.file(package="MEF2C.data", "extdata", "tbl.dhsMotifs.RData"))
   misc.data[["tbl.dhsMotifs"]] <- tbl.dhsMotifs

       #--------------------------------------------------------------------------------
       # load igap variants (228 in this region)
       #--------------------------------------------------------------------------------
   load(system.file(package="MEF2C.data", "extdata", "tbl.igapVariants.mef2c.RData"))
   misc.data[["IGAP.snpChip"]] <- tbl.igapVariants.mef2c

       #--------------------------------------------------------------------------------
       # load all-DNA motifs, jaspar2018 human + mouse
       #--------------------------------------------------------------------------------
   load(system.file(package="MEF2C.data", "extdata",
                    "tbl.motifs.jaspar2018.human.mouse.chr5.88391000.89322000.RData"))
   misc.data[["allDNA-jaspar2018-human-mouse-motifs"]] <- tbl.motifs

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
setMethod('getVariants', 'MEF2C.data',

    function(obj, source.name, roi, score.1.threshold=NA_real_,
                  score.2.threshold=NA_real_, score.3.threshold=NA_real_){

      stopifnot(source.name %in% c("ADNI.WGS", "IGAP.snpChip", "MAYO.eqtl.snps"))

      if(source.name == "ADNI.WGS"){
         tbl.pos <- obj@misc.data[["wgVariants"]]
         tbl.out <- subset(tbl.pos, chrom==roi$chrom & start >= roi$start & end <= roi$end)
         if(nrow(tbl.out) == 0)
            return(data.frame())
         if(!is.na(score.1.threshold)){
            altToRefRatio <- score.1.threshold
            tbl.out <- subset(tbl.out, any.altAD > (altToRefRatio * any.altCTL))
            if(nrow(tbl.out) == 0)
               return(data.frame())
            } # AD/CTL > score.1.threshold
         if(!is.na(score.2.threshold)){
            altAnyCountMinimum <- score.2.threshold
            tbl.out <- subset(tbl.out, any.altAD > altAnyCountMinimum)
            if(nrow(tbl.out) == 0)
               return(data.frame())
            } # AD/CTL

         tbl.out$name <- with(tbl.out, sprintf("%d.%s.%s.%d.%d.%d.%d",
                              start, het.altAD, het.altCTL, hom.altAD, hom.altCTL, any.altAD, any.altCTL))
         tbl.out$score <- with(tbl.out, any.altAD/any.altCTL)

         tbl.out <- tbl.out[, c("chrom", "start", "end", "name", "score")]
         } # adni.wgs
      else if(source.name == "IGAP.snpChip"){
         tbl.var <- obj@misc.data[["IGAP.snpChip"]]
         tbl.out <- subset(tbl.var, chrom==roi$chrom & start >= roi$start & end <= roi$end)
         if(nrow(tbl.out) == 0)
            return(data.frame())
         if(!is.na(score.1.threshold)){
            minusLog10.pval <- score.1.threshold
            tbl.out <- subset(tbl.out, -log10(pval) >= minusLog10.pval)
            } # threshold supplied
         } # igap.snpchip
      else if(source.name == "MAYO.eqtl.snps"){
         tbl.eqtl <- obj@misc.data$MAYO.eqtl.snps
         tbl.out <- subset(tbl.eqtl, chrom==roi$chrom & start >= roi$start & start <= roi$end)
         if(!is.na(score.1.threshold))
            tbl.out <- subset(tbl.out, -log10(CER_P) >= score.1.threshold)
         } # MAYO.eqtl.snps
      tbl.out
      }) # getVariants

#----------------------------------------------------------------------------------------------------
setMethod('getMotifs', 'MEF2C.data',

    function(obj, source.name, roi, score.threshold=NA_real_){

      stopifnot(source.name %in% c("allDNA-jaspar2018-human-mouse-motifs"))
       if(source.name == "allDNA-jaspar2018-human-mouse-motifs"){
          tbl.motifs <- obj@misc.data[["allDNA-jaspar2018-human-mouse-motifs"]]
          tbl.out <- subset(tbl.motifs, chrom==roi$chrom & start >= roi$start & end < roi$end)
          if(!is.na(score.threshold))
             tbl.out <- subset(tbl.out, score >= score.threshold)
          } # allDNA-jaspar2018-human-mouse-motifs

       tbl.out
       })

#----------------------------------------------------------------------------------------------------
# build a model using all footprints in enhancers in the currently displayed region
setMethod('makeModelForRegion', 'MEF2C.data',

    function(obj, expression.matrix.name, region.string=NA, trenaViz=NA){
       stopifnot(expression.matrix.name %in% c("mtx.cer", "mtx.ros","mtx.tcx"))
       mtx <- getExpressionMatrices(obj)[[expression.matrix.name]]
       tss <- obj@misc.data$TSS
       browser()
       if(is.na(region.string)){
          if(is.na(trenaViz)){
             printf("if not supplying explicit genomic region, must supply trenaViz so that can be queried")
             return(NA)
             }
          region.string <- getGenomicRegion(trenaViz)
          }
       region <- parseChromLocString(region.string)

           #--------------------------------------------------------------------------------
           # consider all the footprints in enhancer regions in region
           #--------------------------------------------------------------------------------

       tbl.enhancers <- obj@misc.data$enhancer.locs
       tbl.enhancersClipped <- subset(tbl.enhancers, chrom=region$chrom & start >= region$start & end <= region$end)
       tbl.fp <- getFootprints(mef2c, region)
       tbl.ov <- as.data.frame(findOverlaps(GRanges(tbl.fp), GRanges(tbl.enhancers)))
       colnames(tbl.ov) <- c("fp", "enhancers")
       tbl.fpe <- tbl.fp[unique(tbl.ov$fp),]

       tfClass.candidate.tfs <- unique(tbl.fpe$geneSymbol)
       tfClass.candidate.tfs <- intersect(tfClass.candidate.tfs, rownames(mtx))

       trena <- Trena("hg38")
       solver.names <- c("lasso", "pearson", "randomForest", "ridge", "spearman")
       tbl.model.tfClass <- createGeneModel(trena, "MEF2C", solver.names, tbl.fpe, mtx)

       pfms.human <- as.list(query(query(MotifDb, "jaspar2018"), "hsapiens"))
       pfms.mouse <- as.list(query(query(MotifDb, "jaspar2018"), "mmusculus"))
       pfms <- c(pfms.human, pfms.mouse)
       mm <- MotifMatcher("hg38", pfms)

       tbl.fpe.regions <- unique(tbl.fpe[, c("chrom", "start", "end")])
       tbl.eMotifs.mdb <- obj@misc.data$enhancer.motifs.mdb[, c("chrom", "motifStart", "motifEnd", "motifName", "shortMotif")]
       colnames(tbl.eMotifs.mdb) <- c("chrom", "start", "end", "motifName", "shortMotif")
       tbl.ov <- as.data.frame(findOverlaps(GRanges(tbl.eMotifs.mdb), GRanges(tbl.fpe.regions), type="any"))  # any|within
       colnames(tbl.ov) <- c("motifs", "fpEnhancers")
       tbl.motifs <- tbl.eMotifs.mdb[unique(tbl.ov$motifs),]
       dim(tbl.motifs)
       tbl.motifs <- associateTranscriptionFactors(MotifDb, tbl.motifs, source="MotifDb")
       tbl.model.motifDb <- createGeneModel(trena, "MEF2C", solver.names, tbl.motifs, mtx)

       return(list(tfClassFp=list(model=tbl.model.tfClass, motifs=tbl.fpe),
                     motifDb=list(model=tbl.model.motifDb, motifs=tbl.motifs)))
       }) # makeModelForREgion

#------------------------------------------------------------------------------------------------------------------------

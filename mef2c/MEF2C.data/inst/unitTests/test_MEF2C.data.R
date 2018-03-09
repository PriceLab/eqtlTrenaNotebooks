library(MEF2C.data)
library(RUnit)
library(GenomicRanges)
#----------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#----------------------------------------------------------------------------------------------------
if(!exists("mef2c")){
   mef2c <- MEF2C.data()    # just load once, speeding up the tests
   checkTrue(all(c("SingleGeneData", "MEF2C.data") %in% is(mef2c)))
   }

if(exists("viz")){
   if(viz){
      library(trenaViz)
      if(!exists("tv")){
         PORTS <- 10000:10020   # trenaViz claims a websocket port in this range, over which R, igv and cytoscape communicate
         tv <- trenaViz(portRange=PORTS)
         setGenome(tv, "hg38")
         Sys.sleep(3)
         showGenomicRegion(tv, "MEF2C")
         Sys.sleep(3)
         addBedTrackFromDataFrame(tv, "enhancers", mef2c@misc.data$enhancer.locs, color="black")
         }
     } # if viz
   } # if exists

#----------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_getGenomicBounds()
   test_getExpressionMatrices()
   test_getFootprints()
   test_getModels()
   test_getVariants()
   test_getMotifs()
   test_buildModels()
   test_intersectRegions_findVariantsInModel()

} # runTests
#------------------------------------------------------------------------------------------------------------------------
test_getGenomicBounds <- function()
{
   printf("--- test_getGenomicBounds")

   roi <- getGenomicBounds(mef2c)
   with(roi, checkEquals(roi, list(chrom="chr5", start=88391000, end=89322000)))

   roi.string <- getGenomicBounds(mef2c, asString=TRUE)
   checkEquals(roi.string, "chr5:88391000-89322000")

} # test_getGenomicBounds
#------------------------------------------------------------------------------------------------------------------------
test_getExpressionMatrices <- function()
{
    printf("--- test_getExpressionMatrices")

    x <- getExpressionMatrices(mef2c)
    mtx.names <- c("mtx.cer", "mtx.ros", "mtx.tcx")
    checkEquals(sort(names(x)), mtx.names)
    dims <- lapply(x, dim)
    checkTrue(all(lapply(dims, length) == 2))

} # test_getExpressionMatrices
#----------------------------------------------------------------------------------------------------
test_getFootprints <- function()
{
    printf("--- test_getFootprints")

    roi <- getGenomicBounds(mef2c)
    tbl.fp <- getFootprints(mef2c, roi)
    checkEquals(dim(tbl.fp), c(29397, 13))

} # test_getFootprints
#----------------------------------------------------------------------------------------------------
test_getModels <- function()
{
    printf("--- test_getModels")

    model.list <- getModels(mef2c)
    checkEquals(sort(names(model.list)),
                c("mef2c.cory.wgs.cer.tfClass", "mef2c.cory.wgs.ros.tfClass", "mef2c.cory.wgs.tcx.tfClass"))

    checkTrue(all(lapply(model.list, class) == "data.frame"))
    checkTrue(all(lapply(model.list, ncol) > 5))
    checkTrue(all(lapply(model.list, nrow) > 5))

} # test_getModels
#----------------------------------------------------------------------------------------------------
test_getWholeGenomeVariants <- function()
{
   printf("--- test_getWholeGenomeVariants")
   roi <- getGenomicBounds(mef2c)
   tbl <- getWholeGenomeVariants(mef2c, roi, altToRefRatio=2.5, minAltCount=10)
   checkEquals(dim(tbl), c(6, 5))
   checkEquals(lapply(tbl, class),
               list(chrom="character",
                    start="integer",
                    end="integer",
                    name="character",
                    score="numeric"))
   checkTrue(all(tbl$score >= 2.5))

} # test_getWholeGenomeVariants
#----------------------------------------------------------------------------------------------------
# we currently have three kinds of variants
test_getVariants <- function()
{
   printf("--- test_getVariants")

      # first the IGAP snpChip variants
   roi <- getGenomicBounds(mef2c)
   tbl.igap <- getVariants(mef2c, "IGAP.snpChip", roi, score.1.threshold=2.5)
   checkEquals(dim(tbl.igap), c(25, 11))

      # now the ADNI whole genome sequencing results
   tbl.adni <- getVariants(mef2c, "ADNI.WGS", roi)
   checkTrue(nrow(tbl.adni) > 5000)

   tbl.adni.2 <- getVariants(mef2c, "ADNI.WGS", roi, score.1.threshold=3)
   checkEqualsNumeric(nrow(tbl.adni.2), 1000, tol=10)

   tbl.adni.3 <- getVariants(mef2c, "ADNI.WGS", roi, score.2.threshold=30)  # AD samples, het or hom, >= this
   checkEqualsNumeric(nrow(tbl.adni.3), 730, tol=10)

   tbl.adni.4 <- getVariants(mef2c, "ADNI.WGS", roi, score.1.threshold=2, score.2.threshold=10)  # AD samples, het or hom, >= this
   checkEqualsNumeric(nrow(tbl.adni.4), 13, tol=5)

      # MAYO.eqtl.snps
   tbl.eqtl <- getVariants(mef2c, "MAYO.eqtl.snps", roi, score.1.threshold=2)
   checkTrue(nrow(tbl.eqtl) > 50)
   checkTrue(nrow(tbl.eqtl) < 60)
   checkTrue(all(-log10(tbl.eqtl$CER_P) >= 2))

      #--------------------------------------------------------------------------------
      # find a small region mentioned in all three sources
      #--------------------------------------------------------------------------------

   roi <- getGenomicBounds(mef2c)
   tbl.igap <- getVariants(mef2c, "IGAP.snpChip", roi)
   tbl.adni <- getVariants(mef2c, "ADNI.WGS", roi)
   tbl.eqtl <- getVariants(mef2c, "MAYO.eqtl.snps", roi)
   gr.igap <- GRanges(tbl.igap)
   gr.adni <- GRanges(tbl.adni)
   gr.eqtl <- GRanges(tbl.eqtl)

   tbl.ov <- as.data.frame(findOverlaps(gr.igap, gr.adni))
   colnames(tbl.ov) <- c("igap", "adni")

   tbl.ov2 <- as.data.frame(findOverlaps(gr.igap, gr.eqtl))
   colnames(tbl.ov2) <- c("igap", "eqtl")

   rsid.loc <- 88820502    # rs159950
   shoulder <- 100
   roi.small.shared <- list(chrom="chr5", start=rsid.loc-shoulder, end=rsid.loc+shoulder)
   tbl.igap <- getVariants(mef2c, "IGAP.snpChip",   roi.small.shared)   # 1 11
   tbl.adni <- getVariants(mef2c, "ADNI.WGS",       roi.small.shared)   # 3  5
   tbl.eqtl <- getVariants(mef2c, "MAYO.eqtl.snps", roi.small.shared)   # 1 18

   checkEquals(nrow(tbl.igap), 1)
   checkEquals(nrow(tbl.adni), 3)
   checkEquals(nrow(tbl.eqtl), 1)

} # test_getVariants
#------------------------------------------------------------------------------------------------------------------------
test_getMotifs <- function()
{
   printf("--- test_getMotifs")
   roi <- list(chrom="chr5", start=88881285, end=88885739)

   tbl.motifs <- getMotifs(mef2c, source.name="allDNA-jaspar2018-human-mouse-motifs", roi, 0.95)
   checkTrue(nrow(tbl.motifs) > 300)
   checkTrue(ncol(tbl.motifs) > 10)
   checkEquals(colnames(tbl.motifs)[1:5], c("chrom", "start", "end", "name", "score"))

} # test_getMotifs
#----------------------------------------------------------------------------------------------------
test_buildModels <- function(viz=FALSE)
{
   printf("--- test_buildModels")
   tss <- mef2c@misc.data[["TSS"]]
   #tss <- 88884466
   roi.string <- sprintf("chr5:%d-%d", tss-20000, tss+20000)
   full.study.roi <- getGenomicBounds(mef2c)
   full.study.roi.string <- with(full.study.roi, sprintf("%s:%d-%d", chrom, start, end))

   x <- makeModelForRegion(mef2c, "mtx.tcx", roi.string, maxTfsInModel=15, orderBy="rfScore")
   x2 <- makeModelForRegion(mef2c, "mtx.tcx", "chr5:88,883,421-88,883,900")
   x3 <- makeModelForRegion(mef2c, "mtx.tcx", full.study.roi.string, maxTfsInModel=20, orderBy="rfScore")

   if(exists("viz")){
      showGenomicRegion(tv, full.study.roi.string)
      g <- buildMultiModelGraph("MEF2C", x3)
      g.lo <- addGeneModelLayout(g, xPos.span=1500)
      setGraph(tv, g.lo, names(x3))
      setStyle(tv, system.file(package="FRD3.data", "extdata", "style.js"))
      }

   # motifTrackForTF(mef2c, x3$motifDb$regions, "EGR3", tv)
   # motifTrackForTF(mef2c, x3$tfClass$regions, "EGR3", tv)

   #tbl.egr3 <- subset(x3$tfClassFp$regions, geneSymbol=="EGR3")
   #mm <- MotifMatcher("hg38")
   #tbl.egr3.seq <- do.call(rbind, lapply(seq_len(nrow(tbl.egr3)), function(r) getSequence(mm, tbl[r,])))
   #tbl.egr3.scores <- score.motifs(tbl.egr3)
   #seq <- getSequence(mm, data.frame(chrom="chr5", start=88661291, end=88661305, stringsAsFactors=FALSE))$seq
      # [1] "CACCGCCCACGCCAA"


} # test_buildModels
#----------------------------------------------------------------------------------------------------
test_intersectRegions_findVariantsInModel <- function(viz=FALSE)
{
   printf("--- test_intersectRegions_findVariantsInModel")
   tbl.enhancers <- mef2c@misc.data$enhancer.locs
   study.roi <- list(chrom="chr5", start=min(tbl.enhancers$start) - 10000, end=max(tbl.enhancers$end) + 10000)
   study.roi.string <- with(study.roi, sprintf("%s:%d-%d", chrom, start, end))
   tbl.enhancers <- tbl.enhancers[order(tbl.enhancers[, 2], decreasing=FALSE),]
   if(viz){
      showGenomicRegion(tv, study.roi.string)
      addBedTrackFromDataFrame(tv, "enhancers", tbl.enhancers, color="black")
      }

   tbl.fp <- getFootprints(mef2c, study.roi)

     # MotifDb::associateTranscriptionFactors needs a shortMotif column, with values like "MA0144.2"
     # extracted from the MotifDb long form, "Hsapiens-jaspar2016-STAT3-MA0144.2"
     # add them with the help of this little function

   getLastToken <- function(string){
      tokens <- strsplit(string, "-")[[1]]; return(tokens[length(tokens)])
      }

   tbl.fp$shortMotif <- unlist(lapply(tbl.fp$motifName, getLastToken))
   tbl.fpInEnhancers <- intersectRegions(mef2c, tbl.fp, tbl.enhancers)
   tbl.fpInEnhancers <- tbl.fpInEnhancers[order(tbl.fpInEnhancers$start, decreasing=FALSE),]  # 8259 x 14

      # be inclusive in the motif/TF mapping, using both methods.
      # we will call associateTranscriptionFactors[withMotif] adds a geneSymbol column, so hide the old one
      # note that we don't keep the resulting tables, just the tfs to use as candidates for trena

   colnames(tbl.fpInEnhancers)[grep("geneSymbol", colnames(tbl.fpInEnhancers))] <- "geneSymbol.orig"
   tbl.fp.motifDb.tfs <- associateTranscriptionFactors(MotifDb, tbl.fpInEnhancers, source="MotifDb", expand.rows=TRUE) # 8259
   tbl.fp.tfClass.tfs <- associateTranscriptionFactors(MotifDb, tbl.fpInEnhancers, source="TFClass", expand.rows=TRUE) # 36292
   candidate.tfs <- na.omit(unique(c(tbl.fp.motifDb.tfs$geneSymbol, tbl.fp.tfClass.tfs$geneSymbol)))  # 713

   solverNames <- c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman")
   mtx.tcx <- getExpressionMatrices(mef2c)$mtx.tcx
   solver <- EnsembleSolver(mtx.tcx, targetGene = "MEF2C", candidateRegulators = candidate.tfs, solverNames, geneCutoff = 0.5)
   tbl.model <- run(solver)

   tbl.model <- tbl.model[order(tbl.model$rfScore, decreasing=TRUE),]
   tbl.model <- subset(tbl.model, rfScore >= 1)
   tfs <- tbl.model$gene

   tbl.fpInEnhancersInModel.tfClass.tfs <- subset(tbl.fp.tfClass.tfs, geneSymbol %in% tfs)   # 781
   tbl.fpInEnhancersInModel.motifDb.tfs <- subset(tbl.fp.motifDb.tfs, geneSymbol %in% tfs)   # 194
   tbl.fpInEnhancerInModel <- rbind(tbl.fpInEnhancersInModel.tfClass.tfs,
                                    tbl.fpInEnhancersInModel.motifDb.tfs)  # 975 x 16.  should be 1901

   tbl.igapVariants <- mef2c@misc.data[["IGAP.snpChip"]]
   dim(tbl.igapVariants)

       # add a new score, -log10(pval), saving the old one
   colnames(tbl.igapVariants)[grep("score", colnames(tbl.igapVariants))] <- "oldScore"
   tbl.igapVariants$score <- -log10(tbl.igapVariants$pval)

   tbl.snpsInFpInEnhancerInModel <- intersectRegions(mef2c, tbl.igapVariants, tbl.fpInEnhancerInModel)
   tbl.snpsInFpInEnhancerInModel <- tbl.snpsInFpInEnhancerInModel[order(tbl.snpsInFpInEnhancerInModel$start, decreasing=FALSE),]

   if(viz)
      addBedGraphTrackFromDataFrame(tv, "Isnp",
                                    tbl.snpsInFpInEnhancerInModel[, c("chrom", "start", "end", "id", "score")],
                                    color="red")

   tbl.snps.final <- intersectRegions(mef2c, tbl.igapVariants, tbl.fpInEnhancerInModel)
   tbl.fp.final <- intersectRegions(mef2c, tbl.fpInEnhancerInModel, tbl.igapVariants)
   checkEquals(tbl.snps.final$start, 88883758)
   checkEquals(unique(tbl.fp.final$geneSymbol), "EGR3")

} # test_intersectRegions_findVariantsInModel
#----------------------------------------------------------------------------------------------------
test_intersectRegions_findVariantsCloseToModelTf <- function(viz=FALSE)
{
   printf("--- test_intersectRegions_findVariantsCloseToModelTf")
   tbl.enhancers <- mef2c@misc.data$enhancer.locs
   study.roi <- list(chrom="chr5", start=min(tbl.enhancers$start) - 10000, end=max(tbl.enhancers$end) + 10000)
   study.roi.string <- with(study.roi, sprintf("%s:%d-%d", chrom, start, end))
   tbl.enhancers <- tbl.enhancers[order(tbl.enhancers[, 2], decreasing=FALSE),]
   if(viz){
      showGenomicRegion(tv, study.roi.string)
      addBedTrackFromDataFrame(tv, "enhancers", tbl.enhancers, color="black")
      }

   tbl.fp <- getFootprints(mef2c, study.roi)

     # MotifDb::associateTranscriptionFactors needs a shortMotif column, with values like "MA0144.2"
     # extracted from the MotifDb long form, "Hsapiens-jaspar2016-STAT3-MA0144.2"
     # add them with the help of this little function

   getLastToken <- function(string){
      tokens <- strsplit(string, "-")[[1]]; return(tokens[length(tokens)])
      }

   tbl.fp$shortMotif <- unlist(lapply(tbl.fp$motifName, getLastToken))
   tbl.fpInEnhancers <- intersectRegions(mef2c, tbl.fp, tbl.enhancers)
   tbl.fpInEnhancers <- tbl.fpInEnhancers[order(tbl.fpInEnhancers$start, decreasing=FALSE),]  # 8259 x 14

      # be inclusive in the motif/TF mapping, using both methods.

#------------------------------------------------------------------------------------------------------------------------
if(!interactive())
   runTests()

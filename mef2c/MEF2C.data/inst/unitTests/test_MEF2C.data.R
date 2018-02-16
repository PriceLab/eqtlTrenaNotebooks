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
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_getGenomicBounds()
   test_getExpressionMatrices()
   test_getFootprints()
   test_getModels()
   test_getVariants()
   test_getMotifs()

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
    checkEquals(dim(tbl.fp), c(29397, 12))

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
if(!interactive())
   runTests()

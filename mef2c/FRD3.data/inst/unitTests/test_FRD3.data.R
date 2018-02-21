library(FRD3.data)
library(RUnit)
library(GenomicRanges)
#----------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#----------------------------------------------------------------------------------------------------
if(!exists("frd3")){
   frd3 <- FRD3.data()    # just load once, speeding up the tests
   checkTrue(all(c("SingleGeneData", "FRD3.data") %in% is(frd3)))
   }
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_getGenomicBounds()
   test_getExpressionMatrices()
   test_getMotifs()
   #test_getFootprints()
   #test_getModels()
   #test_getVariants()

} # runTests
#------------------------------------------------------------------------------------------------------------------------
test_getGenomicBounds <- function()
{
   printf("--- test_getGenomicBounds")

   roi <- getGenomicBounds(frd3)
   with(roi, checkEquals(roi, list(chrom="chr3", start=2566277, end=2572151)))

   roi.string <- getGenomicBounds(frd3, asString=TRUE)
   checkEquals(roi.string, "chr3:2566277-2572151")



} # test_getGenomicBounds
#------------------------------------------------------------------------------------------------------------------------
test_getExpressionMatrices <- function()
{
   printf("--- test_getExpressionMatrices")

   mtxs <- getExpressionMatrices(frd3)
   checkTrue("varying.zinc" %in% names(mtxs))
   mtx <- mtxs[["varying.zinc"]]
   checkEquals(dim(mtx), c(22810, 42))

} # test_getExpressionMatrices
#------------------------------------------------------------------------------------------------------------------------
test_getMotifs <- function()
{
   printf("--- test_getMotifs")

   motif.set <- "jaspar2018.athaliana.Chr3.2566277.2572151.ge85"
   checkTrue(motif.set %in% names(frd3@misc.data$motifs))        # violate privacy - fix this.
   roi <- getGenomicBounds(frd3)

   tbl.motifs <- getMotifs(frd3, motif.set, roi)
   checkTrue(nrow(tbl.motifs)  > 2500)
   checkTrue(all(tbl.motifs$score >= 0.85))
   checkEquals(colnames(tbl.motifs)[1:5], c("chrom", "start", "end", "name", "score"))

      # now, only the perfect matches
   tbl.motifs <- getMotifs(frd3, motif.set, roi, score.threshold=1.0)
   checkTrue(nrow(tbl.motifs) < 40)
   checkTrue(all(tbl.motifs$score > 0.99))
   checkEquals(colnames(tbl.motifs)[1:5], c("chrom", "start", "end", "name", "score"))

     # now narrow the roi
   span <- 1 + roi$end - roi$start
   quarter.span <- as.integer(span/4)
   roi$start <- roi$start + quarter.span
   roi$end   <- roi$end   - quarter.span
   tbl.motifs <- getMotifs(frd3, motif.set, roi, score.threshold=1.0)
   checkEquals(nrow(tbl.motifs), 14)

} # test_getExpressionMatrices
#------------------------------------------------------------------------------------------------------------------------
if(!interactive())
   runTests()

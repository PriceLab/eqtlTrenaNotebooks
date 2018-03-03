library(trena)
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

PORTS <- 5548:5599
trenaViz <- NA
if(exists("test.viz")){
   if(test.viz){
      library(trenaViz)
      trenaViz <- trenaViz(portRange=PORTS)
      setGenome(trenaViz, "tair10")
      Sys.sleep(5)
      roi <- getGenomicBounds(frd3)
      roi.string <- with(roi, sprintf("%s:%d-%d", chrom, start, end))
      showGenomicRegion(trenaViz, roi.string)   # about 3kb upstream and downstream of primary tss
      }
   }
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_getGenomicBounds()
   test_getExpressionMatrices()
   test_getMotifs()
   test_getDHS()
   test_findDHSpeaks()
   test_makeModelForRegion()
   test_motifTrackForTF()

} # runTests
#------------------------------------------------------------------------------------------------------------------------
test_getGenomicBounds <- function()
{
   printf("--- test_getGenomicBounds")

   roi <- getGenomicBounds(frd3)
   with(roi, checkEquals(roi, list(chrom="chr3", start=2566277, end=2572151)))

   roi.string <- getGenomicBounds(frd3, asString=TRUE)
   checkEquals(roi.string, "chr3:2566277-2572151")

       # these will be come first-class slots in the SingleGeneData class
   checkEquals(frd3@misc.data$targetGene, "AT3G08040")  # an orf name, the standard for athaliana
   checkEquals(frd3@misc.data$TSS, list(primary=2569502,  secondary=2572149))

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
test_getDHS <- function()
{
   printf("--- test_getDHS")
   tbls.list <- frd3@misc.data$dhs
   checkEquals(sort(names(tbls.list)), c("buds", "leaves"))
   checkEquals(dim(tbls.list$buds), c(23500, 4))
   checkEquals(dim(tbls.list$leaves), c(23500, 4))

} # test_getDHS
#------------------------------------------------------------------------------------------------------------------------
test_findDHSpeaks <- function()
{
   printf("--- test_findDHSPeaks")

   x <- sin(0:20)
   threshold <- 0.6
   minSpan <- 2
   peak.info <- findRegionsAboveThreshold(x, threshold, minSpan)
   checkEquals(peak.info$starts, c(2,8, 15))
   checkEquals(peak.info$ends, c(3, 9, 16))
   checkTrue(all(x[2:3] >= threshold))
   checkTrue(all(x[8:9] >= threshold))
   checkTrue(all(x[15:16] >= threshold))

   x <- sin(0:20)
   threshold <- 0.2
   peak.info <- findRegionsAboveThreshold(x, threshold, 3)
   checkEquals(peak.info$starts, c(8, 14))
   checkEquals(peak.info$ends, c(10, 16))
   checkTrue(all(x[8:10] >= threshold))
   checkTrue(all(x[14:16] >= threshold))

   tbl.frd3buds  <- frd3@misc.data$dhs$buds

   minSpan <- 5
   peaks.8 <- findRegionsAboveThreshold(tbl.frd3buds$score, 8, minSpan)
   checkEquals(peaks.8$starts, 9229)
   checkEquals(peaks.8$ends,   9310)
   threshold <- 3
   peaks.3 <- findRegionsAboveThreshold(tbl.frd3buds$score, threshold, minSpan)
   checkEquals(length(peaks.3$starts), 11)
   checkEquals(length(peaks.3$ends), 11)
      # all at least of minimum length?
   checkTrue(all(1 + peaks.3$ends - peaks.3$starts > minSpan))
   a <- peaks.3$starts[8]
   b <- peaks.3$ends[8]
      # check one region.  all above threshold?
   checkTrue(all(tbl.frd3buds$score[a:b] >= threshold))

} # test_findDHSpeaks
#------------------------------------------------------------------------------------------------------------------------
test_makeModelForRegion <- function(trenaViz=NA)
{
   printf("--- test_makeModelForRegion")

   roi <- "3:2,569,288-2,572,388"   # 3100 bp region, classical proximal promoter of main transcript
   x <- makeModelForRegion(frd3, dhs.cutoff=0.5, region=roi, trenaViz=trenaViz)
   checkEquals(sort(names(x)), c("model", "regions"))

   checkEquals(dim(x$model), c(10, 10))

   checkEquals(ncol(x$regions), 16)

     # todo :usually 30, 32, 34 motif/rgeions.  saw 26 once.
     # todo: not sure why the number of regions varies on each identical call.
     # wherein lies the randomness?

   checkTrue(all(x$model$gene %in% x$regions$geneSymbol))
   checkTrue(all(x$regions$geneSymbol %in% x$model$gene))
   checkTrue(all(x$regions$distance.from.tss < 3100))


} # test_makeModelForRegion
#------------------------------------------------------------------------------------------------------------------------
test_motifTrackForTF <- function(trenaViz=NA)
{
   printf("--- test_motifTrackForTF")
      # first make a model
   roi <- "3:2,569,288-2,572,388"   # 3100 bp region, classical proximal promoter of main transcript
   x <- makeModelForRegion(frd3, dhs.cutoff=1.5, region=roi, trenaViz=trenaViz)
   tbl.model <- x$model
   tbl.motifs <- x$regions
      # choose a tf with multiple binding sites
   max.bindingSites <- max(tbl.model$bindingSites)
   tf.with.max.bindingSites <- which(tbl.model$bindingSites == max.bindingSites)[1]
   tf <- tbl.model$gene[tf.with.max.bindingSites]
   tbl.bed <- motifTrackForTF(frd3, tbl.motifs, tf, trenaViz)
   checkEquals(dim(tbl.bed), c(4, 5))
   checkEquals(lapply(tbl.bed, class),
               list(chrom="character", start="numeric", end="numeric", name="character", score="numeric"))
   checkTrue(all(tbl.bed$motif == "MA0981.1"))

   checkEquals(nrow(motifTrackForTF(frd3, tbl.motifs, "intentional error")), 0)

   tf.withOneBindingSite <- subset(tbl.model, bindingSites==1)$gene[1]
   tbl.bed <- motifTrackForTF(frd3, tbl.motifs, tf.withOneBindingSite)
   checkEquals(dim(tbl.bed), c(1, 5))
   checkEquals(lapply(tbl.bed, class),
               list(chrom="character", start="numeric", end="numeric", name="character", score="numeric"))

} # test_motifTrackForTF
#------------------------------------------------------------------------------------------------------------------------
if(!interactive())
   runTests()

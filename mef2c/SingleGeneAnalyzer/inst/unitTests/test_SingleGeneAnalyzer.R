library(MEF2C.data)
library(SingleGeneAnalyzer)
library(RUnit)
#------------------------------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#------------------------------------------------------------------------------------------------------------------------
if(!exists("sga")){
   MEF2C.data <- MEF2C.data()
   sga <- SingleGeneAnalyzer(genomeName="hg38", targetGene="MEF2C", targetGene.TSS=88904257, MEF2C.data)
   checkEquals(is(sga), "SingleGeneAnalyzer")
   }
#------------------------------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_dataFrameToPandasFriendlyList()
   test_summarizeExpressionMatrices()
   test_getFootprintsForRegion()
   test_getVariantsForRegion()
   test_getDHSForRegion()
   test_getEnhancersForRegion()

} # runTests
#------------------------------------------------------------------------------------------------------------------------
test_dataFrameToPandasFriendlyList <- function()
{
   printf("--- test_dataFrameToPandasFriendlyList")

   row.count <- 5

   tbl <- data.frame(lower.case=letters[1:row.count],
                     upper.case=LETTERS[1:row.count],
                     int=sample(1:100, size=row.count),
                     float=runif(row.count),
                     stringsAsFactors=FALSE)
   rownames(tbl) <- colors()[1:row.count]


   x <- dataFrameToPandasFriendlyList(tbl)
   checkEquals(sort(names(x)), c("colnames", "rownames", "tbl"))
   checkEquals(x$colnames, colnames(tbl))
   checkEquals(x$rownames, rownames(tbl))
   checkTrue(is.null(colnames(x$tbl)))
   checkEquals(rownames(x$tbl), as.character(1:row.count))

} # test_dataFrameToPandasFriendlyList
#------------------------------------------------------------------------------------------------------------------------
test_summarizeExpressionMatrices <- function()
{
    printf("--- test_summarizeExpressionMatrices")

    tbl.summary <- summarizeExpressionMatrices(sga)
    checkEquals(class(tbl.summary), "data.frame")
    checkEquals(dim(tbl.summary), c (3, 7))
    checkEquals(colnames(tbl.summary), c("nrow", "ncol", "min", "q1", "median", "q3", "max"))
    checkTrue(all(tbl.summary$min < 0))
    checkTrue(all(tbl.summary$max > 0))

} # test_summarizeExpressionMatrices
#------------------------------------------------------------------------------------------------------------------------
test_getFootprintsForRegion <- function()
{
    printf("--- test_getFootprintsForRegion")

    roiString <- "chr5:88,883,173-88,884,172"
    tbl.fp <- getFootprintsForRegion(sga, roiString)  # uses default score.threshold of 10
    checkEquals(dim(tbl.fp), c(216, 12))
       # make sure no filtering has taken place
    checkTrue(any(tbl.fp$score < 0))
    checkTrue(any(tbl.fp$score > 0))

       # use a strong filter
    tbl.fp <- getFootprintsForRegion(sga, roiString, score.threshold=20.0)
    checkEquals(dim(tbl.fp), c(17, 12))
    checkTrue(all(tbl.fp$score >= 20.0))

} # test_getFootprintsForRegion
#------------------------------------------------------------------------------------------------------------------------
test_getVariantsForRegion <- function()
{
    printf("--- test_getVariantsForRegion")

    roiString <- "chr5:88727837-88940643"

    tbl.snp <- getVariantsForRegion(sga, roiString)  # no filtering on score (-log10(CER_P))
    checkEquals(dim(tbl.snp), c(94, 18))

       # use a strong filter
    tbl.snp <- getVariantsForRegion(sga, roiString, score.threshold=5)
    checkEquals(dim(tbl.snp), c(12, 18))
    checkTrue(all(-log10(tbl.snp$CER_P) >= 5))

      # a region with no variants
    tbl.snp <- getVariantsForRegion(sga, "chr5:88,883,347-88,884,158")
    checkEquals(nrow(tbl.snp), 0)

} # test_getVariantsForRegion
#------------------------------------------------------------------------------------------------------------------------
test_getDHSForRegion <- function()
{
    printf("--- test_getDHSForRegion")

    roiString <- "chr5:88727837-88940643"
    roiString <- "chr5:88,896,922-88,909,299"

    tbl.dhs <- getDHSForRegion(sga, roiString)  # no filtering on score, which ranges 0:1000
    checkEquals(dim(tbl.dhs), c(4, 5))
    checkEquals(colnames(tbl.dhs), c("chrom", "start", "end", "count", "score"))

       # use a strong filter
    tbl.dhs <- getDHSForRegion(sga, roiString, score.threshold=500)
    checkEquals(dim(tbl.dhs), c(1, 5))
    checkTrue(all(tbl.dhs$score >= 500))

} # test_getDHSForRegion
#------------------------------------------------------------------------------------------------------------------------
test_getEnhancersForRegion <- function()
{
    printf("--- test_getEnhancersForRegion")

    roiString <- "chr5:88727837-88940643"

    tbl <- getEnhancersForRegion(sga, roiString)
    checkEquals(dim(tbl), c(12, 3))
    checkEquals(colnames(tbl), c("chrom", "start", "end"))

    roiString <- "chr5:88727837-88727847"   # just 10 bases, expect 0 enhancers
    tbl <- getEnhancersForRegion(sga, roiString)
    checkEquals(nrow(tbl), 0)

} # test_getEnhancersForRegion
#------------------------------------------------------------------------------------------------------------------------
test_findVariantsInModelForRegion <- function()
{
   printf("--- findVariantsInModelForRegion")
   roiSting <- "chr5:88,821,191-88,821,221"    # includes rs244761, first intro variant, Chromosome: 5:88821210
                                               #    GCTAATAATTGAATATCTTTTCTTT[A/T]TTTATATATAGTTGCAGCTACAGTG
   pfms <- as.list(query(query(MotifDb, "jaspar2018"), "sapiens"))
   tbl <- findVariantsInModelForRegion(sga, roiString, 0, pfms)

} # test_findVariantsInModel
#------------------------------------------------------------------------------------------------------------------------
if(!interactive())
   runTests()


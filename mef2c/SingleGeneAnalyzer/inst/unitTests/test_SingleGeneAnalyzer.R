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
   test_summarizeExpressionMatrices()
   test_getFootprints()

} # runTests
#------------------------------------------------------------------------------------------------------------------------
test_summarizeExpressionMatrices <- function()
{
    printf("--- test_summarizeExpressionMatrices")

    pandas.friendly.list <- summarizeExpressionMatrices(sga)
    checkTrue(all(c("rownames", "colnames","tbl") %in% names(pandas.friendly.list)))
    tbl.summary <- pandas.friendly.list$tbl
    checkEquals(class(tbl.summary), "data.frame")
    checkEquals(nrow(tbl.summary), length(pandas.friendly.list$rownames))

} # test_summarizeExpressionMatrices
#------------------------------------------------------------------------------------------------------------------------
test_getFootprints <- function()
{
   printf("--- test_getFootprints")

   bounds <- getGenomicBounds(MEF2C.data)
   tbl.fp <-

} # test_getFootprints
#------------------------------------------------------------------------------------------------------------------------

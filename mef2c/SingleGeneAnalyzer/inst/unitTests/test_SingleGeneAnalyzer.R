library(SingleGeneAnalyzer)
library(RUnit)
#----------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_basicConstructor()
   test_summarizeExpressionMatrices()

} # runTests
#------------------------------------------------------------------------------------------------------------------------
test_basicConstructor <- function()
{
    printf("--- test_basicConstructor")

    ma <- SingleGeneAnalyzer(genomeName="hg38", targetGene="MEF2C", targetGene.TSS=88904257,
                             "MEF2C.data")
    checkEquals(is(ma), "SingleGeneAnalyzer")

} # test_basicConstructor
#----------------------------------------------------------------------------------------------------
test_summarizeExpressionMatrices <- function()
{
    printf("--- test_summarizeExpressionMatrices")

    ma <- SingleGeneAnalyzer(genomeName="hg38", targetGene="MEF2C", targetGene.TSS=88904257,
                             "MEF2C.data")
    pandas.friendly.list <- summarizeExpressionMatrices(ma)
    checkTrue(all(c("rownames", "colnames","tbl") %in% names(pandas.friendly.list)))
    tbl.summary <- pandas.friendly.list$tbl
    checkEquals(class(tbl.summary), "data.frame")
    checkEquals(nrow(tbl.summary), length(pandas.friendly.list$rownames))

} # test_basicConstructor
#----------------------------------------------------------------------------------------------------

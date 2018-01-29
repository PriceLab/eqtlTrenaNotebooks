library(MEF2C.data)
library(RUnit)
#----------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_basicConstructor()
   test_getGenomicBounds()
   test_getExpressionMatrices()

} # runTests
#------------------------------------------------------------------------------------------------------------------------
test_basicConstructor <- function()
{
    printf("--- test_basicConstructor")

    gd <- MEF2C.data()
    checkTrue("MEF2C.data" %in% is(gd))

} # test_basicConstructor
#------------------------------------------------------------------------------------------------------------------------
test_getGenomicBounds <- function()
{
    printf("--- test_getGenomicBounds")

    gd <- MEF2C.data()
    roi <- getGenomicBounds(gd)
    with(roi, checkEquals(roi, list(chrom="chr5", start=88618415, end=89052071)))

} # test_getGenomicBounds
#------------------------------------------------------------------------------------------------------------------------
test_getExpressionMatrices <- function()
{
    printf("--- test_getExpressionMatrices")

    gd <- MEF2C.data()
    x <- getExpressionMatrices(gd)
    mtx.names <- c("mtx.cer", "mtx.ros", "mtx.tcx")
    checkEquals(sort(names(x)), mtx.names)
    dims <- lapply(x, dim)
    checkTrue(all(lapply(dims, length) == 2))

} # test_getExpressionMatrices
#----------------------------------------------------------------------------------------------------
test_getFootprints <- function()
{
    printf("--- test_getFootprints")

    gd <- MEF2C.data()
    tbl.fp <- getFootprints(gd)

} # test_getPrecalculatedFootprints
#----------------------------------------------------------------------------------------------------

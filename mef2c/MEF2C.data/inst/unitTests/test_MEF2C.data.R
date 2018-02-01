library(MEF2C.data)
library(RUnit)
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
    checkEquals(dim(tbl.fp), c(13712, 12))

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
if(!interactive())
   runTests()

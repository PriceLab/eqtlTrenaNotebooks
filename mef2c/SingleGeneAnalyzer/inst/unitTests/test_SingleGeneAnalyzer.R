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
   test_getRegulatoryModel()
   test_getFootprintsForRegion()
   test_getVariantsForRegion()
   test_getWholeGenomeVariantsForRegion()
   test_getDHSForRegion()
   test_getEnhancersForRegion()
   test_findMotifsInRegion()
   test_findVariantsInModelForRegion()

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
test_getRegulatoryModel <- function()
{
   printf("--- test_getRegulatoryModel")
   names <- getRegulatoryModelNames(sga)
   checkTrue(length(names) > 0)
   checkTrue(nchar(names[1]) > 0)

   tbl <- getRegulatoryModel(sga, names[1])
   checkEquals(class(tbl), "data.frame")
   checkTrue(nrow(tbl) > 10)
   checkTrue(ncol(tbl) > 10)

   tbl.empty <- getRegulatoryModel(sga, "bogus")
   checkEquals(tbl.empty, data.frame())

} # test_getRegulatoryModel
#------------------------------------------------------------------------------------------------------------------------
test_getFootprintsForRegion <- function()
{
    printf("--- test_getFootprintsForRegion")

    roi.string <- "chr5:88,883,173-88,884,172"
    tbl.fp <- getFootprintsForRegion(sga, roi.string)  # uses default score.threshold of 10
    checkEquals(dim(tbl.fp), c(869, 12))
       # make sure no filtering has taken place
    checkTrue(any(tbl.fp$score < 0))
    checkTrue(any(tbl.fp$score > 0))

       # use a strong filter
    tbl.fp <- getFootprintsForRegion(sga, roi.string, score.threshold=20.0)
    checkEquals(dim(tbl.fp), c(25, 12))
    checkTrue(all(tbl.fp$score >= 20.0))

} # test_getFootprintsForRegion
#------------------------------------------------------------------------------------------------------------------------
test_getVariantsForRegion <- function()
{
    printf("--- test_getVariantsForRegion")

    roi.string <- "chr5:88727837-88940643"

    tbl.snp <- getVariantsForRegion(sga, roi.string)  # no filtering on score (-log10(CER_P))
    checkEquals(dim(tbl.snp), c(94, 18))

       # use a strong filter
    tbl.snp <- getVariantsForRegion(sga, roi.string, score.threshold=5)
    checkEquals(dim(tbl.snp), c(12, 18))
    checkTrue(all(-log10(tbl.snp$CER_P) >= 5))

      # a region with no variants
    tbl.snp <- getVariantsForRegion(sga, "chr5:88,883,347-88,884,158")
    checkEquals(nrow(tbl.snp), 0)

} # test_getVariantsForRegion
#------------------------------------------------------------------------------------------------------------------------
test_getWholeGenomeVariantsForRegion <- function()
{
    printf("--- test_getWholeGenomeVariantsForRegion")

    roi.string <- "chr5:88727837-88940643"

       # use strong filters
    tbl.var <- getWholeGenomeVariantsForRegion(sga, roi.string, altToRefRatio=2.5, minAltCount=10)
    checkEquals(dim(tbl.var), c(2, 5))

       # use weaker filter
    tbl.var <- getWholeGenomeVariantsForRegion(sga, roi.string, altToRefRatio=1.5, minAltCount=5)
    checkEquals(dim(tbl.var), c(46, 5))
    checkTrue(all(tbl.var$score >= 1.5))

      # a region and scores with no variants
    tbl.var <- getWholeGenomeVariantsForRegion(sga, "chr5:88,883,347-88,884,158", 2, 20)
    checkEquals(nrow(tbl.var), 0)

      # test a region - "chr1:1-248,956,422" - the default when a trena jupyter notebook starts up -
      # for which there is no data (in the case of mef2c)
    tbl.var <- getWholeGenomeVariantsForRegion(sga, "chr1:1-248,956,422", 0, 1)
    checkEquals(nrow(tbl.var), 0)


} # test_getVariantsForRegion
#------------------------------------------------------------------------------------------------------------------------
test_getDHSForRegion <- function()
{
    printf("--- test_getDHSForRegion")

    roi.string <- "chr5:88727837-88940643"
    roi.string <- "chr5:88,896,922-88,909,299"

    tbl.dhs <- getDHSForRegion(sga, roi.string)  # no filtering on score, which ranges 0:1000
    checkEquals(dim(tbl.dhs), c(4, 5))
    checkEquals(colnames(tbl.dhs), c("chrom", "start", "end", "count", "score"))

       # use a strong filter
    tbl.dhs <- getDHSForRegion(sga, roi.string, score.threshold=500)
    checkEquals(dim(tbl.dhs), c(1, 5))
    checkTrue(all(tbl.dhs$score >= 500))

} # test_getDHSForRegion
#------------------------------------------------------------------------------------------------------------------------
test_getEnhancersForRegion <- function()
{
    printf("--- test_getEnhancersForRegion")

    roi.string <- "chr5:88727837-88940643"

    tbl <- getEnhancersForRegion(sga, roi.string)
    checkEquals(dim(tbl), c(12, 3))
    checkEquals(colnames(tbl), c("chrom", "start", "end"))

    roi.string <- "chr5:88727837-88727847"   # just 10 bases, expect 0 enhancers
    tbl <- getEnhancersForRegion(sga, roi.string)
    checkEquals(nrow(tbl), 0)

} # test_getEnhancersForRegion
#------------------------------------------------------------------------------------------------------------------------
test_findMotifsInRegion <- function()
{
   printf("--- test_findMotifsInRegion")
   roi.string <- "chr5:88391000-89322000"   # the ~1Mb study region for MEF2C
   stat4.tfClass.motifs <-  unique(grep("^MA", geneToMotif(MotifDb, "STAT4", source="tfclass", ignore.case=TRUE)$motif, value=TRUE))
       # [1] "MA0518.1" "MA0137.3" "MA0144.2" "MA0517.1" "MA0519.1" "MA0520.1"
   stat4.mdb.motifs <-  unique(grep("^MA", geneToMotif(MotifDb, "STAT4", source="motifdb", ignore.case=TRUE)$motif, value=TRUE))

   tbl.mdb <- findMotifsInRegion(sga, roi.string, stat4.mdb.motifs, 96)
   checkEquals(dim(tbl.mdb), c(6, 13))
   tbl.tfClass <- findMotifsInRegion(sga, roi.string, stat4.tfClass.motifs, 99)
   checkEquals(dim(tbl.tfClass), c(9, 13))

} # test_findMotifsInRegion
#------------------------------------------------------------------------------------------------------------------------
test_findVariantsInModelForRegion <- function()
{
   printf("--- findVariantsInModelForRegion")
   roi.string <- "chr5:88,821,191-88,821,221"    # includes rs244761, first intro variant, Chromosome: 5:88821210
                                                 #    GCTAATAATTGAATATCTTTTCTTT[A/T]TTTATATATAGTTGCAGCTACAGTG
   roi.string <- "chr5:88,882,292-88,889,457"
   tfs.from.all.models <- unique(unlist(lapply(sga@singleGeneData@models, function(tbl) tbl$gene), use.names=FALSE))

   tbl.1 <- findVariantsInModelForRegion(sga,
                                         roi.string,
                                         motif.track="footprints",
                                         variants.track="eqtl.snps",
                                         candidate.tfs=tfs.from.all.models,
                                         tfMotifMappingName="MotifDb",
                                         shoulder=0)
   checkTrue(nrow(tbl.1) > 10)
   tfs <- sort(unique(tbl.1$geneSymbol))
   checkEquals(tfs, c("SP3", "SP4", "SP8", "ZIC1", "ZNF740"))

   tbl.2 <- findVariantsInModelForRegion(sga,
                                         roi.string,
                                         motif.track="enhancer.motifs",
                                         variants.track="wgVariants",
                                         candidate.tfs=tfs.from.all.models,
                                         tfMotifMappingName="MotifDb",
                                         shoulder=0)
   tfs.2 <- sort(unique(tbl.2$geneSymbol))
   checkEquals(tfs.2, c("BARHL2", "CUX2", "DLX6", "FOXO6", "LBX2", "MEF2A", "MEF2D", "SP3", "ZNF740"))

   tbl.3 <- findVariantsInModelForRegion(sga,
                                         roi.string,
                                         motif.track="DHS.motifs",
                                         variants.track="wgVariants",
                                         candidate.tfs=tfs.from.all.models,
                                         tfMotifMappingName="MotifDb",
                                         shoulder=0)
   tfs.3 <- sort(unique(tbl.3$geneSymbol))
   checkEquals(tfs.2, c("BARHL2", "CUX2", "DLX6", "FOXO6", "LBX2", "MEF2A", "MEF2D", "SP3", "ZNF740"))

   tbl.4 <- findVariantsInModelForRegion(sga,
                                         roi.string,
                                         motif.track="DHS.motifs",
                                         variants.track="eqtl.snps",
                                         candidate.tfs=tfs.from.all.models,
                                         tfMotifMappingName="MotifDb",
                                         shoulder=0)
   tfs.4 <- sort(unique(tbl.4$geneSymbol))
   checkTrue(nrow(tbl.4) >= 2)
   checkEquals(tfs.4, c("SP3", "ZNF740"))

   tbl.5 <- findVariantsInModelForRegion(sga,
                                         roi.string,
                                         motif.track="DHS.motifs",
                                         variants.track="eqtl.snps",
                                         candidate.tfs=tfs.from.all.models,
                                         tfMotifMappingName="TFClass",
                                         shoulder=0)

      # note especially that there are more predicted tfs from TFClass mapping
   tfs.5 <- sort(unique(tbl.5$geneSymbol))
   checkTrue(nrow(tbl.5) >= 5)
   checkEquals(tfs.5, c("EGR3", "ELK4", "SP3", "STAT4", "ZNF740"))
      # with shoulder of 0, all snps should be directly within a motif
   snps.in.motifs <- unlist(with(tbl.5, lapply(seq_len(nrow(tbl.5)),
                                               function(r) start[r] %in% motifStart[r]:motifEnd[r])))
   checkTrue(all(snps.in.motifs))

      # now expand the number of eqtl.snps found by setting a generous shoulder

   tbl.6 <- findVariantsInModelForRegion(sga,
                                         roi.string,
                                         motif.track="DHS.motifs",
                                         variants.track="eqtl.snps",
                                         candidate.tfs=tfs.from.all.models,
                                         tfMotifMappingName="TFClass",
                                         shoulder=10)

      # note especially that there are more predicted tfs from TFClass mapping
   tfs.6 <- sort(unique(tbl.6$geneSymbol))
      # with shoulder of 10, not all snps should be directly within a motif
   snps.in.motifs <- unlist(with(tbl.6, lapply(seq_len(nrow(tbl.6)),
                                               function(r) start[r] %in% motifStart[r]:motifEnd[r])))
   checkTrue(any(snps.in.motifs))
   checkTrue(!all(snps.in.motifs))

} # test_findVariantsInModelForRegion
#------------------------------------------------------------------------------------------------------------------------
doofus <- function(tf="FOXP1")
{
  tbl <- findMotifsInRegion(sga, "chr5:88391000-89322000", "MA0518.1", 90)
  checkEquals(dim(tbl), c(7, 12))

  model.name <- "mef2c.cory.wgs.tcx.tfClass"
  tbl.model <- sga@singleGeneData@models[[model.name]]
  tfs <- head(tbl.model$gene, n=10)

   lookup <- function(gene){
     tbl.mdb <- subset(geneToMotif(MotifDb, gene, "MotifDb", ignore.case=TRUE), dataSource=="jaspar2018")
     motifs.mdb <- unique(grep("^MA", tbl.mdb$motif, value=TRUE))
     tbl.tfc <- subset(geneToMotif(MotifDb, gene, "TFClass", ignore.case=TRUE))
     motifs.tfc <- c() #unique(grep("^MA", tbl.tfc$motif, value=TRUE))
     list(motifdb=motifs.mdb, tfclass=motifs.tfc)
     }


   for(tf in tfs){
      printf("---------- tf: %s", tf)
      motifs <- unique(unlist(lookup(tf), use.names=FALSE))   # lowest of 5 footprints is for MEF2C
      motif.assay <- lapply(motifs, function(motif) assessSnp(trena, as.list(query(query(MotifDb, "jaspar2018"), motif)), "rs244761", 15, 80))
      names(motif.assay) <- motifs
      print(motif.assay)
      }

} # doofus
#------------------------------------------------------------------------------------------------------------------------
if(!interactive())
   runTests()

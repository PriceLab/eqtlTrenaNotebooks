library(trena)
library(trenaViz); PORT_RANGE <- 10000:10200
library(RUnit)
#------------------------------------------------------------------------------------------------------------------------
explore.enhancers <- function(targetGene="MEF2C")
{
   targetGene <- "MEF2C"
   TSS <- 88904257

   data.dir <- "~/github/dockerizedMicroservices/trena.mef2c/dev/trns/geneHancer"
   tbl.elite <- read.table(sprintf("%s/enhancer_elite_ids.txt", data.dir), sep="\t", as.is=TRUE, header=TRUE) #  243281      6
   tbl.geneScoresAll <- read.table(sprintf("%s/enhancer_gene_scores.txt", data.dir),
                                   sep="\t", as.is=TRUE, header=TRUE)  #  934287      9
   tbl.mef2cScores <- subset(tbl.geneScoresAll, symbol==targetGene)
   clusters.mef2c <- tbl.mef2cScores$cluster_id
   tbl.mef2c.eLocs <- subset(tbl.elite, cluster_id %in% clusters.mef2c)[, c("chrom", "start", "end")]
   tbl.mef2c.eLocs$chrom <- paste("chr", tbl.mef2c.eLocs$chrom, sep="")

   tbl.roi <- data.frame(chrom="chr5", start=TSS-5000, end=TSS+5000, stringsAsFactors=FALSE)

    #--------------------------------------------------------------------------------
    # use mtx.tcx, pre-calculated footprints, cory's 10kb region, TFclass mapping
    #--------------------------------------------------------------------------------

   x0 <- makeGeneModel(bindingSiteSource="footprints",
                       tbl.roi,
                       mtx=mtx.tcx,
                       tfMappingSource="TFclass",
                       targetGene="MEF2C",
                       orderByColumn="pcaMax",
                       solverInclusivenessCutoff=1.0)


   tbl.roi.enhanced <- rbind(tbl.roi, tbl.mef2c.eLocs)
   x1 <- makeGeneModel(bindingSiteSource="footprints",
                       tbl.roi.enhanced,
                       mtx=mtx.tcx,
                       tfMappingSource="TFclass",
                       targetGene="MEF2C",
                       orderByColumn="pcaMax",
                       solverInclusivenessCutoff=1.0)


   tbl.mef2c.eLocs$size <-  1 + tbl.mef2c.eLocs$end - tbl.mef2c.eLocs$start
   tbl.roi.enhanced.small <- rbind(tbl.roi, subset(tbl.mef2c.eLocs, size <500)[, 1:3])

   x2 <- makeGeneModel(bindingSiteSource="footprints",
                       tbl.roi.enhanced.small,
                       mtx=mtx.tcx,
                       tfMappingSource="TFclass",
                       targetGene="MEF2C",
                       orderByColumn="pcaMax",
                       solverInclusivenessCutoff=1.0)





} # explore.enhancers
#------------------------------------------------------------------------------------------------------------------------
getEnhancerMotifs <- function(tbl.mef2c.eLocs)
{
   pfms <- as.list(query(query(MotifDb, "jaspar2018"), "hsapiens"))
   mm <- MotifMatcher(genomeName="hg38", pfms)
   tbl.motifs <- findMatchesByChromosomalRegion(mm, tbl.mef2c.eLocs, pwmMatchMinimumAsPercentage=85)

    shortMotifNames <- unlist(lapply(tbl.motifs$motifName,
                              function(name) {tokens <- strsplit(name, "-")[[1]]; return(tokens[length(tokens)])}))
    tbl.motifs$shortMotif <- shortMotifNames
    tbl.eMotifs.tfc <- associateTranscriptionFactors(MotifDb, tbl.motifs, source="TFClass", expand.rows=TRUE)
    dim(tbl.eMotifs.tfc)  # [1] 891914     15
    tbl.eMotifs.mdb <- associateTranscriptionFactors(MotifDb, tbl.motifs, source="MotifDb", expand.rows=TRUE)
    dim(tbl.eMotifs.mdb)  # [1]  43076    15

     save(tbl.eMotifs.mdb, tbl.eMotifs.tfc, file="~/github/eqtlTrenaNotebooks/mef2c/trena/data/tbl.eMotifs.RData")


} # getEnhancerMotifs
#------------------------------------------------------------------------------------------------------------------------


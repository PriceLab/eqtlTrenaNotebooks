library(trena)
library(trenaViz)
library(FRD3.data)
#------------------------------------------------------------------------------------------------------------------------
PORTS <- 10000:10020
#------------------------------------------------------------------------------------------------------------------------
if(!exists("trena")){
   #trena <- Trena("mm10")   # a dodge, one which works if getProximalPromoter method is avoided
   tv <- trenaViz(portRange=PORTS)
   setGenome(tv, "tair10")
   frd3 <- FRD3.data()
   target.orf <- frd3@misc.data$targetGene
   setBrowserWindowTitle(tv, "FRD3")
   }
#------------------------------------------------------------------------------------------------------------------------
make.model.in.currentRegion <- function(dhs.cutoff)
{
   tss <- frd3@misc.data$TSS

   current.roi <- parseChromLocString(getGenomicRegion(tv))
   tbl.croi <- with(current.roi, data.frame(chrom=chrom, start=start, end=end, stringsAsFactors=FALSE))
   tbl.croi$chrom <- sub("chr", "Chr", tbl.croi$chrom)

   tbl.frd3buds  <- frd3@misc.data$dhs$buds
   minSpan <- 5
   peaks <- findRegionsAboveThreshold(tbl.frd3buds$score, threshold=dhs.cutoff, minSpan)

   peak.position.offset <- tbl.frd3buds$start[1]

   tbl.dhs <- with(peaks, data.frame(chrom="Chr3",
                                     start=starts+peak.position.offset,
                                     end=ends+peak.position.offset,
                                     stringsAsFactors=FALSE))

   track.title <- sprintf("bud dhs > %f", dhs.cutoff)
   addBedTrackFromDataFrame(tv, track.title, tbl.dhs, color="darkGreen", trackHeight=50)

   tbl.ov <- as.data.frame(findOverlaps(GRanges(tbl.croi), GRanges(tbl.dhs), type="any"))
   colnames(tbl.ov) <- c("croi", "dhs")
   stopifnot(nrow(tbl.ov) > 0)
   tbl.dhs.currentView <- tbl.dhs[tbl.ov$dhs,]

   pfms <- as.list(query(query(MotifDb, "jaspar2018"), "athaliana"))
   mm <- MotifMatcher("tair10", pfms)
   tbl.motifs <- findMatchesByChromosomalRegion(mm, tbl.dhs.currentView, pwmMatchMinimumAsPercentage=90)

   dim(tbl.motifs)
   uniprotIDs <- mcols(MotifDb[tbl.motifs$motifName])$proteinId
   tbl.motifs$uniprot <- uniprotIDs
   match.indices <- match(tbl.motifs$uniprot, tbl.xref$proteinId)
   tbl.motifs$tf <- tbl.xref$orfs[match.indices]

   mtx <- getExpressionMatrices(frd3)$varying.zinc

   stopifnot(target.orf %in% rownames(mtx))

   tf.candidates <- tbl.motifs$tf
   colnames(tbl.motifs)[grep("tf", colnames(tbl.motifs))] <- "geneSymbol"
   colnames(tbl.motifs)[grep("motifStart", colnames(tbl.motifs))] <- "start"
   colnames(tbl.motifs)[grep("motifEnd", colnames(tbl.motifs))] <- "end"

   solver.names <- c("lasso", "pearson", "randomForest", "ridge", "spearman")

   tbl.model <- createGeneModel(trena, target.orf, solver.names, tbl.motifs, mtx)
   match.indices <- match(tbl.model$gene, tbl.orfSym$TAIR)
   tbl.model$sym <- tbl.orfSym$SYMBOL[match.indices]
   mapping.failures <- which(is.na(tbl.model$sym))

   if(length(mapping.failures) > 0)
      tbl.model$sym[mapping.failures] <- tbl.model$gene[mapping.failures]
   tbl.model <- tbl.model[order(tbl.model$rfScore, decreasing=TRUE),]
   tbl.model.strong <- head(tbl.model, n=10)
   tbl.motifs.strong <- subset(tbl.motifs, geneSymbol %in% tbl.model.strong$gene)

   distance <- tbl.motifs.strong$start - tss$primary
   direction <- rep("upstream", length(distance))
   direction[which(distance < 0)] <- "downstream"
   tbl.motifs.strong$distance.from.tss <- distance
   tbl.motifs.strong$id <- sprintf("%s.dhs.%s.%06d.%s", target.orf, direction, abs(distance), tbl.motifs.strong$motifName)

   return(list(model=tbl.model.strong, regions=tbl.motifs.strong))

} # make.model.in.currentRegion
#------------------------------------------------------------------------------------------------------------------------
run <- function()
{
   roi <- getGenomicBounds(frd3)
   roi.string <- with(roi, sprintf("%s:%d-%d", chrom, start, end))
   showGenomicRegion(tv, roi.string)
   x1 <- make.model.in.currentRegion(dhs.cutoff=0.5)
   x2 <- make.model.in.currentRegion(1.0)

   showGenomicRegion(tv, "3:2,569,270-2,571,414")   # main transcript
   x3 <- make.model.in.currentRegion(0.5)
   showGenomicRegion(tv, "3:2,571,841-2,573,985")
   x4 <- make.model.in.currentRegion(0.5)

   models <- list(dhs.gt05=x1, dhs.gt10=x2, dhs.left=x3, dhs.right=x4)
   #models <- list(a=x1)

   g <- buildMultiModelGraph(target.orf, models)
   g.lo <- addGeneModelLayout(g, xPos.span=1500)
   setGraph(tv, g.lo, names(models))
   setStyle(tv, system.file(package="FRD3.data", "extdata", "style.js"))

} # run
#------------------------------------------------------------------------------------------------------------------------

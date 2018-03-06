library(trena)
library(trenaViz)
library(MEF2C.data)

PORTS <- 10000:10020   # trenaViz claims a websocket port in this range, over which R, igv and cytoscape communicate

tv <- trenaViz(portRange=PORTS)
setGenome(tv, "hg38")
mef2c <- MEF2C.data()
target.gene <- "MEF2C"
setBrowserWindowTitle(tv, "MEF2C")
showGenomicRegion(tv, "MEF2C")
tbl.enhancers <- mef2c@misc.data$enhancer.locs
# roi: region of interest
study.roi <- list(chrom="chr5", start=min(tbl.enhancers$start) - 10000, end=max(tbl.enhancers$end) + 10000)
study.roi.string <- with(study.roi, sprintf("%s:%d-%d", chrom, start, end))
showGenomicRegion(tv, study.roi.string)
tbl.enhancers <- tbl.enhancers[order(tbl.enhancers[, 2], decreasing=FALSE),]
addBedTrackFromDataFrame(tv, "enhancers", tbl.enhancers, color="black")

tbl.fp <- getFootprints(mef2c, study.roi)
tbl.fp <- tbl.fp[order(tbl.fp$start, decreasing=FALSE),]


tbl.ov <- as.data.frame(findOverlaps(GRanges(tbl.fp), GRanges(tbl.enhancers)))
colnames(tbl.ov) <- c("fp", "enhancers")
tbl.fpe <- tbl.fp[unique(tbl.ov$fp),]
addBedTrackFromDataFrame(tv, "fpe", tbl.fpe, color="darkGreen")



addBedTrackFromDataFrame(tv, "fp", tbl.fp, color="darkGreen")


tbl.eMotifs.mdb <- mef2c@misc.data$enhancer.motifs.mdb
name <- paste(tbl.eMotifs.mdb$shortMotif, tbl.eMotifs.mdb$geneSymbol, sep="-")
tbl.eMotifs.mdb <- mef2c@misc.data$enhancer.motifs.mdb[, c("chrom", "motifStart", "motifEnd", "motifRelativeScore")]
tbl.eMotifs.mdb$name <- name
tbl.eMotifs.mdb <- tbl.eMotifs.mdb[, c("chrom", "motifStart", "motifEnd", "name", "motifRelativeScore")]
colnames(tbl.eMotifs.mdb) <- c("chrom", "start", "end", "name", "score")
tbl.eMotifs.mdb <- tbl.eMotifs.mdb[order(tbl.eMotifs.mdb[, 2], decreasing=FALSE),]
addBedTrackFromDataFrame(tv, "eMoMDB", tbl.eMotifs.mdb, color="brown")


tbl.eMotifs.tfc <- mef2c@misc.data$enhancer.motifs.tfc
name <- paste(tbl.eMotifs.tfc$shortMotif, tbl.eMotifs.tfc$geneSymbol, sep="-")
tbl.eMotifs.tfc <- mef2c@misc.data$enhancer.motifs.tfc[, c("chrom", "motifStart", "motifEnd", "motifRelativeScore")]
tbl.eMotifs.tfc$name <- name
tbl.eMotifs.tfc <- tbl.eMotifs.tfc[, c("chrom", "motifStart", "motifEnd", "name", "motifRelativeScore")]
colnames(tbl.eMotifs.tfc) <- c("chrom", "start", "end", "name", "score")
tbl.eMotifs.tfc <- tbl.eMotifs.tfc[order(tbl.eMotifs.tfc[, 2], decreasing=FALSE),]
addBedTrackFromDataFrame(tv, "eMoTFC", tbl.eMotifs.tfc, color="brown")

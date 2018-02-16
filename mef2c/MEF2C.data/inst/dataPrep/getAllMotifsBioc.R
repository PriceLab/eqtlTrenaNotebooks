library(trena)
library(MEF2C.data)
pfms.human <- as.list(query(query(MotifDb, "jaspar2018"), "hsapiens"))
pfms.mouse <- as.list(query(query(MotifDb, "jaspar2018"), "mmusculus"))
pfms <- c(pfms.human, pfms.mouse)
mm <- MotifMatcher(genomeName="hg38", pfms)
mef2c.data <-MEF2C.data()
roi <- getGenomicBounds(mef2c.data)
tbl.roi <- as.data.frame(roi, stringsAsFactor=FALSE)
#tbl.roi$end <- tbl.roi$start + 10000
printf("region size: %d", with(tbl.roi, 1+end-start))
tbl.motifs <- findMatchesByChromosomalRegion(mm, tbl.roi, pwmMatchMinimumAsPercentage=85)

getLastToken <- function(longMotifName){
  tokens <- strsplit(longMotifName, "-")[[1]];
  return(tokens[length(tokens)])
  }

names <-unlist(lapply(tbl.motifs$motifName, getLastToken))
tbl.motifs$name <- names

colnames(tbl.motifs)[grep("motifStart", colnames(tbl.motifs))] <- "start"
colnames(tbl.motifs)[grep("motifEnd", colnames(tbl.motifs))] <- "end"
colnames(tbl.motifs)[grep("motifRelativeScore", colnames(tbl.motifs))] <- "score"
coi.bed <- c("chrom", "start", "end", "name", "score")
coi.other <- setdiff(colnames(tbl.motifs), coi.bed)
tbl.motifs <- tbl.motifs[, c(coi.bed, coi.other)]

filename <- sprintf("tbl.motifs.jaspar2018.human.mouse.%s.%d.%d.RData",  tbl.roi$chrom, tbl.roi$start, tbl.roi$end)
save(tbl.motifs, file=filename)

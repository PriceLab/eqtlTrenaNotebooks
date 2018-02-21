library(trena)
library(FRD3.data)
frd3 <- FRD3.data()

mdb <- query(query(MotifDb, "jaspar2018"), "athaliana")
pfms <- as.list(mdb)   # 452
mm <- MotifMatcher(genomeName="tair10", pfms)

roi <- getGenomicBounds(frd3)   # 5874 bp
tbl.roi <- as.data.frame(roi, stringsAsFactor=FALSE)
tbl.roi$chrom <- sub("chr", "Chr", tbl.roi$chrom)

printf("region size: %d", with(tbl.roi, 1+end-start))
tbl.motifs <- findMatchesByChromosomalRegion(mm, tbl.roi, pwmMatchMinimumAsPercentage=85)

# make this into a standard bed table format data.frame
colnames(tbl.motifs)[grep("motifStart", colnames(tbl.motifs))] <- "start"
colnames(tbl.motifs)[grep("motifEnd",   colnames(tbl.motifs))] <- "end"
colnames(tbl.motifs)[grep("motifRelativeScore", colnames(tbl.motifs))] <- "score"

getLastToken <- function(longMotifName){
  tokens <- strsplit(longMotifName, "-")[[1]];
  return(tokens[length(tokens)])
  }

names <-unlist(lapply(tbl.motifs$motifName, getLastToken))
tbl.motifs$name <- names

coi.bed <- c("chrom", "start", "end", "name", "score")
all(coi.bed %in% colnames(tbl.motifs))
coi.other <- setdiff(colnames(tbl.motifs), coi.bed)
tbl.motifs <- tbl.motifs[, c(coi.bed, coi.other)]

filename <- sprintf("../extdata/tbl.motifs.jaspar2018.athaliana.%s.%d.%d.ge85.RData",  tbl.roi$chrom, tbl.roi$start, tbl.roi$end)
save(tbl.motifs, file=filename)

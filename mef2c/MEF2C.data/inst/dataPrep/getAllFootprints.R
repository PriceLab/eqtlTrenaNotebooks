library(trena)
library(trenaViz)
library(MotifDb)
#------------------------------------------------------------------------------------------------------------------------
# takes about 4 minutes (5 feb 2018)
#
# previous version:
#   print(load("../extdata/tbl.fp.chr5.88615025-89052115.4sources.noDups.RData")) #[1] "tbl.fp"
#   dim(tbl.fp) # [1] 13712    12
# this version:
#    29397    12
#------------------------------------------------------------------------------------------------------------------------
if(!exists("trena"))
   trena <- Trena("hg38")

roi <- list(chrom='chr5', start= 88391000, end=89322000)
targetGene <- "MEF2C"
# targetGene.tss <- 88904257
targetGene.tss <- 88884466   # from igv loaded with genecode v24

source.1 <- "postgres://bddsrds.globusgenomics.org/skin_wellington_16"
source.2 <- "postgres://bddsrds.globusgenomics.org/skin_wellington_20"
source.3 <- "postgres://bddsrds.globusgenomics.org/skin_hint_16"
source.4 <- "postgres://bddsrds.globusgenomics.org/skin_hint_20"
sources <- c(source.1, source.2, source.3, source.4)
source.short.names <-  c("well_16", "well_20", "hint_16", "hint_20")
names(sources) <-source.short.names

# roi$end <- roi$start + 2000   # for a quick experiment

x <- getRegulatoryChromosomalRegions(trena, roi$chrom, roi$start, roi$end, sources, targetGene, targetGene.tss)
names(x) <- source.short.names
x2 <- x

for(source.name in source.short.names){
   tbl <- x2[[source.name]]
   if(nrow(tbl) > 0)
      tbl$db <- source.name
   x2[[source.name]] <- tbl
   }

tbl.fp0 <- do.call(rbind, x2)
tbl.fp2 <- associateTranscriptionFactors(MotifDb, tbl.fp0, source="MotifDb", expand.rows=TRUE)
rownames(tbl.fp2) <- NULL
colnames(tbl.fp2)[grep("motifStart", colnames(tbl.fp2))] <- "start"
colnames(tbl.fp2)[grep("motifEnd", colnames(tbl.fp2))] <- "end"
preferred.column.order <- c("chrom", "start", "end", "motifName", "strand", "score", "length",
                            "distance.from.tss", "id", "db", "geneSymbol", "pubmedID")
stopifnot(sort(colnames(tbl.fp2)) == sort(preferred.column.order))
tbl.fp <- tbl.fp2[, preferred.column.order]

# crude elimination of duplicates, with
dups <- which(duplicated(tbl.fp[, c("chrom", "start", "end")]))
tbl.fp <- tbl.fp[-dups,]    # 23397 x 12
filename <- sprintf("../extdata/tbl.fp.%s.%d-%d.4sources.noDups.RData", roi$chrom, roi$start, roi$end)
save(tbl.fp, file=filename)


#------------------------------------------------------------------------------------------------------------------------

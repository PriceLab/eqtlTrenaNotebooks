library(rtracklayer)
#------------------------------------------------------------------------------------------------------------------------
buds.bw <- "Ath_buds_DNase.bw"
file.exists(buds.bw)
chr3.bw <- import(buds.bw, which=seqinfo(BigWigFile(buds.bw))["Chr3"])
frd3.roi.gr <- GRanges(seqnames="Chr3", IRanges(start=2557465, end=2580964))
frd3.data <- subsetByOverlaps(chr3.bw, frd3.roi.gr)
frd3.data <- renameSeqlevels(frd3.data, tolower(seqlevels(chr3.bw)))
export(frd3.data, "frd3-buds.bw", format="bigwig")

tbl.frd3buds <- as.data.frame(frd3.data)[, c("seqnames", "start", "end", "score")]
colnames(tbl.frd3buds)[1] <- "chrom"
tbl.frd3buds$chrom <- as.character(tbl.frd3buds$chrom)

leaf.bw <- "Ath_leaf_DNase.bw"
chr3.bw <- import(leaf.bw, which=seqinfo(BigWigFile(leaf.bw))["Chr3"])
frd3.roi.gr <- GRanges(seqnames="Chr3", IRanges(start=2557465, end=2580964))
frd3.data <- subsetByOverlaps(chr3.bw, frd3.roi.gr)
frd3.data <- renameSeqlevels(frd3.data, tolower(seqlevels(chr3.bw)))
export(frd3.data, "frd3-leaf.bw", format="bigwig")

tbl.frd3leaf <- as.data.frame(frd3.data)[, c("seqnames", "start", "end", "score")]
colnames(tbl.frd3leaf)[1] <- "chrom"
tbl.frd3leaf$chrom <- as.character(tbl.frd3leaf$chrom)

save(tbl.frd3buds, tbl.frd3leaf, file="tbls.frd3.dhs.budAndLeaf.RData")



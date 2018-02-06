library(VariantAnnotation)
vcf.filename <- "~/github/projects/priceLab/cory/mef2c-vcf/SCH_11923_B01_GRM_WGS_2017-04-27_5.recalibrated_variants.vcf.gz"
file.exists(vcf.filename)

# current project bounds: chr5:88391000-89322000
# https://genome.ucsc.edu/cgi-bin/hgLiftOver
hg38.roi <- "chr5:88391000-89322000"
hg19.roi <- "chr5:87686817-88617817"
hg19.start <- 87686817
hg19.end   <- 88617817

mef2c.gr <- GRanges(5, IRanges(hg19.start, hg19.end))
params <- ScanVcfParam(which=mef2c.gr)

tabixFile <- TabixFile(vcf.filename)
vcf <- readVcf(tabixFile, "hg19", params)
vcf.geno <- geno(vcf)
mtx.geno <- vcf.geno$GT   # 1687  349
mtx.012 <- matrix(0, nrow=nrow(mtx.geno), ncol=ncol(mtx.geno), dimnames=list(rownames(mtx.geno), colnames(mtx.geno)))
mtx.012[which(mtx.geno=="0/1")] <- 1
mtx.012[which(mtx.geno=="1/1")] <- 2
mtx.geno <- mtx.012
dim(mtx.geno)  # 5824 349
sum(mtx.geno)  # [1] 299636


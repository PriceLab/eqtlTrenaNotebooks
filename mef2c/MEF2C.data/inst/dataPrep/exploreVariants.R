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

#--------------------------------------------------------------------------------
# load (among other things) sample classification
#--------------------------------------------------------------------------------
print(load("~/github/dockerizedMicroservices/trena.mef2c/trena/data/mtx.tcx.pheno.geno.RData"))
length(samples.ad)   # 79
length(samples.ctl)  # 73

#--------------------------------------------------------------------------------
# read in and clean up the variants matrix: 0 for wt/wt, 1 for snp/wt, 2 for snp/snp
#--------------------------------------------------------------------------------
vcf.geno <- geno(vcf)
mtx.geno <- vcf.geno$GT   # 1687  349



mtx.012 <- matrix(0, nrow=nrow(mtx.geno), ncol=ncol(mtx.geno), dimnames=list(rownames(mtx.geno), colnames(mtx.geno)))
mtx.012[which(mtx.geno=="0/1")] <- 1
mtx.012[which(mtx.geno=="1/1")] <- 2
mtx.geno <- mtx.012
dim(mtx.geno)   # 5824 349

#--------------------------------------------------------------------------------
# define a function to parse the rownames into chrom, pos, ref, alt
#--------------------------------------------------------------------------------
parse <- function(string, regex, numeric.fields){
   match <- regexpr(regex, string, perl=TRUE)
   if(match == 1){
      info <- attributes(match)
      result <- list()
      for(name in info$capture.names){
         start <- as.integer(info$capture.start[1, name])
         length <- as.integer(info$capture.length[1, name])
         result[[name]] <- substring(string, start, start+length-1)
         if(name %in% numeric.fields)
            result[[name]] <- as.numeric(result[[name]])
         } # for name
      } # if match
   result
   }
#--------------------------------------------------------------------------------
# call that function, create a data.frame of the results
#--------------------------------------------------------------------------------
regex <- "(?<chrom>.*):(?<pos>\\d+)_(?<ref>.*)\\/(?<alt>.*)"
x <- lapply(rownames(mtx.geno), function(rowname) parse(rowname, regex, "pos"))
tbl.pos <- as.data.frame(do.call(rbind, x))
tbl.pos$chrom <- paste("chr", tbl.pos$chrom, sep="")

#--------------------------------------------------------------------------------
# add two columns, giving the count of variant samples in AD and CONTROL
#--------------------------------------------------------------------------------
ad.alt.count <- as.numeric(apply(mtx.geno, 1, function(row) length(which(row[samples.ad] == 1))))
ctl.alt.count <- as.numeric(apply(mtx.geno, 1, function(row) length(which(row[samples.ctl] == 1))))
tbl.pos$het.altAD <- ad.alt.count
tbl.pos$het.altCTL <- ctl.alt.count

ad.alt.count <- as.numeric(apply(mtx.geno, 1, function(row) length(which(row[samples.ad] == 2))))
ctl.alt.count <- as.numeric(apply(mtx.geno, 1, function(row) length(which(row[samples.ctl] == 2))))
tbl.pos$hom.altAD <- ad.alt.count
tbl.pos$hom.altCTL <- ctl.alt.count

ad.alt.count <- as.numeric(apply(mtx.geno, 1, function(row) length(which(row[samples.ad] > 0))))
ctl.alt.count <- as.numeric(apply(mtx.geno, 1, function(row) length(which(row[samples.ctl] > 0))))
tbl.pos$any.altAD <- ad.alt.count
tbl.pos$any.altCTL <- ctl.alt.count

subset(tbl.pos, any.altAD > (2*any.altCTL) &  any.altAD > 10)[, c(1,2,5:10)]

#--------------------------------------------------------------------------------
# translate to hg38
#--------------------------------------------------------------------------------
chromLocs <- sprintf("%s:%d-%d", tbl.pos$chrom, unlist(tbl.pos$pos), unlist(tbl.pos$pos))
write.table(data.frame(chromLoc=chromLocs, stringsAsFactors=FALSE), file="tbl.pos.hg19.bed",
            sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

#--------------------------------------------------------------------------------
# do liftover, hg19 -> hg38  https://genome.ucsc.edu/cgi-bin/hgLiftOver
# when complete: e.g. curl -O https://genome.ucsc.edu/trash/hglft_genome_496c_90e830.bed
#--------------------------------------------------------------------------------
hg38.coords <- read.table("tbl.pos.hg38.bed", sep="\t", as.is=TRUE)$V1
hg38.pos <- unlist(lapply(strsplit(hg38.coords, "-"), function(tokens) tokens[length(tokens)]))
tbl.pos$pos <- as.integer(hg38.pos)
tbl.pos$ref <- as.character(tbl.pos$ref)
tbl.pos$alt <- as.character(tbl.pos$alt)

save(tbl.pos, file="tbl.vcf.chr5.88391000.89322000.79AD.73CTL.RData")
save(tbl.pos, file="../extdata/tbl.vcf.chr5.88391000.89322000.79AD.73CTL.RData")

# class: CollapsedVCF
# dim: 3 349
# rowRanges(vcf):
#   GRanges with 5 metadata columns: paramRangeID, REF, ALT, QUAL, FILTER
# info(vcf):
#   DataFrame with 24 columns: AC, AF, AN, BaseQRankSum, DP, DS, END, ExcessHet, FS, HaplotypeScore, InbreedingCoeff, MLEAC, MLEAF, MQ, MQ0, MQRankSum, NEGATIVE_TRAIN_SITE, POSITIVE_TRAIN_SITE, QD, ReadPosRankSum, SOR, VQSLOD, VariantType, culprit
# info(header(vcf)):
#                        Number Type    Description
#    AC                  A      Integer Allele count in genotypes, for each ALT allele, in the same order as listed
#    AF                  A      Float   Allele Frequency, for each ALT allele, in the same order as listed
#    AN                  1      Integer Total number of alleles in called genotypes
#    BaseQRankSum        1      Float   Z-score from Wilcoxon rank sum test of Alt Vs. Ref base qualities
#    DP                  1      Integer Approximate read depth; some reads may have been filtered
#    DS                  0      Flag    Were any of the samples downsampled?
#    END                 1      Integer Stop position of the interval
#    ExcessHet           1      Float   Phred-scaled p-value for exact test of excess heterozygosity
#    FS                  1      Float   Phred-scaled p-value using Fisher's exact test to detect strand bias
#    HaplotypeScore      1      Float   Consistency of the site with at most two segregating haplotypes
#    InbreedingCoeff     1      Float   Inbreeding coefficient as estimated from the genotype likelihoods per-sample when compared against the Hardy-Weinberg expectation
#    MLEAC               A      Integer Maximum likelihood expectation (MLE) for the allele counts (not necessarily the same as the AC), for each ALT allele, in the same order as listed
#    MLEAF               A      Float   Maximum likelihood expectation (MLE) for the allele frequency (not necessarily the same as the AF), for each ALT allele, in the same order as listed
#    MQ                  1      Float   RMS Mapping Quality
#    MQ0                 1      Integer Total Mapping Quality Zero Reads
#    MQRankSum           1      Float   Z-score From Wilcoxon rank sum test of Alt vs. Ref read mapping qualities
#    NEGATIVE_TRAIN_SITE 0      Flag    This variant was used to build the negative training set of bad variants
#    POSITIVE_TRAIN_SITE 0      Flag    This variant was used to build the positive training set of good variants
#    QD                  1      Float   Variant Confidence/Quality by Depth
#    ReadPosRankSum      1      Float   Z-score from Wilcoxon rank sum test of Alt vs. Ref read position bias
#    SOR                 1      Float   Symmetric Odds Ratio of 2x2 contingency table to detect strand bias
#    VQSLOD              1      Float   Log odds of being a true variant versus being false under the trained gaussian mixture model
#    VariantType         1      String  Variant type description
#    culprit             1      String  The annotation which was the worst performing in the Gaussian mixture model, likely the reason why the variant was filtered out
# geno(vcf):
#   SimpleList of length 12: GT, AB, AD, DP, GQ, MIN_DP, MQ0, PGT, PID, PL, RGQ, SB
# geno(header(vcf)):
#           Number Type    Description
#    GT     1      String  Genotype
#    AB     1      Float   Allele balance for each het genotype
#    AD     .      Integer Allelic depths for the ref and alt alleles in the order listed
#    DP     1      Integer Approximate read depth (reads with MQ=255 or with bad mates are filtered)
#    GQ     1      Integer Genotype Quality
#    MIN_DP 1      Integer Minimum DP observed within the GVCF block
#    MQ0    1      Integer Number of Mapping Quality Zero Reads per sample
#    PGT    1      String  Physical phasing haplotype information, describing how the alternate alleles are phased in relation to one another
#    PID    1      String  Physical phasing ID information, where each unique ID within a given sample (but not across samples) connects records within a phasing group
#    PL     G      Integer Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification
#    RGQ    1      Integer Unconditional reference genotype confidence, encoded as a phred quality -10*log10 p(genotype call is wrong)
#    SB     4      Integer Per-sample component statistics which comprise the Fisher's Exact Test to detect strand bias.
#
#
# vcf.geno <- geno(vcf)
# mtx.geno <- vcf.geno$GT   # 3  349
# mtx.012 <- matrix(0, nrow=nrow(mtx.geno), ncol=ncol(mtx.geno), dimnames=list(rownames(mtx.geno), colnames(mtx.geno)))
# mtx.012[which(mtx.geno=="0/1")] <- 1
# mtx.012[which(mtx.geno=="1/1")] <- 2
# mtx.geno <- mtx.012
# dim(mtx.geno)  # 5824 349
# sum(mtx.geno)  # [1] 299636


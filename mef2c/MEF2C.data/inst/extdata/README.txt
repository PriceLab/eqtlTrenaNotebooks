---- bounds

   currently chr5:88391000-89322000 based on enhancers, padded, rounded to nearest 1000 bp


---- enhancer regions:  covers 920kb
   tbl.mef2c.eLocs.RData
   head(tbl.mef2c.eLocs)
           chrom    start      end
     66007  chr5 88822470 88835034
     66197  chr5 89300676 89303600
     66370  chr5 88659269 88669144
     66755  chr5 88891684 88895274
     67465  chr5 88396337 88397974
     68176  chr5 88753320 88753716

   min.enhancer.loc <- with(tbl.mef2c.eLocs, min(c(start, end))) # 88396337
   max.enhancer.loc <- with(tbl.mef2c.eLocs, max(c(start, end))) # 89316026
   1 + with(tbl.mef2c.eLocs, max(c(start, end))) - with(tbl.mef2c.eLocs, min(c(start, end))) #  919690
   enhancer.roiString <- sprintf("chr5:%d-%d", min.enhancer.loc-5000, max.enhancer.loc + 5000)
   tv <- trenaViz(10000:10200)
   setGeneome(tv, "hg38")
   showGenomicRegion(tv, enhancer.roiString)
   provenance:
      data (not in repo) ~/github/dockerizedMicroservices/trena.mef2c/dev/trns/geneHancer/
      code moved to ~/github/eqtlTrenaNotebooks/mef2c/MEF2C.data/inst/dataPrep/getEnhancersAndTheirMotifs.R

   enhancers, compared to variants: extra 244kb upstream, 314kb downstream, roughtly triple the size

   ---> get variants from the vcf file

--- enhancer motifs: mapped to TFs both with MotifDb and TFClass
   print(load("~/github/eqtlTrenaNotebooks/mef2c/trena/data/tbl.eMotifs.RData"))  # [1] "tbl.eMotifs.mdb" "tbl.eMotifs.tfc"
   max(tbl.eMotifs.mdb$motifEnd) - min(tbl.eMotifs.mdb$motifStart)  # [1] 919689
   max(tbl.eMotifs.tfc$motifEnd) - min(tbl.eMotifs.tfc$motifStart)  # [1] 919689
   range(tbl.eMotifs.tfc$motifRelativeScore)  # [1] 0.85 1.00
   range(tbl.eMotifs.mdb$motifRelativeScore)  # [1] 0.85 1.00


--- variants from eqtl study: covers 361kb
   f <- ~/github/dockerizedMicroservices/trena.mef2c/trena/data/mef2c.tsv
   provenance:
      ~/github/projects/external/mayo.AD.collaboration/020516_TableForCorySeth_AD_eQTL_Loci.xlsx
        liftover to  ~/github/projects/external/mayo.AD.collaboration/mef2c-snps-hg38.tsv
        and then to ~/github/projects/external/mayo.AD.collaboration/tbl.snp.hg38.score-ref-alt.RData

       print(load("~/github/projects/external/mayo.AD.collaboration/tbl.snp.hg38.score-ref-alt.RData"))
       dim(tbl.snp) # 155 18
       min.snp.loc <- min(tbl.snp$pos)
       max.snp.loc <- max(tbl.snp$pos)
       snp.roiString <- sprintf("chr5:%d-%d", min.snp.loc-5000, max.snp.loc+5000)


--- vcf data reduction
* associating mef2c gene models, expression data, vcf genotype (20 dec 2017)

  --- sample names in the vcf file:  integers between 129 and 18239

     cd ~/github/projects/priceLab/cory/mef2c-vcf
     start.loc <- 88010000
     end.loc   <- 88201000
     mef2c.gr <- GRanges(5, IRanges(start.loc, end.loc))
     params <- ScanVcfParam(which=mef2c.gr)

     vcfFilename <- "SCH_11923_B01_GRM_WGS_2017-04-27_5.recalibrated_variants.vcf.gz"
     tabixFile <- TabixFile(vcfFilename)
     vcf <- readVcf(tabixFile, "hg19", params)
     vcf.geno <- geno(vcf)
     mtx.geno <- vcf.geno$GT   # 1687  349
     mtx.012 <- matrix(0, nrow=nrow(mtx.geno), ncol=ncol(mtx.geno), dimnames=list(rownames(mtx.geno), colnames(mtx.geno)))
     mtx.012[which(mtx.geno=="0/1")] <- 1
     mtx.012[which(mtx.geno=="1/1")] <- 2
     mtx.geno <- mtx.012

    sort(as.integer(colnames(mtx.geno)))
       [1]   129   134   135   141   142   149   716   719   731   732   736   742   744   746   751
      [16]   766   775   777   783   785   786   787   791   797   812   813   816   828   834   843
      [31]   849   850   851   859   871   886   892   894   896   931   933   937   946   948   952
      ...
     [316] 18201 18202 18203 18204 18206 18209 18210 18211 18212 18213 18214 18216 18217 18218 18219
     [331] 18220 18221 18222 18223 18224 18225 18226 18228 18229 18230 18231 18232 18233 18234 18235
     [346] 18236 18237 18238 18239

    samples.geno <- colnames(mtx.geno)

  --- sample names in a promising expression matrix
     ~/s/work/priceLab/AD/expression.matrix.prep
     prepareMatrices.R starts off like this:
        load("ampADMayo.64253genes.278samples.RData")
        stopifnot(dim(mtx) == c(64253, 278))
        covariates.file <- "AMP-AD_MayoBB_UFL-Mayo-ISB_IlluminaHiSeq2000_TCX_Covariates_Flowcell_1.csv"
     and produces
        print(load("~/s/work/priceLab/AD/expression.matrix.prep/prepped.tcx.matrices.RData"))
           "mtx.tcx"     "mtx.tcx.ctl" "mtx.tcx.ad"
	fivenum(mtx.tcx); fivenum(mtx.tcx.ctl); fivenum(mtx.tcx.ad)
          [1] 0.000000e+00 2.269386e+00 1.384377e+01 4.185981e+01 6.446054e+05
          [1] -9.965784  1.152255  3.758768  5.350882 19.298057
          [1] -9.965784  1.233155  3.814262  5.387233 17.440484
       dim(mtx.tcx); dim(mtx.tcx.ctl); dim(mtx.tcx.ad)
          [1] 18281   278
          [1] 18281    80
          [1] 18281    84
       mtx.tcx.normalized <- asinh(mtx.tcx)
       fivenum(mtx.tcx.normalized)
          [1]  0.000000  1.558003  3.322285  4.427616 14.069541
       head(colnames(mtx.tcx.normalized))
          [1] "S11344_TCX" "S11316_TCX" "S11431_TCX" "S11341_TCX" "S11289_TCX" "S11327_TCX"
       samples.mtx <- sub("^S", "", colnames(mtx.tcx.normalized))
       samples.mtx <- sub("_TCX", "", samples.mtx)
       colnames(mtx.tcx.normalized) <- samples.mtx

       length(samples.mtx)  # [1] 278
       length(samples.geno) # [1] 349
       length(intersect(samples.mtx, samples.geno)) # [1] 263

    covariates.file <- "~/s/work/priceLab/AD/expression.matrix.prep/AMP-AD_MayoBB_UFL-Mayo-ISB_IlluminaHiSeq2000_TCX_Covariates_Flowcell_1.csv"
    tbl.covariates <- read.table(covariates.file, sep=",", as.is=TRUE, header=TRUE)
    dim(tbl.covariates) # [1] 278   8
    as.data.frame(table(tbl.covariates$Diagnosis))
                     Var1 Freq
       1               AD   84
       2          Control   80
       3 Pathologic Aging   30
       4              PSP   84
    tbl.covariates$sample <- gsub("_TCX", "", tbl.covariates$ID)
    length(intersect(samples.geno, tbl.covariates$sample))  # [1] 263
    tbl.pheno <- subset(tbl.covariates, sample %in% samples.geno)  # 263 9
    as.data.frame(table(tbl.pheno$Diagnosis))
                  Var1 Freq
    1               AD   79
    2          Control   73
    3 Pathologic Aging   30
    4              PSP   81
    samples.ad <- subset(tbl.pheno, Diagnosis=="AD")$sample
    samples.ctl <- subset(tbl.pheno, Diagnosis=="Control")$sample
    cor(mtx.tcx.normalized["MEF2C", samples.ctl], mtx.tcx.normalized["NFATC3", samples.ctl])  # [1] -0.8153529
    cor(mtx.tcx.normalized["MEF2C", samples.ad], mtx.tcx.normalized["NFATC3", samples.ad])    # [1] -0.7929382

    --- create hg38 data.frame version of mtx.geno, adding chrom start end columns


   save(tbl.pheno, tbl.geno, samples.ad, samples.ctl, mtx.tcx.normalized, file="~/github/dockerizedMicroservices/trena.mef2c/trena/data/mtx.tcx.pheno.geno.RData")
    load("~/github/dockerizedMicroservices/trena.mef2c/trena/data/mtx.tcx.pheno.geno.RData")


---- eqtl snps

   cd ~/github/eqtlTrenaNotebooks/mef2c/MEF2C.data/inst/extdata/
   cp    ~/github/projects/external/mayo.AD.collaboration/tbl.snp.hg38.score-ref-alt.RData    .
   print(load("~/github/projects/external/mayo.AD.collaboration/tbl.snp.hg38.score-ref-alt.RData"))  # [1] "tbl.snp"
   dim(tbl.snp) # [1] 155  18
     head(tbl.snp)
               rsid chrom      pos      score iupac ref alt A1 CER_Beta    CER_P TX_Beta  TX_P IGAP_A1 IGAP_OR IGAP_Pvalue RegulomeDB Rsquared_rs254776 Dprime_rs254776
     149  rs7721099  chr5 88640561 0.13018179     Y   C   T  C     0.01 0.722000    0.01 0.756       C    0.99     0.74100          6             0.021           0.213
     144  rs7703782  chr5 88642739 0.07417243     W   A   T  A     0.01 0.722000    0.01 0.756       A    1.00     0.84300          7             0.021           0.213
     7   rs10514301  chr5 88643836 0.07262964     Y   C   T  T     0.01 0.722000    0.01 0.756       T    1.00     0.84600          7             0.021           0.213

--- dhs regions: get just the regions
   chrom <- "chr5"; start <- 88391000; end <- 89322000
   library(trena)
   tbl.dhs <- getEncodeDHSRegions("hg38", "wgEncodeRegDnaseClustered", chrom, start, end, quiet=FALSE)
      [1] connecting to genome-mysql.cse.ucsc.edu/hg38/wgEncodeRegDnaseClustered as genome
      [1] query: select * from wgEncodeRegDnaseClustered where chrom = 'chr5' and chromStart >= 88391000 and chromEnd <= 89322000
      [1] 776 DHS regions reported in 931001 bases, start:end unmodified
   save(tbl.dhs, file="/github/eqtlTrenaNotebooks/mef2c/MEF2C.data/inst/extdata/tbl.dhs.RData")

--- get all the motifs - very lengthy!
   library(trena)
   trena <- Trena("hg38")
      # this next step takes a LONG time, maybe an hour, sice it also finds motifs
   x <- getRegulatoryChromosomalRegions(trena, "chr5", 88391000, 89322000, "encodeHumanDHS", "MEF2C", 88904257)
   [1] calling HumanDHSFilter, span: 931001
    > > names(x) # [1] "encodeHumanDHS"
    dim(x$encodeHumanDHS) # [1] 83829     9
    head(x$encodeHumanDHS)
           chrom motifStart motifEnd                           motifName strand     score length distance.from.tss                                                              id
     11528  chr5   88540032 88540047   Hsapiens-jaspar2016-SOX8-MA0868.1      + 0.9166201     16           -364225   MEF2C.dhs.downstream.364225.Hsapiens-jaspar2016-SOX8-MA0868.1
     37170  chr5   88676322 88676342 Hsapiens-jaspar2016-ZNF263-MA0528.1      - 0.9993665     21           -227935 MEF2C.dhs.downstream.227935.Hsapiens-jaspar2016-ZNF263-MA0528.1
     38168  chr5   88676325 88676345 Hsapiens-jaspar2016-ZNF263-MA0528.1      - 0.9993665     21           -227932 MEF2C.dhs.downstream.227932.Hsapiens-jaspar2016-ZNF263-MA0528.1
     39168  chr5   88676328 88676348 Hsapiens-jaspar2016-ZNF263-MA0528.1      - 0.9993665     21           -227929 MEF2C.dhs.downstream.227929.Hsapiens-jaspar2016-ZNF263-MA0528.1
     40164  chr5   88676331 88676351 Hsapiens-jaspar2016-ZNF263-MA0528.1      - 0.9993665     21           -227926 MEF2C.dhs.downstream.227926.Hsapiens-jaspar2016-ZNF263-MA0528.1
     41203  chr5   88676334 88676354 Hsapiens-jaspar2016-ZNF263-MA0528.1      - 0.9993665     21           -227923 MEF2C.dhs.downstream.227923.Hsapiens-jaspar2016-ZNF263-MA0528.1
     tbl.dhsMotifs <- x$encodeHumanDHS
     save(tbl.dhsMotifs, file="tbl.dhsMotifs.RData")

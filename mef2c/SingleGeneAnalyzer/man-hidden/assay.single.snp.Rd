\name{assay.single.snp}
\alias{assay.single.snp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
assay.single.snp(tbl.snp.tf, r, target.gene, tf, tbl.geno, tbl.pheno, mtx)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tbl.snp.tf}{
%%     ~~Describe \code{tbl.snp.tf} here~~
}
  \item{r}{
%%     ~~Describe \code{r} here~~
}
  \item{target.gene}{
%%     ~~Describe \code{target.gene} here~~
}
  \item{tf}{
%%     ~~Describe \code{tf} here~~
}
  \item{tbl.geno}{
%%     ~~Describe \code{tbl.geno} here~~
}
  \item{tbl.pheno}{
%%     ~~Describe \code{tbl.pheno} here~~
}
  \item{mtx}{
%%     ~~Describe \code{mtx} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (tbl.snp.tf, r, target.gene, tf, tbl.geno, tbl.pheno, 
    mtx) 
{
    browser()
    rsid <- tbl.snp.tf$rsid[r]
    printf("\%s \%s:\%s", target.gene, rsid, tf)
    snp.chrom <- tbl.snp.tf$chrom[r]
    snp.loc <- tbl.snp.tf$pos[r]
    tbl.geno.sub <- subset(tbl.geno, chrom == snp.chrom & start == 
        snp.loc)
    mtx.samples <- colnames(mtx)
    wt.samples <- colnames(tbl.geno)[which(subset(tbl.geno, chrom == 
        snp.chrom & start == snp.loc) == 0)]
    wt.samples <- intersect(wt.samples, mtx.samples)
    het.samples <- colnames(tbl.geno)[which(subset(tbl.geno, 
        chrom == snp.chrom & start == snp.loc) == 1)]
    het.samples <- intersect(het.samples, mtx.samples)
    hom.samples <- colnames(tbl.geno)[which(subset(tbl.geno, 
        chrom == snp.chrom & start == snp.loc) == 2)]
    hom.samples <- intersect(hom.samples, mtx.samples)
    mut.samples <- c(het.samples, hom.samples)
    ad.samples <- sub("_TCX", "", subset(tbl.pheno, Diagnosis == 
        "AD")$ID)
    ctl.samples <- sub("_TCX", "", subset(tbl.pheno, Diagnosis == 
        "Control")$ID)
    pathAging.samples <- sub("_TCX", "", subset(tbl.pheno, Diagnosis == 
        "Pathologic Aging")$ID)
    psp.samples <- sub("_TCX", "", subset(tbl.pheno, Diagnosis == 
        "PSP")$ID)
    ad.samples <- intersect(ad.samples, mtx.samples)
    ctl.samples <- intersect(ctl.samples, mtx.samples)
    pathAging.samples <- intersect(pathAging.samples, mtx.samples)
    psp.samples <- intersect(psp.samples, mtx.samples)
    tbl.plot <- data.frame(target = mtx[target.gene, ], tf = mtx[tf, 
        ], geno = "black", pheno = 1, stringsAsFactors = FALSE)
    rownames(tbl.plot) <- colnames(mtx)
    tbl.plot[ctl.samples, "pheno"] <- 15
    tbl.plot[ad.samples, "pheno"] <- 16
    tbl.plot[pathAging.samples, "pheno"] <- 17
    tbl.plot[psp.samples, "pheno"] <- 18
    tbl.plot[wt.samples, "geno"] <- "darkgreen"
    tbl.plot[het.samples, "geno"] <- "pink"
    tbl.plot[hom.samples, "geno"] <- "red"
    with(tbl.plot, plot(tf, target, col = geno, pch = pheno))
    cor(mtx[tf, wt.samples], mtx[target.gene, wt.samples])
    cor(mtx[tf, het.samples], mtx[target.gene, het.samples])
    cor(mtx[tf, hom.samples], mtx[target.gene, hom.samples])
    cor(mtx[tf, mut.samples], mtx[target.gene, mut.samples])
    cor(mtx[tf, ctl.samples], mtx[target.gene, ctl.samples])
    cor(mtx[tf, ad.samples], mtx[target.gene, ad.samples])
    cor(mtx[tf, pathAging.samples], mtx[target.gene, pathAging.samples])
    cor(mtx[tf, psp.samples], mtx[target.gene, psp.samples])
    ctl.wt.samples <- intersect(ctl.samples, wt.samples)
    ctl.mut.samples <- intersect(ctl.samples, mut.samples)
    ctl.het.samples <- intersect(ctl.samples, het.samples)
    ctl.hom.samples <- intersect(ctl.samples, hom.samples)
    length(ctl.wt.samples)
    length(ctl.mut.samples)
    length(ctl.het.samples)
    length(ctl.hom.samples)
    ad.wt.samples <- intersect(ad.samples, wt.samples)
    ad.mut.samples <- intersect(ad.samples, mut.samples)
    ad.het.samples <- intersect(ad.samples, het.samples)
    ad.hom.samples <- intersect(ad.samples, hom.samples)
    length(ad.wt.samples)
    length(ad.mut.samples)
    length(ad.het.samples)
    length(ad.hom.samples)
    pathAging.wt.samples <- intersect(pathAging.samples, wt.samples)
    pathAging.mut.samples <- intersect(pathAging.samples, mut.samples)
    pathAging.het.samples <- intersect(pathAging.samples, het.samples)
    pathAging.hom.samples <- intersect(pathAging.samples, hom.samples)
    length(pathAging.wt.samples)
    length(pathAging.mut.samples)
    length(pathAging.het.samples)
    length(pathAging.hom.samples)
    psp.wt.samples <- intersect(psp.samples, wt.samples)
    psp.mut.samples <- intersect(psp.samples, mut.samples)
    psp.het.samples <- intersect(psp.samples, het.samples)
    psp.hom.samples <- intersect(psp.samples, hom.samples)
    length(psp.wt.samples)
    length(psp.mut.samples)
    length(psp.het.samples)
    length(psp.hom.samples)
    cor(mtx[tf, ctl.wt.samples], mtx[target.gene, ctl.wt.samples])
    cor(mtx[tf, ctl.mut.samples], mtx[target.gene, ctl.mut.samples])
    cor(mtx[tf, ctl.het.samples], mtx[target.gene, ctl.het.samples])
    cor(mtx[tf, ctl.hom.samples], mtx[target.gene, ctl.hom.samples])
    cor(mtx[tf, ad.wt.samples], mtx[target.gene, ad.wt.samples])
    cor(mtx[tf, ad.mut.samples], mtx[target.gene, ad.mut.samples])
    cor(mtx[tf, ad.het.samples], mtx[target.gene, ad.het.samples])
    cor(mtx[tf, ad.hom.samples], mtx[target.gene, ad.hom.samples])
    cor(mtx[tf, pathAging.wt.samples], mtx[target.gene, pathAging.wt.samples])
    cor(mtx[tf, pathAging.mut.samples], mtx[target.gene, pathAging.mut.samples])
    cor(mtx[tf, pathAging.het.samples], mtx[target.gene, pathAging.het.samples])
    cor(mtx[tf, pathAging.hom.samples], mtx[target.gene, pathAging.hom.samples])
    cor(mtx[tf, psp.wt.samples], mtx[target.gene, psp.wt.samples])
    cor(mtx[tf, psp.mut.samples], mtx[target.gene, psp.mut.samples])
    cor(mtx[tf, psp.het.samples], mtx[target.gene, psp.het.samples])
    cor(mtx[tf, psp.hom.samples], mtx[target.gene, psp.hom.samples])
    vec.tf <- mtx[tf, ]
    vec.target <- mtx[target.gene, ]
    printf("genotype for samples in mtx: \%d wt \%d het \%d hom", 
        length(wt.samples), length(het.samples), length(hom.samples))
    vec.tf.wt <- as.numeric(mtx[tf, wt.samples])
    vec.tf.mut <- mtx[tf, c(het.samples, hom.samples)]
    vec.target.wt <- mtx[target.gene, wt.samples]
    vec.target.mut <- mtx[target.gene, c(het.samples, hom.samples)]
    correlation <- cor(vec.tf.wt, vec.target.wt)
    title <- sprintf("\%s-\%s vs \%s wt: \%5.3f", tf, rsid, target.gene, 
        correlation)
    plot(vec.tf.wt, vec.target.wt, ylim = c(0, 8), xlim = c(0, 
        8), main = title)
    correlation <- cor(vec.tf.mut, vec.target.mut)
    title <- sprintf("\%s vs \%s mut: \%5.3f", tf, target.gene, 
        cor(vec.tf.mut, vec.target.mut))
    plot(vec.tf.mut, vec.target.mut, ylim = c(0, 8), xlim = c(0, 
        8), main = title)
    boxplot(vec.tf.wt, vec.tf.mut, vec.target.wt, vec.target.mut)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

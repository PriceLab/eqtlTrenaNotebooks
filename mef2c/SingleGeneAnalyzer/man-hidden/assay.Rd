\name{assay}
\alias{assay}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
assay(tbl.snp.tf, target.gene, tf, tbl.geno, tbl.pheno, mtx)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tbl.snp.tf}{
%%     ~~Describe \code{tbl.snp.tf} here~~
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
function (tbl.snp.tf, target.gene, tf, tbl.geno, tbl.pheno, mtx) 
{
    browser()
    stopifnot(tf \%in\% tbl.snp.tf$geneSymbol)
    tbl.snp.locs <- subset(tbl.snp.tf, geneSymbol == tf)
    all.samples <- colnames(mtx)
    par(mfrow = c(1, 3))
    for (r in 1:nrow(tbl.snp.locs)) {
        rsid <- tbl.snp.locs$rsid[r]
        printf("\%s \%s \%d/\%d", tf, rsid, r, nrow(tbl.snp.locs))
        snp.chrom <- tbl.snp.locs$chrom[r]
        snp.loc <- tbl.snp.locs$pos[r]
        tbl.geno.sub <- subset(tbl.geno, chrom == snp.chrom & 
            start == snp.loc)
        if (nrow(tbl.geno.sub) == 0) 
            next
        wt.samples <- colnames(tbl.geno)[which(subset(tbl.geno, 
            chrom == snp.chrom & start == snp.loc) == 0)]
        wt.samples <- intersect(wt.samples, all.samples)
        het.samples <- colnames(tbl.geno)[which(subset(tbl.geno, 
            chrom == snp.chrom & start == snp.loc) == 1)]
        het.samples <- intersect(het.samples, all.samples)
        hom.samples <- colnames(tbl.geno)[which(subset(tbl.geno, 
            chrom == snp.chrom & start == snp.loc) == 2)]
        hom.samples <- intersect(hom.samples, all.samples)
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
        Sys.sleep(5)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

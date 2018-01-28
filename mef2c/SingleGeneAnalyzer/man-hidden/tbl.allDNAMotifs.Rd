\name{tbl.allDNAMotifs}
\alias{tbl.allDNAMotifs}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data("tbl.allDNAMotifs")}
\format{
  A data frame with 233600 observations on the following 7 variables.
  \describe{
    \item{\code{chrom}}{a factor with levels \code{chr5}}
    \item{\code{start}}{a numeric vector}
    \item{\code{end}}{a numeric vector}
    \item{\code{motifName}}{a character vector}
    \item{\code{strand}}{a character vector}
    \item{\code{motifRelativeScore}}{a numeric vector}
    \item{\code{geneSymbol}}{a character vector}
  }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(tbl.allDNAMotifs)
## maybe str(tbl.allDNAMotifs) ; plot(tbl.allDNAMotifs) ...
}
\keyword{datasets}

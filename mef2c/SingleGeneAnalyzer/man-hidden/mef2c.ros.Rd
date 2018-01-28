\name{mef2c.ros}
\alias{mef2c.ros}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data("mef2c.ros")}
\format{
  A data frame with 20 observations on the following 12 variables.
  \describe{
    \item{\code{gene}}{a character vector}
    \item{\code{beta.lasso}}{a numeric vector}
    \item{\code{beta.ridge}}{a numeric vector}
    \item{\code{rf.score}}{a numeric vector}
    \item{\code{beta.sqrtlasso}}{a numeric vector}
    \item{\code{lasso.p.value}}{a numeric vector}
    \item{\code{pearson.coeff}}{a numeric vector}
    \item{\code{spearman.coeff}}{a numeric vector}
    \item{\code{concordance}}{a numeric vector}
    \item{\code{pcaMax}}{a numeric vector}
    \item{\code{target.gene}}{a character vector}
    \item{\code{Rank}}{a numeric vector}
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
data(mef2c.ros)
## maybe str(mef2c.ros) ; plot(mef2c.ros) ...
}
\keyword{datasets}

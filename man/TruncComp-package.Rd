\name{TruncComp-package}
\alias{TruncComp-package}
\alias{TruncComp}
\docType{package}

\title{\packageTitle{TruncComp}}

\description{\packageDescription{TruncComp}}

\details{
The DESCRIPTION file:
\packageDESCRIPTION{TruncComp}
\packageIndices{TruncComp}
~~ An overview of how to use the package, including the most important functions ~~
}

\author{
\packageAuthor{TruncComp}

Maintainer: \packageMaintainer{TruncComp}
}


\references{
~~ Literature or other references for background information ~~
}


\keyword{ package }

\examples{
library(TruncComp)
d <- simTruncData(n = 100, mu0 = 3, mu1 = 3.5, pi0 = 0.6, pi1 = 0.5)
truncComp(d$Y, d$A, d$Z, method="LRT")
truncComp(d$Y, d$A, d$Z, method="SPLRT")
}

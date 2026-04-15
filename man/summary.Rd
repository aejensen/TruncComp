\name{summary.TruncComp}
\alias{summary.TruncComp}
\title{Summarize a TruncComp fit}
\usage{
\method{summary}{TruncComp}(object, ...)
}
\description{
Prints the estimated treatment contrasts, their confidence intervals, and the
joint likelihood-ratio test result for a fitted \code{TruncComp} object.
}
\arguments{
\item{object}{A \code{TruncComp} object returned by \code{\link{truncComp}}.}
\item{...}{Unused additional arguments.}
}

\value{
The input object, returned invisibly.
}

\examples{
library(TruncComp)
d <- simTruncData(n = 40, mu0 = 3, mu1 = 3.5, pi0 = 0.6, pi1 = 0.5)
fit <- truncComp(Y ~ R, atom = 0, data = d, method = "SPLRT")
summary(fit)
}

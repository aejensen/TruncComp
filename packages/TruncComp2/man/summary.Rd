\name{summary.TruncComp2}
\alias{summary.TruncComp2}
\title{Summarize a TruncComp2 fit}
\usage{
\method{summary}{TruncComp2}(object, ...)
}
\description{
Prints the estimated treatment contrasts, their confidence intervals, and the
joint likelihood-ratio test result for a fitted \code{TruncComp2} object. For
successful unadjusted fits, this includes the stored \code{DeltaMarginalCI}
and the default grid-based \code{DeltaProfileCI}. The projected interval
\code{DeltaProjectedCI} remains available on demand through
\code{\link{confint.TruncComp2}}.
}
\arguments{
\item{object}{A \code{TruncComp2} object returned by \code{\link{truncComp}}.}
\item{...}{Unused additional arguments.}
}

\value{
The input object, returned invisibly.
}

\examples{
library(TruncComp2)
f0 <- function(n) stats::rnorm(n, 3, 1)
f1 <- function(n) stats::rnorm(n, 3.5, 1)
d <- simulateTruncatedData(n = 40, f0 = f0, f1 = f1, pi0 = 0.6, pi1 = 0.5)
fit <- truncComp(Y ~ R, atom = 0, data = d, method = "SPLRT")
summary(fit)
}

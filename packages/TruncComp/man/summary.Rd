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
f0 <- function(n) stats::rnorm(n, 3, 1)
f1 <- function(n) stats::rnorm(n, 3.5, 1)
d <- simulateTruncatedData(n = 40, f0 = f0, f1 = f1, pi0 = 0.6, pi1 = 0.5)
fit <- truncComp(Y ~ R, atom = 0, data = d, method = "SPLRT")
summary(fit)
}

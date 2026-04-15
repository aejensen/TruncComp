\name{print.TruncComp}
\alias{print.TruncComp}

\title{Print a TruncComp fit}

\usage{
\method{print}{TruncComp}(x, ...)
}

\description{
Prints a fitted \code{TruncComp} object using the package summary format.
}

\arguments{
\item{x}{A \code{TruncComp} object returned by \code{\link{truncComp}}.}
\item{...}{Unused additional arguments.}
}

\author{Andreas Kryger Jensen <aeje@sund.ku.dk>}

\examples{
library(TruncComp)
f0 <- function(n) stats::rnorm(n, 3, 1)
f1 <- function(n) stats::rnorm(n, 3.5, 1)
d <- simulateTruncatedData(n = 40, f0 = f0, f1 = f1, pi0 = 0.6, pi1 = 0.5)
fit <- truncComp(Y ~ R, atom = 0, data = d, method = "SPLRT")
print(fit)
}

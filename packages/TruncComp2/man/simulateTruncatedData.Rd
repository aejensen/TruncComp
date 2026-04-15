\name{simulateTruncatedData}
\alias{simulateTruncatedData}

\title{Simulate two-group truncated outcome data}

\usage{
simulateTruncatedData(n, f0, f1, pi0, pi1)
}

\description{
Simulates two groups of data with a binary observation indicator and a
continuous outcome among observed subjects. Unobserved outcomes are encoded as
zero in the returned \code{Y} column.
}

\arguments{
\item{n}{Number of observations per group.}
\item{f0}{Function that generates observed outcomes for group \code{R = 0}.
It may return a finite numeric vector of length \code{n} when called with
\code{n}, or a single finite numeric value when called repeatedly with
\code{1}.}
\item{f1}{Function that generates observed outcomes for group \code{R = 1}.
It may return a finite numeric vector of length \code{n} when called with
\code{n}, or a single finite numeric value when called repeatedly with
\code{1}.}
\item{pi0}{Probability that the outcome is observed in group \code{R = 0}.}
\item{pi1}{Probability that the outcome is observed in group \code{R = 1}.}
}

\value{
A \code{data.frame} with columns:
\item{R}{Binary treatment indicator.}
\item{A}{Binary indicator for whether the continuous outcome is observed.}
\item{Y}{Observed outcome, with unobserved values encoded as zero.}
}

\examples{
f0 <- function(n) stats::rnorm(n, 3, 1)
f1 <- function(n) stats::rnorm(n, 3.5, 1)

d <- simulateTruncatedData(25, f0 = f0, f1 = f1, pi0 = 0.35, pi1 = 0.6)
head(d)
}

\name{simTruncData}
\alias{simTruncData}

\title{Simulate some truncated data}


\description{
hej
}


\usage{simTruncData(n, mu0, mu1, pi0, pi1, sigma = 1, dist = "norm", CV = 1.5)}
\arguments{
  \item{n}{?}
  \item{mu0}{?}
  \item{mu1}{?}
  \item{pi0}{?}
  \item{pi1}{?}
  \item{sigma}{?}
  \item{dist}{?}
  \item{CV}{?}
}


\value{
Returns a \code{data.frame} with column names
\item{Y}{?}
\item{Z}{?}
\item{A}{?}
}

\author{Andreas Kryger Jensen <aeje@sund.ku.dk>}

\examples{
d <- simTruncData(50, 3, 1, 0.5, 0.5)
}

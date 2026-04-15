\name{TruncComp2Example}
\alias{TruncComp2Example}
\docType{data}

\title{Example truncated-outcome data}

\description{
Example data set used throughout the package documentation and manuscript to
illustrate the two-sample comparison procedures.
}

\usage{
data(TruncComp2Example)
}

\format{
A \code{data.frame} with 50 observations and 2 variables:
\describe{
  \item{\code{R}}{Binary treatment indicator coded \code{0} or \code{1}.}
  \item{\code{Y}}{Observed combined outcome, where the atom value \code{0}
  represents an unobserved or undefined continuous outcome.}
}
}

\examples{
data(TruncComp2Example)
head(TruncComp2Example)
}

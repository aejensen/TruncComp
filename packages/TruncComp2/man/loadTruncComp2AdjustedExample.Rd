\name{loadTruncComp2AdjustedExample}
\alias{loadTruncComp2AdjustedExample}

\title{Load the packaged adjusted TruncComp2 example data}

\usage{
loadTruncComp2AdjustedExample()
}

\description{
Loads the packaged example data used to illustrate covariate adjustment in the
parametric likelihood-ratio test.
}

\value{
A \code{data.frame} with 50 observations and 3 variables:
\describe{
  \item{\code{R}}{Binary treatment indicator coded \code{0} or \code{1}.}
  \item{\code{L}}{Categorical baseline covariate with levels \code{low},
  \code{mid}, and \code{high}.}
  \item{\code{Y}}{Observed combined outcome, where the atom value \code{0}
  represents an unobserved or undefined continuous outcome.}
}
}

\details{
The example data are stored as an \code{.rds} file under
\code{inst/extdata/TruncComp2AdjustedExample.rds}. This helper provides a
stable public way to load that file without relying on \code{data()}.
}

\examples{
d <- loadTruncComp2AdjustedExample()
head(d)
}

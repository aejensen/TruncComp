\name{loadTruncComp2Example}
\alias{loadTruncComp2Example}

\title{Load the packaged TruncComp2 example data}

\usage{
loadTruncComp2Example()
}

\description{
Loads the packaged example data used in the \code{TruncComp2} documentation,
tests, and examples.
}

\value{
A \code{data.frame} with 50 observations and 2 variables:
\describe{
  \item{\code{R}}{Binary treatment indicator coded \code{0} or \code{1}.}
  \item{\code{Y}}{Observed combined outcome, where the atom value \code{0}
  represents an unobserved or undefined continuous outcome.}
}
}

\details{
The example data are stored as an \code{.rds} file under
\code{inst/extdata/TruncComp2Example.rds}. This helper provides a stable public
way to load that file without relying on \code{data()}.
}

\examples{
d <- loadTruncComp2Example()
head(d)
}

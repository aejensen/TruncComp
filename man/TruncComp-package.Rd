\name{TruncComp-package}
\alias{TruncComp-package}
\alias{TruncComp}
\docType{package}

\title{\packageTitle{TruncComp}}

\description{\packageDescription{TruncComp}}

\details{
The package provides two-sample comparison procedures for continuous outcomes
with a distinguished atom value representing an unobserved or undefined outcome.

The main user entry points are:
\itemize{
  \item{\code{\link{truncComp}} for fitting the parametric and semi-parametric models}
  \item{\code{\link{summary.truncComp}} and \code{\link{print.truncComp}} for reporting results}
  \item{\code{\link{confint.TruncComp}} for marginal confidence intervals and simultaneous confidence regions}
}

The semi-parametric method uses an internal empirical-likelihood implementation
for the two-sample mean-difference problem and no longer requires the external
\code{EL} package at runtime.
}

\author{
\packageAuthor{TruncComp}

Maintainer: \packageMaintainer{TruncComp}
}


\references{
Andreas Kryger Jensen and Theis Lange.
\emph{A novel high-power test for continuous outcomes truncated by death}.
}


\keyword{ package }

\examples{
library(TruncComp)
d <- simTruncData(n = 100, mu0 = 3, mu1 = 3.5, pi0 = 0.6, pi1 = 0.5)
truncComp(Y ~ R, atom = 0, data = d, method = "LRT")
truncComp(Y ~ R, atom = 0, data = d, method = "SPLRT")
}

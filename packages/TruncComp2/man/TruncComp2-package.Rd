\name{TruncComp2-package}
\alias{TruncComp2-package}
\alias{TruncComp2}
\docType{package}

\title{\packageTitle{TruncComp2}}

\description{\packageDescription{TruncComp2}}

\details{
The package provides two-sample comparison procedures for continuous outcomes
with a distinguished atom value representing an unobserved or undefined outcome.

The main user entry points are:
\itemize{
  \item{\code{\link{truncComp}} for fitting the parametric and semi-parametric models}
  \item{\code{\link{summary.TruncComp2}} and \code{\link{print.TruncComp2}} for reporting results}
  \item{\code{\link{confint.TruncComp2}} for marginal confidence intervals and simultaneous confidence regions}
  \item{\code{\link{simulateTruncatedData}} for generating example two-group data}
  \item{\code{\link{loadTruncComp2Example}} for loading the packaged example data}
  \item{\code{\link{loadTruncComp2AdjustedExample}} for loading the packaged adjusted example data}
}

The semi-parametric method uses an internal empirical-likelihood implementation
for the two-sample mean-difference problem and no longer requires the external
\code{EL} package at runtime.

The helper \code{\link{jointContrastCI}} is intended for fitted
\code{method = "SPLRT"} objects and computes a simultaneous confidence-region
surface for the semi-parametric procedure.

For illustrating covariate adjustment in the parametric model, the package also
ships \code{\link{loadTruncComp2AdjustedExample}}, a fixed example with a
3-level categorical covariate.
}

\author{
\packageAuthor{TruncComp2}

Maintainer: \packageMaintainer{TruncComp2}
}


\references{
Andreas Kryger Jensen and Theis Lange.
\emph{A novel high-power test for continuous outcomes truncated by death}.
}


\keyword{ package }

\examples{
library(TruncComp2)
f0 <- function(n) stats::rnorm(n, 3, 1)
f1 <- function(n) stats::rnorm(n, 3.5, 1)
d <- simulateTruncatedData(n = 100, f0 = f0, f1 = f1, pi0 = 0.6, pi1 = 0.5)
truncComp(Y ~ R, atom = 0, data = d, method = "LRT")
truncComp(Y ~ R, atom = 0, data = d, method = "SPLRT")
}

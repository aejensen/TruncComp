\name{truncComp}
\alias{truncComp}
\title{Two-sample comparison of truncated continuous outcome}
\usage{
truncComp(y, a, z, method, conf.level = 0.95, init = NULL)
}
\description{
Two-sample comparison of truncated continuous outcome by parametric of semi-parametric maximum likelihood.
}
\arguments{
\item{y}{outcome}
\item{a}{truncation indicator}
\item{z}{binary treatment indicator}
\item{method}{\code{LRT} for parametric maximum likelihood and \code{SPLRT} for semi-parametric maximum likelihood}
\item{conf.level}{The level for the confidence intervals (default = 0.95)}
\item{init}{An optinal named list of initial values for the parametric maximum likelihood optimaztion. If uninitialized the function tried to calculate sensible initial values by moment estimation.}
}

\value{
Returns an S3 object of class \code{truncComp}. The object includes the following elements:
\item{muDelta}{hej}
\item{muDeltaCI}{hej}
\item{alphaDelta}{hej}
\item{alphaDeltaCI}{hej}
\item{W}{hej}
\item{p}{hej}
\item{method}{hej}
\item{conf.level}{hej}
\item{success}{hej}
\item{error}{hej}
\item{init}{hej}
}

\author{Andreas Kryger Jensen <aeje@sund.ku.dk>}

\seealso{
\code{\link{print.truncComp}}, \code{\link{summary.truncComp}}
}

\examples{
library(TruncComp)
d <- simTruncData(n = 100, mu0 = 3, mu1 = 3.5, pi0 = 0.6, pi1 = 0.5)
truncComp(d$Y, d$A, d$Z, method="LRT")
truncComp(d$Y, d$A, d$Z, method="SPLRT")
}

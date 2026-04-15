# TruncComp2
Development version of the R package TruncComp2 for two-sample comparison of truncated continuous outcomes.

The package implements:

- a parametric likelihood-ratio test (`method = "LRT"`)
- a semi-parametric likelihood-ratio test (`method = "SPLRT"`)

The current implementation keeps the same public API for both methods, but now computes both paths internally:

- the parametric method uses `glm` and observed-outcome `lm` fits with ML log-likelihood comparison, plus explicit fallbacks for singular boundary cases
- both methods can adjust for baseline covariates through a shared `adjust = ~ ...` specification, returning conditional treatment effects
- the semi-parametric method uses an internal pure-R empirical-likelihood implementation for both the unadjusted two-sample mean difference and the adjusted observed-outcome treatment coefficient

As a result, the package no longer depends on `EL` or `bbmle` at runtime.

The currently supported scope is:

- one binary treatment indicator coded `0/1`
- one atom value representing the unobserved or undefined outcome
- optional additive baseline-covariate adjustment for both `method = "LRT"` and `method = "SPLRT"` through `adjust = ~ ...`
- adjusted `SPLRT` provides fitted tests and marginal confidence intervals, but not simultaneous confidence regions

To install the development version of TruncComp2 run the following commands from within R

```r
library(remotes)
install_github("aejensen/TruncComp", subdir = "packages/TruncComp2")
```

# Documentation

- Statistical model specification: [MODEL.md](MODEL.md)
- Implementation walkthrough: [IMPLEMENTATION.md](IMPLEMENTATION.md)
- Adjusted semi-parametric methodology note: [ADJUSTED_SPLRT.md](ADJUSTED_SPLRT.md)
- Package-local development guide: [DEVELOPMENT.md](DEVELOPMENT.md)
- Packaged example data loader: `loadTruncComp2Example()`
- Packaged adjusted example data loader: `loadTruncComp2AdjustedExample()`

# Main Interface

The primary entry point is:

```r
truncComp(Y ~ R, atom = 0, data = d, method = "LRT")
truncComp(Y ~ R, atom = 0, data = d, method = "LRT", adjust = ~ age + sex)
truncComp(Y ~ R, atom = 0, data = d, method = "SPLRT")
truncComp(Y ~ R, atom = 0, data = d, method = "SPLRT", adjust = ~ age + sex)
```

The fitted object reports:

- `muDelta`: difference in means among the observed
- `alphaDelta`: odds ratio of being observed
- `W`: joint likelihood-ratio test statistic
- `p`: joint p-value

When `adjust` is supplied with either method, `muDelta` and `alphaDelta` are
conditional treatment effects from the adjusted observed-outcome and logistic
submodels. In that adjusted setting, `Delta` is not reported and remains `NA`.

For both `method = "LRT"` and `method = "SPLRT"`, simultaneous
confidence-region surfaces are available for unadjusted fits through
`confint(..., type = "simultaneous")`. Adjusted fits currently support only
marginal intervals.

# Example
```r
library(TruncComp2)

#Define the two distributions for the observed data
f0 <- function(n) stats::rnorm(n, 3, 1)
f1 <- function(n) stats::rnorm(n, 3.5, 1)

#Define probabilities of being observed
pi0 <- 0.35
pi1 <- 0.6

#Simulate data
d <- TruncComp2::simulateTruncatedData(25, f0, f1, pi0, pi1)

#Estimate parameters using the parametric method
fit_lrt <- truncComp(Y ~ R, atom = 0, data = d, method = "LRT")
summary(fit_lrt)
confint(fit_lrt, type = "marginal")

#Load the fixed adjusted example and compare unadjusted vs adjusted LRT
d_adjusted <- loadTruncComp2AdjustedExample()
fit_lrt_unadjusted_example <- truncComp(Y ~ R, atom = 0, data = d_adjusted[, c("Y", "R")], method = "LRT")
fit_lrt_adjusted <- truncComp(Y ~ R, atom = 0, data = d_adjusted, method = "LRT", adjust = ~ L)
summary(fit_lrt_unadjusted_example)
summary(fit_lrt_adjusted)

#The same adjusted example also works with the semi-parametric method
fit_splrt_unadjusted_example <- truncComp(Y ~ R, atom = 0, data = d_adjusted[, c("Y", "R")], method = "SPLRT")
fit_splrt_adjusted <- truncComp(Y ~ R, atom = 0, data = d_adjusted, method = "SPLRT", adjust = ~ L)
summary(fit_splrt_unadjusted_example)
summary(fit_splrt_adjusted)

#Estimate parameters using the semi-parametric method
fit_splrt <- truncComp(Y ~ R, atom = 0, data = d, method = "SPLRT")
summary(fit_splrt)

#Get simultaneous confidence region for an unadjusted fit
confint(fit_lrt, type = "simultaneous", plot = TRUE, resolution = 10)
confint(fit_splrt, type = "simultaneous", plot = TRUE, resolution = 10)
```

# Development

`TruncComp2` is intended to be developed from its own package subtree. For
ordinary work, no edits outside `packages/TruncComp2/` should be needed.

From the package root, the package-local verification entry point is:

```sh
Rscript tools/check-package.R
```

Manual verification is also available:

```sh
R CMD build .
R CMD check --no-manual --no-build-vignettes TruncComp2_*.tar.gz
```

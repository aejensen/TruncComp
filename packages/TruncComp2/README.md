# TruncComp2
Development version of the R package TruncComp2 for two-sample comparison of truncated continuous outcomes.

The package implements:

- a parametric likelihood-ratio test (`method = "LRT"`)
- a semi-parametric likelihood-ratio test (`method = "SPLRT"`)

The current implementation keeps the same public API for both methods, but now computes both paths internally:

- the parametric method uses `glm` and observed-outcome `lm` fits with ML log-likelihood comparison, plus explicit fallbacks for singular boundary cases
- the parametric method can also adjust for baseline covariates through a shared `adjust = ~ ...` specification, returning conditional treatment effects
- the semi-parametric method uses an internal pure-R empirical-likelihood implementation for the observed-outcome mean difference

As a result, the package no longer depends on `EL` or `bbmle` at runtime.

The currently supported scope is:

- one binary treatment indicator coded `0/1`
- one atom value representing the unobserved or undefined outcome
- optional baseline-covariate adjustment for `method = "LRT"` through `adjust = ~ ...`
- no covariate-adjusted `SPLRT` implementation yet

To install the development version of TruncComp2 run the following commands from within R

```r
library(remotes)
install_github("aejensen/TruncComp", subdir = "packages/TruncComp2")
```

# Documentation

- Statistical model specification: [MODEL.md](MODEL.md)
- Implementation walkthrough: [IMPLEMENTATION.md](IMPLEMENTATION.md)
- Package-local development guide: [DEVELOPMENT.md](DEVELOPMENT.md)
- Packaged example data loader: `loadTruncComp2Example()`
- Packaged adjusted example data loader: `loadTruncComp2AdjustedExample()`

# Main Interface

The primary entry point is:

```r
truncComp(Y ~ R, atom = 0, data = d, method = "LRT")
truncComp(Y ~ R, atom = 0, data = d, method = "LRT", adjust = ~ age + sex)
truncComp(Y ~ R, atom = 0, data = d, method = "SPLRT")
```

The fitted object reports:

- `muDelta`: difference in means among the observed
- `alphaDelta`: odds ratio of being observed
- `W`: joint likelihood-ratio test statistic
- `p`: joint p-value

When `adjust` is supplied with `method = "LRT"`, `muDelta` and `alphaDelta`
are conditional treatment effects from the adjusted linear and logistic
submodels. In that adjusted setting, `Delta` is not reported and remains `NA`.

For `method = "SPLRT"`, simultaneous confidence-region surfaces are available through `confint(..., type = "simultaneous")`.

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

#Estimate parameters using the semi-parametric method
fit_splrt <- truncComp(Y ~ R, atom = 0, data = d, method = "SPLRT")
summary(fit_splrt)

#Get simultaneous confidence region
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

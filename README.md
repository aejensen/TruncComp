# TruncComp

`TruncComp` is an R package for comparing two treatment groups when a
patient-centered continuous outcome is undefined after death or another
substantive event. Examples include quality-of-life, functional, or cognitive
scores measured at a fixed follow-up time when some participants have died.

The package represents the observed endpoint by two components

- whether the continuous outcome is defined and
- the value of that outcome when it is defined.

Its primary frequentist procedures jointly test whether treatment groups differ
in either component. They provide one joint $p$-value while retaining separate,
clinically interpretable estimates of the conditional non-atom mean difference
and the odds ratio for being outside the atom. This is an observed-data
comparison. It is neither a survivor-only analysis nor a principal-stratum
analysis.

The package accompanies the manuscript [*Two-component treatment comparisons
for continuous outcomes truncated by a substantive
event*](manuscript/manuscript.pdf). The frequentist procedures constitute its
primary interface. The Bayesian interface is experimental.

## Main features

- `method = "lrt"` combines a Bernoulli likelihood for the non-atom indicator
  with a homoscedastic Gaussian working likelihood for the conditional outcome.
- `method = "splrt"` replaces the Gaussian criterion with empirical likelihood
  mean restrictions and otherwise leaves the conditional non-atom distributions
  unspecified.
- Both methods support additive adjustment for baseline covariates and provide
  component confidence intervals and a joint component confidence region.
- Unadjusted analyses can report three confidence intervals for a mean contrast
  obtained after assigning the atom a numerical code.
- An experimental Bayesian interface estimates the two arm-specific observed
  data laws using flexible finite mixture models.

When baseline covariates define finitely many categories, they can be combined
into a joint stratum factor and included with a saturated representation. The
semiparametric nuisance model is then fully nonparametric across strata. The
probabilities of being outside the atom and the conditional non-atom means may
vary freely by stratum. The two-degree-of-freedom formulation still assumes a
common conditional log odds ratio and a common conditional non-atom mean
difference across strata.

## Installation

Install from GitHub with

```r
install.packages("remotes")
remotes::install_github("aejensen/TruncComp")
```

The package contains Stan models, so installation requires the compilation
toolchain needed by `rstan`. For a local clone, the repository root is also the
package root and the source can be loaded without installing it.

```r
pkgload::load_all(".")
```

## Quick start

The formula response is the coded endpoint. It equals the continuous outcome
when that outcome is defined and equals `atom` after the substantive event. The
treatment indicator must be coded `0/1`.

```r
library(TruncComp)
data("trunc_comp_example", package = "TruncComp")

fit <- trunc_comp(
  Y ~ R,
  atom = 0,
  data = trunc_comp_example,
  method = "splrt"
)

coef(fit)
#>    mu_delta alpha_delta       delta
#>   1.8565084   0.5238095   0.0180069

c(statistic = fit$statistic, p.value = fit$p.value)
#>    statistic      p.value
#> 3.109545e+01 1.768924e-07
```

In this illustrative dataset, the estimated non-atom mean is about 1.86 units
higher under treatment, while the odds of being outside the atom are about 0.52
times those under control. These opposing component differences nearly cancel
on the numerically coded endpoint, whose estimated mean difference is 0.018.
The joint test nevertheless detects the component differences because it does
not collapse them onto that single coded scale.

Component intervals and the simultaneous two-component region are available
through `confint()`.

```r
confint(fit)
confint(fit, parameter = "mu_delta")
confint(fit, parameter = "alpha_delta")
confint(fit, parameter = "joint", plot = TRUE)
```

The component intervals are separate one-parameter likelihood ratio intervals.
The joint region gives a simultaneous statement about the non-atom mean
difference and the log odds ratio. The joint test uses a two-degree-of-freedom
reference distribution.

## Covariate adjustment

Supply an additive adjustment formula through `adjust`. The same covariates are
used in the conditional non-atom mean model and the logistic model for the
probability of being outside the atom.

```r
data("trunc_comp_adjusted_example", package = "TruncComp")

fit_adjusted <- trunc_comp(
  Y ~ R,
  atom = 0,
  data = trunc_comp_adjusted_example,
  method = "splrt",
  adjust = ~ L
)
```

For adjusted fits, the reported component effects are conditional treatment
coefficients. The package does not report a coded-endpoint mean contrast because
those coefficients alone do not determine a standardized marginal contrast.
Component intervals and the joint component region remain available.

For categorical covariates, construct a factor identifying every observed joint
covariate stratum and use `adjust = ~ stratum` to obtain the saturated
representation described above.

## Coded-endpoint inference

For an unadjusted analysis, `delta` is the treatment-minus-control mean
difference after assigning the atom its supplied numerical value. Its point
estimate is stored in the fitted object, and confidence intervals are computed
on demand.

```r
confint(fit, parameter = "delta", method = "welch")
confint(fit, parameter = "delta", method = "profile")
confint(fit, parameter = "delta", method = "projected")
```

- `welch` gives the ordinary Welch interval for the coded endpoint.
- `profile` profiles the full likelihood criterion for the coded contrast.
- `projected` projects a simultaneous region for the four arm-specific
  component parameters.

These intervals answer a coding-dependent question and are distinct from the
primary two-component analysis. They are unavailable for adjusted fits because
a marginal standardization target has not been specified.

The default interface also accepts separate outcome-storage, non-atom-indicator,
and treatment vectors.

```r
fit_default <- trunc_comp(
  trunc_comp_example$Y,
  trunc_comp_example$A,
  trunc_comp_example$R,
  atom = 0,
  method = "splrt"
)
```

Only outcomes whose non-atom indicator equals one contribute to the conditional
outcome model. Values stored for atom observations are ignored when `atom` is
supplied explicitly.

## Experimental Bayesian interface

`trunc_comp_bayes()` fits separate atom probabilities and arm-specific finite
stick-breaking mixtures for the conditional non-atom laws. The available
kernels cover outcomes on the real line, positive outcomes, bounded continuous
outcomes, and bounded scores recorded on a discrete reporting grid. Bayesian
fits also provide posterior component and coded-scale summaries, density or
score-mass plots, sampler diagnostics, and posterior predictive checks.

The current Bayesian interface does not support covariate adjustment. See the
[Bayesian worked example](vignettes/bayesian-example-data.Rmd) for model choices,
prior specification, and computation. Posterior predictive $p$-values supplied
by the package are model-checking summaries rather than frequentist test
$p$-values.

## Documentation

- The [frequentist worked example](vignettes/likelihood-example-data.Rmd) shows
  the likelihood ratio procedures and their confidence summaries.
- [MODEL.md](docs/MODEL.md) specifies the estimands and statistical procedures.
- [ADJUSTED_SPLRT.md](docs/ADJUSTED_SPLRT.md) describes the adjusted
  semiparametric estimating equations.
- [The manuscript directory](manuscript/README.md) documents the article sources
  and submission builds.

## Citation, license, and support

Until publication details are available, please cite the accompanying
manuscript as

> Jensen AK, Lange T. *Two-component treatment comparisons for continuous
> outcomes truncated by a substantive event*. Manuscript.

Software citation information for an installed version is available with
`citation("TruncComp")`.

`TruncComp` is distributed under the GPL-3 license. Questions, bug reports, and
feature requests can be submitted through the [GitHub issue
tracker](https://github.com/aejensen/TruncComp/issues).

## Repository structure

This is a single-package repository. The package occupies the repository root.
The `manuscript/` directory contains the article and its generated submission
variants, `simulation-study/` contains the simulation workflow, and `docs/`
contains statistical documentation.

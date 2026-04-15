# Statistical Model Implemented in TruncComp

## Purpose and Scope

`TruncComp` implements a two-sample comparison for a continuous outcome that includes a special atom value representing an unobserved or undefined outcome. The motivating example in the manuscript is "truncation by death", where deceased patients are assigned the lowest possible score, but the same structure applies whenever one distinguished outcome value has a qualitatively different meaning than the rest of an otherwise continuous scale.

The current package implementation is narrower than the general framework discussed in the manuscript:

- one binary treatment indicator only
- one atom value only
- no additional covariates in the public interface
- two estimation methods: parametric likelihood-ratio (`method = "LRT"`) and semi-parametric likelihood-ratio (`method = "SPLRT"`)

This document describes the model as it is implemented in the package source, with the manuscript used for interpretation and motivation.

For an engineering-focused walkthrough of the current code paths and the internal empirical-likelihood routine, see [IMPLEMENTATION.md](/Users/czv146/Documents/GitHub/TruncComp/packages/TruncComp/IMPLEMENTATION.md).

## Observed Data and Notation

For each subject `i = 1, ..., n`, the package works with:

- `Y_i`: the observed outcome value on the combined scale
- `A_i`: an indicator that the continuous outcome is observed
- `R_i`: the binary treatment indicator, coded `0` or `1`
- `atom`: the distinguished outcome value used when the continuous outcome is unobserved

The package interface reconstructs `A` from the combined outcome:

```text
A_i = 1 if Y_i != atom
A_i = 0 if Y_i == atom
```

The combined outcome can be written as

```math
\tilde Y_i = Y_i \mathbb{1}(A_i = 1) + \mathcal{E}\mathbb{1}(A_i = 0),
```

where `\mathcal{E}` is the atom supplied to `truncComp()` as `atom`. In the current code, the user typically supplies only `Y` and `R`; `A` is then derived internally from whether `Y` equals the atom.

### Input restrictions enforced by `truncComp()`

The exported formula method in `R/truncComp.R` enforces the following:

- the model must be specified with a formula
- the formula must contain exactly one covariate
- the covariate must be binary with unique values `0` and `1`
- `atom` must be a single numeric value
- each treatment group must contain at least one observed outcome
- the method is not intended for data where all outcomes are observed

An additional internal data check requires at least two observed outcomes in each treatment group before estimation proceeds.

## Target Estimands

The package reports two primary treatment-effect components and one derived contrast.

### 1. Difference in means among observed outcomes: `muDelta`

This is the contrast in the continuous component among subjects with `A = 1`:

```math
\mu_\delta = E[Y \mid A = 1, R = 1] - E[Y \mid A = 1, R = 0].
```

In user-facing output this is labeled:

`Difference in means among the observed`

### 2. Odds ratio of being observed: `alphaDelta`

The binary component models the probability that the continuous outcome is observed. The package reports the treatment effect as an odds ratio:

```math
\alpha_\delta = \exp\left(
\operatorname{logit} P(A = 1 \mid R = 1) -
\operatorname{logit} P(A = 1 \mid R = 0)
\right).
```

In user-facing output this is labeled:

`Odds ratio of being observed`

### 3. Derived combined-outcome mean contrast: `Delta`

The package also computes

```math
\Delta = E[\tilde Y \mid R = 1] - E[\tilde Y \mid R = 0],
```

implemented as the empirical difference in sample means of the observed combined outcome values in the two groups.

Important: the main package hypothesis test is **not** a direct test of `\Delta = 0`. A treatment can change the observation probability and the observed-outcome mean in offsetting directions, making `\Delta` close to zero even when the treatment clearly affects the joint distribution.

## Why the Joint Test Has Two Degrees of Freedom

The core modeling idea is that the combined outcome mixes:

- a binary component: whether the continuous outcome is observed
- a continuous component: the outcome value among those with `A = 1`

Because these components can each carry a treatment effect, the null hypothesis implemented by the package is:

```math
H_0:\ \mu_\delta = 0 \text{ and } \alpha_\delta = 1,
```

or equivalently, no treatment effect on either the continuous component or the observation/survival component.

This is stronger and more relevant than testing only

```math
E[\tilde Y \mid R = 1] = E[\tilde Y \mid R = 0].
```

The manuscript shows why a one-dimensional contrast of the combined mean can miss scientifically important effects when the treatment shifts the two components in opposite directions. That is why the joint tests in `TruncComp` are two-degree-of-freedom procedures.

## Parametric Model (`method = "LRT"`)

The parametric likelihood-ratio implementation lives in `R/LRT.R`.

### Model factorization

The code factorizes the likelihood into:

```math
A_i \mid R_i \sim \operatorname{Bernoulli}(\pi_i)
```

with logistic regression

```math
\operatorname{logit}(\pi_i) = \alpha + \alpha_\delta R_i,
```

and

```math
Y_i \mid (A_i = 1, R_i) \sim N(\mu_i, \sigma^2)
```

with

```math
\mu_i = \mu + \mu_\delta R_i.
```

### Null and alternative models

The implemented nested models are:

- Under `H0`: common observation probability `pi`, common observed-outcome mean `mu`, and common `sigma`
- Under `HA`: arm-specific observation probabilities `pi_0`, `pi_1`, arm-specific observed-outcome means `mu_0`, `mu_1`, and common `sigma`

The package still reports the treatment effects as:

- `muDelta = mu_1 - mu_0`
- `alphaDelta`, the odds ratio implied by `pi_1` and `pi_0`

### Observation-level likelihood contribution

For one observation, the implemented density contribution is

```math
L_i = f_Y(Y_i \mid \mu_i, \sigma)^{A_i}\pi_i^{A_i}(1-\pi_i)^{1-A_i},
```

where `f_Y` is the normal density. This means:

- if `A_i = 1`, the observation contributes both a Bernoulli term and a normal-density term
- if `A_i = 0`, the normal term is raised to the power zero and the contribution is only `(1 - \pi_i)`

The current implementation evaluates this likelihood through its Bernoulli and Normal sufficient statistics rather than by fitting a generic optimizer object.

### Test statistic

The package computes the maximized likelihood-ratio statistic directly:

```math
W = 2(\ell_{HA} - \ell_{H0}),
```

which decomposes into:

```math
W = W_A + W_Y,
```

where:

- `W_A` is the grouped Bernoulli likelihood-ratio statistic for the observation indicator
- `W_Y` is the Normal likelihood-ratio statistic for the observed outcomes with common variance

The returned p-value is computed from a chi-squared reference distribution with 2 degrees of freedom.

### Reported quantities

The parametric method returns:

- `muDelta` and its confidence interval
- `alphaDelta` as an odds ratio and its confidence interval
- `Delta` as the empirical mean difference in the combined outcome
- `W` and `p` for the joint likelihood-ratio test

## Semi-Parametric Model (`method = "SPLRT"`)

The semi-parametric implementation lives in `R/SPLRT.R`.

### Continuous component

Among subjects with `A = 1`, the package compares the two treatment groups using an internal empirical-likelihood routine for the mean difference. The current implementation profiles the two-sample empirical likelihood directly inside the package rather than delegating to an external package.

This produces:

- the estimate of `muDelta`
- a confidence interval for `muDelta`
- an empirical likelihood-ratio statistic for the continuous component, denoted here by `W_mu`

The continuous outcome distribution is otherwise left unspecified beyond the mean constraint.

### Binary component

For the observation indicator, the package fits:

- `m0: A ~ 1`
- `m1: A ~ R`

using logistic regression via `glm(..., family = binomial())`.

From this step it obtains:

- `alphaDelta = exp(coef(m1)["R"])`
- a confidence interval for `alphaDelta`
- a likelihood-ratio statistic for the observation model, denoted here by `W_alpha`

### Joint statistic

The package forms the joint semi-parametric test statistic additively:

```math
W = W_\mu + W_\alpha.
```

It then computes the p-value from a chi-squared reference distribution with 2 degrees of freedom:

```math
p = 1 - F_{\chi^2_2}(W).
```

This additive construction is also used by the simultaneous confidence-region code in `R/CI.R`, where the joint surface is evaluated over a grid of candidate values for the two treatment-effect parameters.

### Implementation notes

The internal empirical-likelihood engine is implemented in `R/empiricalLikelihood.R` and is intentionally narrower than the former external dependency. It only supports the two-sample mean-difference problem needed by `TruncComp`.

At a high level, the implementation:

- validates that both observed-outcome samples are numeric, finite, and contain at least two observations
- treats `mean(x) - mean(y)` as the exact empirical-likelihood estimator
- profiles the nuisance mean parameter `theta` for any candidate difference `delta`
- solves the one-sample empirical-likelihood lambda equations with `stats::uniroot()`
- falls back to `stats::optimize()` for the nuisance parameter if the `theta` gradient is numerically awkward
- constructs the confidence interval by solving `W_mu(delta) = qchisq(conf.level, 1)` over the feasible parameter range
- uses a cached logistic baseline when evaluating simultaneous confidence-region surfaces so the unconstrained logistic fit is not recomputed at every grid point

This design keeps the `SPLRT` output compatible with the previous package behavior while removing the runtime dependency on `EL`.

## Estimation and Confidence Intervals

For both methods, the returned model object contains:

- `muDelta`
- `muDeltaCI`
- `alphaDelta`
- `alphaDeltaCI`
- `Delta`
- `DeltaCI`
- `W`
- `p`

### What is currently implemented

- `muDelta` and `alphaDelta` are estimated for both `LRT` and `SPLRT`
- marginal confidence intervals for `muDelta` and `alphaDelta` are implemented for both methods
- simultaneous confidence regions are implemented only for the semi-parametric method
- `jointContrastCI()` is exported for compatibility but is only defined for successful `SPLRT` fits
- `confint()` uses the fitted model's `conf.level` by default and treats a different requested marginal confidence level as requiring a refit

### Current implementation detail: `DeltaCI`

`Delta` is computed in both methods as the empirical difference in sample means of the combined outcome, but `DeltaCI` is not currently implemented. In both `R/LRT.R` and `R/SPLRT.R`, it is set to:

```r
c(NA, NA)
```

The summary and confidence-interval methods also omit `Delta` from the displayed treatment-contrast table.

## Code-Level Mapping

| Concept | Symbol in this document | Object/output name | Where implemented |
| --- | --- | --- | --- |
| Difference in means among observed outcomes | `\mu_\delta` | `muDelta` | `R/LRT.R`, `R/SPLRT.R` |
| Confidence interval for observed-outcome mean difference | `CI(\mu_\delta)` | `muDeltaCI` | `R/LRT.R`, `R/SPLRT.R` |
| Odds ratio of being observed | `\alpha_\delta` | `alphaDelta` | `R/LRT.R`, `R/SPLRT.R` |
| Confidence interval for observation odds ratio | `CI(\alpha_\delta)` | `alphaDeltaCI` | `R/LRT.R`, `R/SPLRT.R` |
| Combined-outcome mean difference | `\Delta` | `Delta` | `R/LRT.R`, `R/SPLRT.R` |
| Combined-outcome CI placeholder | `CI(\Delta)` | `DeltaCI` | `R/LRT.R`, `R/SPLRT.R` |
| Joint test statistic | `W` | `W` | `R/LRT.R`, `R/SPLRT.R` |
| Joint p-value | `p` | `p` | `R/LRT.R`, `R/SPLRT.R` |
| Parametric likelihood-ratio test | parametric joint LRT | `method = "LRT"` | `R/LRT.R`, called from `R/truncComp.R` |
| Semi-parametric likelihood-ratio test | additive empirical/logistic LRT | `method = "SPLRT"` | `R/SPLRT.R`, called from `R/truncComp.R` |
| Simultaneous confidence-region surface | joint contrast surface | `surface` in `jointContrastCI()` output | `R/CI.R` |

## Assumptions and Limitations

The following restrictions are part of the current implementation, regardless of broader possibilities discussed in the manuscript:

- treatment must be binary and coded `0/1`
- the formula interface accepts exactly one covariate
- only one atom value can be specified
- there must be at least one observed outcome in each treatment group
- the internal data check effectively requires at least two observed outcomes in each treatment group
- the method is not applicable when all outcomes are observed
- normality of the observed continuous outcomes is assumed only for `method = "LRT"`
- `method = "SPLRT"` leaves the observed-outcome distribution unspecified beyond mean restrictions
- no covariate adjustment is currently available in the public package interface
- `DeltaCI` is not implemented
- simultaneous confidence regions are only implemented for `SPLRT`

## Worked Example

The package example data contains the combined outcome `Y` and the treatment indicator `R`. The observation indicator is inferred from the atom value.

```r
library(TruncComp)
data("TruncCompExample")

atom <- 0
model <- truncComp(Y ~ R, atom = atom, data = TruncCompExample, method = "SPLRT")
summary(model)
```

Interpretation:

- `muDelta` estimates how much the mean outcome differs between treatment groups among subjects with observed non-atom outcomes
- `alphaDelta` estimates the treatment effect on the odds of being observed rather than taking the atom value
- `W` and `p` provide the joint test of no treatment effect on either component

This differs from a standard two-sample t-test or Wilcoxon test on `Y` alone. Those procedures work directly on the combined outcome scale and can miss important effects when the treatment changes the observation probability and the observed-outcome distribution in different directions.

## Source-of-Truth Notes

This document was cross-checked against:

- `R/truncComp.R` for the exported interface and input restrictions
- `R/LRT.R` for the parametric likelihood specification
- `R/SPLRT.R` for the semi-parametric construction
- `R/CI.R` for the simultaneous confidence-region machinery
- `R/classFunctions.R` for user-facing labels
- `manuscript/manuscript.Rmd` and `manuscript/supplementary.Rmd` for motivation and theoretical background

When the manuscript describes a more general framework than the package currently exposes, this document follows the package implementation.

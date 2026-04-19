# Statistical Model Implemented in TruncComp2

## Purpose and Scope

`TruncComp2` implements a two-sample comparison for a continuous outcome that includes a special atom value representing an unobserved or undefined outcome. The motivating example in the manuscript is "truncation by death", where deceased patients are assigned the lowest possible score, but the same structure applies whenever one distinguished outcome value has a qualitatively different meaning than the rest of an otherwise continuous scale.

The current package implementation is narrower than the general framework discussed in the manuscript:

- one binary treatment indicator only
- one atom value only
- optional additive baseline-covariate adjustment for both `method = "LRT"` and `method = "SPLRT"` through `adjust = ~ ...`
- adjusted `method = "SPLRT"` currently implements fitted tests and component confidence intervals, but not joint confidence regions or `Delta` intervals
- two estimation methods: parametric likelihood-ratio (`method = "LRT"`) and semi-parametric likelihood-ratio (`method = "SPLRT"`)

The package also includes an experimental Bayesian pathway through
`trunc_comp_bayes()`. That path currently supports only the no-covariate
setting and uses a two-part model with:

- arm-specific atom probabilities
- arm-specific truncated stick-breaking mixtures for `Y | A = 1`

The continuous component can be specified with:

- `continuous_support = "real_line"` for Gaussian kernels on the real line
- `continuous_support = "positive_real"` for Gamma kernels on the positive real
  line, requiring all non-atom outcomes to be strictly positive
- `continuous_support = "bounded_continuous"` for an experimental
  Beta-mixture survivor model on `(score_min, score_max)`
- `continuous_support = "bounded_score"` for an experimental discretized/heaped
  Beta-mixture survivor-score model on the reported grid in `[score_min,
  score_max]`

The bounded-score model keeps posterior predictive survivor outcomes on the
same reported support as the observed scores. The bounded-continuous model is
for genuinely continuous or finely measured bounded outcomes; exact survivor
values at `score_min` or `score_max` are rejected in the MVP rather than handled
with boundary inflation.

Its primary posterior outputs are:

- `delta_atom`: difference in atom probability
- `mu_delta`: difference in mean among non-atom outcomes
- `alpha_delta`: odds ratio of being observed
- `delta`: combined-outcome contrast

Unlike the likelihood-ratio paths, the Bayesian implementation reports
posterior summaries and diagnostics rather than p-values or LR statistics.

This document describes the model as it is implemented in the package source, with the manuscript used for interpretation and motivation.

For an engineering-focused walkthrough of the current code paths and the internal empirical-likelihood routine, see [IMPLEMENTATION.md](IMPLEMENTATION.md).

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
- the main formula must contain exactly one binary treatment indicator
- the treatment indicator must be binary with unique values `0` and `1`
- `atom` must be a single numeric value
- each treatment group must contain at least one observed outcome
- the method is not intended for data where all outcomes are observed
- when `adjust` is supplied, it must be a one-sided formula and may not include the outcome or treatment variable

An additional internal data check requires at least two observed outcomes in each treatment group before estimation proceeds.

## Target Estimands

The package reports two primary treatment-effect components and one derived contrast.

### 1. Difference in means among observed outcomes: `muDelta`

Without covariate adjustment, this is the contrast in the continuous component among subjects with `A = 1`:

```math
\mu_\delta = E[Y \mid A = 1, R = 1] - E[Y \mid A = 1, R = 0].
```

For adjusted fits, `muDelta` is instead the conditional treatment coefficient from the adjusted observed-outcome model described below.

In user-facing output this is labeled:

`Difference in means among the observed`

### 2. Odds ratio of being observed: `alphaDelta`

The binary component models the probability that the continuous outcome is observed. Without covariate adjustment, the package reports the treatment effect as an odds ratio:

```math
\alpha_\delta = \exp\left(
\operatorname{logit} P(A = 1 \mid R = 1) -
\operatorname{logit} P(A = 1 \mid R = 0)
\right).
```

For adjusted fits, `alphaDelta = exp(\beta_\delta)` is the conditional odds ratio from the adjusted logistic model.

In user-facing output this is labeled:

`Odds ratio of being observed`

### 3. Derived combined-outcome mean contrast: `Delta`

For unadjusted fits, the package also computes

```math
\Delta = E[\tilde Y \mid R = 1] - E[\tilde Y \mid R = 0],
```

using the fitted atom value `atom = \mathcal{E}` stored on the model object. In
the unadjusted case this model-implied contrast agrees exactly with the
empirical difference in sample means of the combined outcome values in the two
groups.

Important: the main package hypothesis test is **not** a direct test of `\Delta = 0`. A treatment can change the observation probability and the observed-outcome mean in offsetting directions, making `\Delta` close to zero even when the treatment clearly affects the joint distribution.

For adjusted fits, `Delta` is returned as `NA` because the implemented treatment effects are conditional regression coefficients rather than standardized marginal combined-outcome contrasts.

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

The manuscript shows why a one-dimensional contrast of the combined mean can miss scientifically important effects when the treatment shifts the two components in opposite directions. That is why the joint tests in `TruncComp2` are two-degree-of-freedom procedures.

## Parametric Model (`method = "LRT"`)

The parametric likelihood-ratio implementation lives in `R/LRT.R`.

### Model factorization

Without covariate adjustment, the code factorizes the likelihood into:

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

When `adjust = ~ L` is supplied to `truncComp(..., method = "LRT")`, the package extends both submodels with the same additive covariate specification:

```math
\operatorname{logit} P(A_i = 1 \mid R_i, L_i) = \beta_0 + \beta_\delta R_i + s_\beta(L_i)
```

```math
E[Y_i \mid A_i = 1, R_i, L_i] = \mu_0 + \mu_\delta R_i + s_\mu(L_i).
```

The adjusted joint null is still two-dimensional:

```math
H_0:\ \beta_\delta = 0 \text{ and } \mu_\delta = 0.
```

### Null and alternative models

The implemented nested models are:

- Under `H0`: common observation probability `pi`, common observed-outcome mean `mu`, and common `sigma`
- Under `HA`: arm-specific observation probabilities `pi_0`, `pi_1`, arm-specific observed-outcome means `mu_0`, `mu_1`, and common `sigma`

The package still reports the treatment effects as:

- `muDelta = mu_1 - mu_0`
- `alphaDelta`, the odds ratio implied by `pi_1` and `pi_0`

With covariate adjustment, the fitted null and alternative models become:

- logistic null: `A ~ L`
- logistic alternative: `A ~ R + L`
- linear null on observed outcomes: `Y ~ L`
- linear alternative on observed outcomes: `Y ~ R + L`

and the reported treatment effects are the conditional treatment coefficients:

- `muDelta = \mu_\delta`
- `alphaDelta = exp(\beta_\delta)`

### Observation-level likelihood contribution

For one observation, the implemented density contribution is

```math
L_i = f_Y(Y_i \mid \mu_i, \sigma)^{A_i}\pi_i^{A_i}(1-\pi_i)^{1-A_i},
```

where `f_Y` is the normal density. This means:

- if `A_i = 1`, the observation contributes both a Bernoulli term and a normal-density term
- if `A_i = 0`, the normal term is raised to the power zero and the contribution is only `(1 - \pi_i)`

The current implementation evaluates this likelihood by fitting the Bernoulli and observed-outcome components separately:

- `glm(A ~ 1, family = binomial())` and `glm(A ~ R, family = binomial())` for the observation model
- `lm(Y ~ 1)` and `lm(Y ~ R)` on the `A = 1` subset for the observed outcomes

For adjusted fits, the same pattern becomes:

- `glm(A ~ L, family = binomial())` and `glm(A ~ R + L, family = binomial())`
- `lm(Y ~ L)` and `lm(Y ~ R + L)` on the `A = 1` subset

### Test statistic

In regular cases the package computes the maximized likelihood-ratio statistic from model log-likelihoods:

```math
W = 2(\ell_{HA} - \ell_{H0}),
```

which decomposes into:

```math
W = W_A + W_Y,
```

where:

- `W_A` is obtained from the difference in `glm` log-likelihoods for the observation indicator
- `W_Y` is obtained from the difference in ML `lm` log-likelihoods for the observed outcomes with common variance

When the fitted models hit boundary cases such as zero cells in the Bernoulli table or zero residual variance in the observed outcomes, the implementation falls back to explicit grouped-count or residual-sum-of-squares logic so the returned statistic remains stable and interpretable.

The returned p-value is computed from a chi-squared reference distribution with 2 degrees of freedom.

Adjusted parametric fits deliberately require regular model estimation. If the adjusted treatment effect is aliased, has non-finite variance, or the fitted likelihood contribution is non-finite, the package returns a failed `TruncComp2` object rather than attempting separation-specific or penalized fallback logic.

### Reported quantities

The parametric method returns:

- `muDelta` and its confidence interval
- `alphaDelta` as an odds ratio and its confidence interval
- `Delta` as the empirical mean difference in the combined outcome
- `W` and `p` for the joint likelihood-ratio test

For adjusted parametric fits, the reported quantities are:

- `muDelta` and `muDeltaCI` from the treatment coefficient in the adjusted observed-outcome `lm`
- `alphaDelta` and `alphaDeltaCI` from the treatment coefficient in the adjusted logistic `glm`
- `Delta = NA`, and `confint()` rejects `parameter = "Delta"` and `parameter = "joint"`
- `W` and `p` from the summed logistic and linear likelihood-ratio statistics

## Semi-Parametric Model (`method = "SPLRT"`)

The semi-parametric implementation lives in `R/SPLRT.R`.

### Continuous component

Among subjects with `A = 1`, the package compares the two treatment groups using an internal empirical-likelihood routine for the mean difference. The current implementation profiles the two-sample empirical likelihood directly inside the package rather than delegating to an external package.

This produces:

- the estimate of `muDelta`
- a confidence interval for `muDelta`
- an empirical likelihood-ratio statistic for the continuous component, denoted here by `W_mu`

The continuous outcome distribution is otherwise left unspecified beyond the mean constraint.

When `adjust = ~ L` is supplied, the semi-parametric continuous component switches to an empirical-likelihood regression-coefficient test on the observed outcomes. Let

```math
X_i = (1, R_i, h(L_i)^T)
```

and define the estimating equations

```math
g_i(\beta) = X_i \left(Y_i - X_i^T \beta\right).
```

The adjusted semi-parametric implementation:

- uses the OLS treatment coefficient as `muDelta`
- profiles over the nuisance coefficients when testing a candidate treatment coefficient `\delta`
- solves the empirical-likelihood dual problem for the constrained estimating equations
- inverts the profiled statistic to obtain `muDeltaCI`

This is a genuine extension beyond the original manuscript and is documented further in [ADJUSTED_SPLRT.md](ADJUSTED_SPLRT.md).

### Binary component

For the observation indicator, the package fits:

- `m0: A ~ 1`
- `m1: A ~ R`

using logistic regression via `glm(..., family = binomial())`.

From this step it obtains:

- `alphaDelta = exp(coef(m1)["R"])`
- a confidence interval for `alphaDelta`
- a likelihood-ratio statistic for the observation model, denoted here by `W_alpha`

With covariate adjustment, the binary component becomes:

- logistic null: `A ~ L`
- logistic alternative: `A ~ R + L`

and the reported treatment effect remains `alphaDelta = exp(\beta_\delta)`.

### Joint statistic

The package forms the joint semi-parametric test statistic additively:

```math
W = W_\mu + W_\alpha.
```

It then computes the p-value from a chi-squared reference distribution with 2 degrees of freedom:

```math
p = 1 - F_{\chi^2_2}(W).
```

This additive construction is also used by the joint confidence-region code in
`R/CI.R`, where the joint surface is evaluated over a grid of candidate values
for the two treatment-effect parameters.

For adjusted `SPLRT`, the additive statistic is still used, but joint
confidence regions are not implemented in this version.

### Implementation notes

The internal empirical-likelihood engine is implemented in `R/empiricalLikelihood.R` and now supports two related problems needed by `TruncComp2`:

- the original two-sample mean-difference empirical likelihood
- an empirical-likelihood regression-coefficient profile used by adjusted `SPLRT`

At a high level, the implementation:

- validates that both observed-outcome samples are numeric, finite, and contain at least two observations
- treats `mean(x) - mean(y)` as the exact empirical-likelihood estimator
- profiles the nuisance mean parameter `theta` for any candidate difference `delta`
- solves the one-sample empirical-likelihood lambda equations with `stats::uniroot()`
- falls back to `stats::optimize()` for the nuisance parameter if the `theta` gradient is numerically awkward
- constructs the confidence interval by solving `W_mu(delta) = qchisq(conf.level, 1)` over the feasible parameter range
- uses a cached logistic baseline when evaluating simultaneous confidence-region surfaces so the unconstrained logistic fit is not recomputed at every grid point

This design keeps the unadjusted `SPLRT` output compatible with the previous package behavior while extending the semi-parametric method to additive covariate adjustment without a runtime dependency on `EL`.

## Estimation and Confidence Intervals

For both methods, the returned model object contains:

- `muDelta`
- `muDeltaCI`
- `alphaDelta`
- `alphaDeltaCI`
- `Delta`
- `W`
- `p`

The fitted object stores the component intervals `muDeltaCI` and
`alphaDeltaCI`, but it does not store any `Delta` confidence interval. Instead,
all `Delta` intervals are computed on demand through `confint()`.

### What is currently implemented

- `muDelta` and `alphaDelta` are estimated for both `LRT` and `SPLRT`
- component intervals for `muDelta` and `alphaDelta` are implemented for both methods and returned by `confint(fit)` or `confint(fit, parameter = ...)`
- joint confidence regions are implemented for successful unadjusted `LRT` and `SPLRT` fits through `confint(fit, parameter = "joint")`
- `jointContrastCI()` is exported for compatibility and is defined for successful unadjusted fits from either method
- for successful unadjusted fits, the package also reports `Delta` together with three on-demand interval constructions selected by `method = "welch" | "projected" | "profile"`
- `confint()` uses the fitted model's `conf.level` by default and treats a different requested component confidence level as requiring a refit

### Component intervals

The package's first layer of inference is componentwise.

For `muDelta`, the component confidence interval is a one-dimensional interval
for the observed-outcome treatment effect:

```math
CI_\mu = \{ \delta : W_\mu(\delta) \le q_{\chi^2_1}(\texttt{conf.level}) \}
```

for the semi-parametric path, and the corresponding fitted-model interval for
the parametric path.

In practice:

- unadjusted `LRT` uses the treatment coefficient and fitted variance from the
  observed-outcome `lm`
- unadjusted `SPLRT` uses empirical-likelihood inversion for the two-sample mean
  difference
- adjusted `LRT` uses the treatment coefficient from the adjusted observed-outcome
  `lm`
- adjusted `SPLRT` uses empirical-likelihood inversion for the adjusted
  regression coefficient

For `alphaDelta`, the package works on the log-odds-ratio scale
`\psi = \log(\alpha_\delta)` and reports the interval on the odds-ratio scale:

```math
CI_\alpha = \exp\{ CI_\psi \}.
```

In practice:

- `LRT` uses a fitted logistic-regression interval for the treatment coefficient
  when the information matrix is regular
- unadjusted `LRT` returns `NA` for `alphaDeltaCI` in boundary cases where the
  odds ratio is still estimable as `0` or `Inf` but the regular interval is not
- unadjusted `SPLRT` obtains the logistic interval from the fitted `glm`
- adjusted `SPLRT` uses the same adjusted logistic coefficient logic as the
  adjusted parametric path

Interpretation:

- `muDeltaCI` answers "which values of the observed-outcome treatment effect are
  supported if I focus only on the continuous component?"
- `alphaDeltaCI` answers "which values of the observation-probability treatment
  effect are supported if I focus only on the binary component?"

These are not joint intervals. They do not control the joint uncertainty in the
pair `(muDelta, alphaDelta)`.

These intervals are returned by `confint(fit)`,
`confint(fit, parameter = "muDelta")`, `confint(fit, parameter = "alphaDelta")`,
or `confint(fit, parameter = c("muDelta", "alphaDelta"))`.

### Joint confidence region

For successful unadjusted fits, the package also computes a joint region for
the two primary treatment-effect components:

```math
C_{\mathrm{joint}}
= \{ (\delta_\mu, \delta_\psi) :
W(\delta_\mu, \delta_\psi) \le q_{\chi^2_2}(\texttt{conf.level}) \},
```

where `\delta_\psi = \log(\alpha_\delta)`.

This region is returned by `confint(..., parameter = "joint")` and
`jointContrastCI()`. It is evaluated on a grid because the user-facing goal is a
two-dimensional contour or surface, not a one-dimensional bound.

Interpretation:

- this is the correct object when the inferential target is the pair
  `(muDelta, logORdelta)`
- it is the right starting point for graphical joint inference
- it is also the object from which the projected `Delta` interval is derived

### The atom-dependent combined-outcome contrast `Delta`

For unadjusted fits, `Delta` is defined from the fitted atom `e = atom` as

```math
\Delta(e)
= [p_1 \mu_1 + (1 - p_1)e] - [p_0 \mu_0 + (1 - p_0)e].
```

This makes the dependence on the atom explicit: changing `atom` changes the
meaning of `Delta` and therefore changes every interval for `Delta`.

In the unadjusted case, the fitted point estimate agrees exactly with the raw
combined-outcome mean difference:

```math
\hat\Delta
= \bar Y_{R=1} - \bar Y_{R=0}.
```

For adjusted fits, `Delta` is not implemented because the fitted treatment
effects are conditional coefficients rather than a standardized marginal
contrast.

### Three different intervals for `Delta`

For successful unadjusted fits, `TruncComp2` supports three distinct confidence
interval constructions for `Delta`. They share the same point estimate but have
different mathematical definitions and different interpretations.

#### 1. `confint(fit, parameter = "Delta", method = "welch")`

This is the ordinary two-sample interval for the raw combined-outcome mean
difference, using the observed combined outcome values directly:

```math
\hat\Delta = \bar Y_1 - \bar Y_0,
```

with Welch standard error

```math
\operatorname{SE}(\hat\Delta)
= \sqrt{\frac{s_1^2}{n_1} + \frac{s_0^2}{n_0}},
```

Welch-Satterthwaite degrees of freedom, and the usual

```math
\hat\Delta \pm t_{\mathrm{df},\,1-(1-\texttt{conf.level})/2}\operatorname{SE}(\hat\Delta).
```

Interpretation:

- this is the interval for the raw combined-outcome mean difference treated as
  an ordinary two-sample summary
- it is fast, transparent, and largely agnostic to the joint model structure
- it is often the best descriptive interval when the user wants a simple summary
  of the combined scale

#### 2. `confint(fit, parameter = "Delta", method = "projected")`

This interval is defined by projecting the joint two-parameter region onto the
`Delta` scale:

```math
CI_{\Delta,\mathrm{proj}}
= \left[
\inf_{(\delta_\mu,\delta_\psi)\in C_{\mathrm{joint}}}\Delta(\delta_\mu,\delta_\psi),
\sup_{(\delta_\mu,\delta_\psi)\in C_{\mathrm{joint}}}\Delta(\delta_\mu,\delta_\psi)
\right].
```

Interpretation:

- this is the range of `Delta` values compatible with the full joint region for
  `(muDelta, logORdelta)`
- it inherits the two-parameter joint guarantee and is therefore usually more
  conservative than a one-dimensional profile interval
- it is the right interval when the user wants a `Delta` statement that remains
  tied to the whole joint two-component treatment effect

Numerical implementation in the current package:

- by default, `TruncComp2` evaluates the joint surface on an adaptive grid,
  keeps all points satisfying the `\chi^2_2` cutoff, maps those accepted
  points to the implied `Delta`, and returns the range of the accepted
  `Delta` values
- this default grid-based version is returned by
  `confint(..., parameter = "Delta", method = "projected")`
- a slower direct constrained-optimization alternative remains available
  through `confint(..., parameter = "Delta", method = "projected",
  algorithm = "optimize")`

#### 3. `confint(fit, parameter = "Delta", method = "profile")`

This is the one-dimensional profile interval for `Delta` itself. Define the
profiled statistic

```math
W_\Delta(d)
= \inf_{\eta:\,\Delta(\eta)=d} W(\eta),
```

where `\eta` denotes the nuisance and treatment-effect parameters needed by the
chosen unadjusted model. Then

```math
CI_{\Delta,\mathrm{profile}}
= \{ d : W_\Delta(d) \le q_{\chi^2_1}(\texttt{conf.level}) \}.
```

Interpretation:

- this is the model-based one-dimensional interval when `Delta` itself is the
  inferential target
- unlike the projected interval, it is not a projection of the 2 d.f. joint
  region
- unlike the Welch interval, it is explicitly tied to the fitted truncated-outcome
  model

Numerical implementation in the current package:

- by default, the package approximates this interval from the evaluated joint
  surface by using the `\chi^2_1` cutoff instead of the joint `\chi^2_2`
  cutoff
- a slower direct optimization-based alternative remains available through
  `confint(..., parameter = "Delta", method = "profile",
  algorithm = "optimize")`; that version profiles directly in `\Delta` by
  minimizing the model statistic subject to the atom-dependent constraint
  `\Delta = d`

### Which interval should be used?

- Use `muDeltaCI` if the scientific focus is the observed-outcome component.
- Use `alphaDeltaCI` if the scientific focus is the observation/survival
  component.
- Use the joint region when you want joint inference for both treatment
  effects together.
- Use `method = "welch"` for a simple, descriptive interval on the combined
  outcome scale.
- Use `method = "profile"` when `Delta` is the main inferential target and you
  want the model-based one-dimensional interval. The default public version is
  grid-based; request `algorithm = "optimize"` if you specifically want the
  direct optimization-based profile interval.
- Use `method = "projected"` when you want the range of `Delta` values
  compatible with the full joint statement. The default public version is
  grid-based; use `algorithm = "optimize"` when you explicitly want the slower
  constrained-optimization version.

For adjusted fits, `Delta` is stored as `NA`, and `confint()` rejects
`parameter = "Delta"` and `parameter = "joint"` because no standardized
marginal combined-outcome contrast is implemented for the conditional model.

## Code-Level Mapping

| Concept | Symbol in this document | Object/output name | Where implemented |
| --- | --- | --- | --- |
| Difference in means among observed outcomes | `\mu_\delta` | `muDelta` | `R/LRT.R`, `R/SPLRT.R` |
| Confidence interval for observed-outcome mean difference | `CI(\mu_\delta)` | `muDeltaCI` | `R/LRT.R`, `R/SPLRT.R` |
| Odds ratio of being observed | `\alpha_\delta` | `alphaDelta` | `R/LRT.R`, `R/SPLRT.R` |
| Confidence interval for observation odds ratio | `CI(\alpha_\delta)` | `alphaDeltaCI` | `R/LRT.R`, `R/SPLRT.R` |
| Combined-outcome mean difference | `\Delta` | `Delta` | `R/LRT.R`, `R/SPLRT.R` |
| Welch interval for combined-outcome mean difference | Welch `CI(\Delta)` | `confint(..., parameter = "Delta", method = "welch")` | `R/delta.R`, `R/CI.R` |
| Projected interval from the joint surface | projected `CI(\Delta)` | `confint(..., parameter = "Delta", method = "projected")` | `R/delta.R`, `R/CI.R` |
| Profile interval for combined-outcome mean difference | profile `CI(\Delta)` | `confint(..., parameter = "Delta", method = "profile")` | `R/delta.R`, `R/CI.R` |
| Joint test statistic | `W` | `W` | `R/LRT.R`, `R/SPLRT.R` |
| Joint p-value | `p` | `p` | `R/LRT.R`, `R/SPLRT.R` |
| Parametric likelihood-ratio test | parametric joint LRT | `method = "LRT"` | `R/LRT.R`, called from `R/truncComp.R` |
| Semi-parametric likelihood-ratio test | additive empirical/logistic LRT | `method = "SPLRT"` | `R/SPLRT.R`, called from `R/truncComp.R` |
| Simultaneous confidence-region surface | joint contrast surface | `surface` in `jointContrastCI()` output | `R/CI.R` |

## Assumptions and Limitations

The following restrictions are part of the current implementation, regardless of broader possibilities discussed in the manuscript:

- treatment must be binary and coded `0/1`
- the main formula interface accepts exactly one treatment variable
- only one atom value can be specified
- there must be at least one observed outcome in each treatment group
- the internal data check effectively requires at least two observed outcomes in each treatment group
- the method is not applicable when all outcomes are observed
- normality of the observed continuous outcomes is assumed only for `method = "LRT"`
- `method = "SPLRT"` leaves the observed-outcome distribution unspecified beyond mean restrictions
- the same additive covariate specification is used in both adjusted `LRT` submodels and in the adjusted `SPLRT` logistic and observed-outcome components
- adjusted treatment effects are conditional coefficients with no `R * L` interactions
- `Delta` inference is implemented only for unadjusted fits
- joint confidence regions are only implemented for unadjusted fits

## Worked Example

The package example data contains the combined outcome `Y` and the treatment indicator `R`. The observation indicator is inferred from the atom value.

```r
library(TruncComp2)
d <- loadTruncComp2Example()

atom <- 0
model <- truncComp(Y ~ R, atom = atom, data = d, method = "SPLRT")
summary(model)
```

Interpretation:

- `muDelta` estimates how much the mean outcome differs between treatment groups among subjects with observed non-atom outcomes
- `alphaDelta` estimates the treatment effect on the odds of being observed rather than taking the atom value
- `W` and `p` provide the joint test of no treatment effect on either component

If you supply `adjust = ~ L1 + L2`, the interpretation changes slightly:

- `muDelta` is the conditional treatment coefficient from the adjusted observed-outcome model
- `alphaDelta` is the conditional odds ratio `exp(coef_R)` from the adjusted logistic model
- `W` and `p` test the joint null that both adjusted treatment coefficients are zero

The packaged adjusted example data make this concrete:

```r
d_adjusted <- loadTruncComp2AdjustedExample()

fit_unadjusted <- truncComp(Y ~ R, atom = 0, data = d_adjusted[, c("Y", "R")], method = "LRT")
fit_adjusted <- truncComp(Y ~ R, atom = 0, data = d_adjusted, method = "LRT", adjust = ~ L)
fit_splrt_unadjusted <- truncComp(Y ~ R, atom = 0, data = d_adjusted[, c("Y", "R")], method = "SPLRT")
fit_splrt_adjusted <- truncComp(Y ~ R, atom = 0, data = d_adjusted, method = "SPLRT", adjust = ~ L)
```

In that fixed example, `L` is prognostic for both being observed and for the observed outcome value, and the treatment groups are imbalanced in `L`. As a result, both the parametric and semi-parametric fits attenuate after adjustment. In the semi-parametric case, the unadjusted fit is below the 5% level while the adjusted fit is above it, but the adjusted effects remain in the same direction and are only moderately attenuated rather than disappearing completely.

This differs from a standard two-sample t-test or Wilcoxon test on `Y` alone. Those procedures work directly on the combined outcome scale and can miss important effects when the treatment changes the observation probability and the observed-outcome distribution in different directions.

## Source-of-Truth Notes

This document was cross-checked against:

- `R/truncComp.R` for the exported interface and input restrictions
- `R/LRT.R` for the parametric likelihood specification
- `R/SPLRT.R` for the semi-parametric construction
- `R/CI.R` for the simultaneous confidence-region machinery
- `R/classFunctions.R` for user-facing labels
- the original project manuscript for optional historical motivation and broader theoretical background when available

When older project material describes a more general framework than the package currently exposes, this document follows the package implementation.

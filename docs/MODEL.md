# Statistical model

## Data structure

Let `R_i` be a binary treatment indicator, with `R_i = 0` for control and
`R_i = 1` for treatment. Let `A_i` indicate whether subject `i` is outside the
substantive atom. Thus `A_i = 1` when the continuous outcome is defined and
`A_i = 0` when the substantive event has occurred.

The continuous outcome `Y_i` is defined only when `A_i = 1`. Let `e_0` denote
the numerical atom used to code the substantive event. The coded endpoint used
for inference is

```math
Z_i =
\begin{cases}
Y_i, & A_i = 1, \\
e_0, & A_i = 0.
\end{cases}
```

This piecewise definition does not evaluate `Y_i` when `A_i = 0`. In the
formula interface, the variable on the left-hand side is `Z_i`, although the
examples call that variable `Y` to match the R interface. The package
reconstructs `A_i` as `1(Z_i != e_0)`.

The frequentist default interface supplies an outcome storage vector and `A_i`
separately. Only entries with `A_i = 1` are used as continuous outcomes. When
the atom is supplied explicitly, finite entries in the storage vector with
`A_i = 0` are arbitrary placeholders and are ignored. The package uses `e_0`
to construct `Z_i` for coded-endpoint inference. If the atom is omitted, it can
be inferred only when all placeholder entries share one finite value.

The package requires

- a finite numeric atom
- one treatment variable coded `0/1`
- a default-interface outside-atom indicator coded `0/1`
- at least one atom observation
- at least two non-atom outcomes in each treatment arm
- no missing values in the variables used by the fit
- an additive adjustment formula without treatment interactions when
  covariates are supplied.

## Estimands

For an unadjusted analysis, define

```math
p_r = \operatorname{P}(A_i = 1 \mid R_i = r),
\qquad
\mu_r = \operatorname{E}(Y_i \mid A_i = 1, R_i = r).
```

The two component estimands are

```math
\mu_\delta = \mu_1 - \mu_0
```

and

```math
\alpha_\delta
= \frac{p_1/(1-p_1)}{p_0/(1-p_0)}.
```

The corresponding fitted-object fields are `mu_delta` and `alpha_delta`.
Values of `alpha_delta` above one indicate larger odds of being outside the
atom under treatment.

The coded-endpoint contrast is

```math
\Delta
= \operatorname{E}(Z_i \mid R_i = 1) - \operatorname{E}(Z_i \mid R_i = 0)
= p_1\mu_1 + (1-p_1)e_0
  - p_0\mu_0 - (1-p_0)e_0.
```

Its fitted-object field is `delta`. The joint component test is not a direct
test that the coded-endpoint contrast is zero. The two components can change in
opposite directions and offset on the coded scale even when treatment changes
both parts of the outcome law.

## Covariate-adjusted estimands

Let `L_i` collect the supplied baseline covariates. Both frequentist methods use
the additive models

```math
\operatorname{logit}\{\operatorname{P}(A_i = 1 \mid R_i,L_i)\}
= \beta_0 + \beta_\delta R_i + \beta_L^\mathsf{T}L_i
```

and

```math
\operatorname{E}(Y_i \mid A_i = 1,R_i,L_i)
= \gamma_0 + \gamma_\delta R_i + \gamma_L^\mathsf{T}L_i.
```

The adjusted component estimands are `mu_delta = gamma_delta` and
`alpha_delta = exp(beta_delta)`. They are conditional treatment contrasts that
are constant across covariate values under the additive models. All adjustment
coefficients are nuisance parameters.

These conditional coefficients do not determine a marginal coded-endpoint
contrast without a standardization distribution and an additional marginal
estimand. Adjusted fits therefore return `delta = NA`, and coded-endpoint
intervals are unavailable. Joint inference for `(gamma_delta, beta_delta)` is
available and profiles all adjustment coefficients.

## Parametric likelihood ratio procedure

For `method = "lrt"`, the unadjusted model is

```math
A_i \mid R_i=r \sim \operatorname{Bernoulli}(p_r)
```

and

```math
Y_i \mid A_i=1,R_i=r \sim
\mathcal{N}(\mu_r,\sigma^2),
```

where the non-atom variance is common across treatment arms. The contribution
of subject `i` is defined without referring to an undefined outcome:

```math
L_i =
\begin{cases}
1-p_{R_i}, & A_i=0, \\
p_{R_i}\,\phi(Y_i;\mu_{R_i},\sigma^2), & A_i=1.
\end{cases}
```

The likelihood separates into a Bernoulli component and a Gaussian component.
The null model imposes both `p_1 = p_0` and `mu_1 = mu_0`. The alternative
allows both arm differences. The statistic is

```math
W = W_A + W_Y
  = 2\{\ell_A(\widehat\eta_A)-\ell_A(\widehat\eta_{A,0})\}
  + 2\{\ell_Y(\widehat\eta_Y)-\ell_Y(\widehat\eta_{Y,0})\}.
```

Under the joint null and the usual likelihood ratio regularity conditions,
`W` has an asymptotic chi-squared distribution with two degrees of freedom.
The package reports `statistic = W` and

```math
\texttt{p.value}=\operatorname{P}(\chi^2_2\ge W).
```

With covariates, the Bernoulli component is a logistic regression and the
Gaussian component is a linear model fitted to the non-atom outcomes. The null
sets the treatment coefficient to zero in both models while retaining all
adjustment coefficients. The Gaussian comparison uses the maximized log
likelihood rather than the restricted log likelihood. Rank-deficient adjusted
designs and separated logistic fits fail rather than producing nonregular
component estimates.

The component intervals invert the one-parameter profile likelihood ratio
statistics. The interval for `alpha_delta` is profiled on the log odds ratio
scale and transformed to the odds ratio scale. Explicit grouped-count and
residual-sum-of-squares calculations handle unadjusted boundary cases. A
component interval can be `NA` when a finite regular interval does not exist.

## Semiparametric likelihood ratio procedure

For `method = "splrt"`, the Bernoulli component is the same logistic likelihood
used by the parametric procedure. The non-atom component replaces the Gaussian
likelihood with empirical likelihood.

In an unadjusted analysis, the empirical likelihood profiles the two arm means
subject to a fixed value of `mu_delta = mu_1 - mu_0`. Its likelihood ratio
statistic is zero at the difference in sample non-atom means. Inverting that
statistic at the one-degree-of-freedom cutoff gives `mu_delta_ci`.

In an adjusted analysis, let `X_i` be the non-atom design row containing the
treatment indicator and the columns generated by the adjustment formula,
including an intercept when specified. The empirical likelihood uses the
estimating equations

```math
\sum_{i:A_i=1} w_i X_i\{Y_i-X_i^\mathsf{T}\gamma\}=0,
\qquad
\sum_{i:A_i=1}w_i=1,
\qquad
w_i>0.
```

For a candidate treatment coefficient, the nuisance coefficients and empirical
weights are profiled out. The unrestricted estimate equals the ordinary least
squares treatment coefficient, while inference uses the empirical likelihood
profile statistic.

The joint statistic is

```math
W = W_A + W_\mu.
```

Under the joint null and the corresponding empirical likelihood and logistic
regularity conditions, `W` has an asymptotic chi-squared distribution with two
degrees of freedom. The fitted object uses the same field names and null values
as the parametric fit.

## Frequentist result object

A successful `trunc_comp_fit` contains

| Field | Meaning |
|---|---|
| `mu_delta` | Treatment-minus-control non-atom mean contrast |
| `mu_delta_ci` | Component interval for `mu_delta` |
| `alpha_delta` | Treatment-to-control odds ratio for being outside the atom |
| `alpha_delta_ci` | Component interval for `alpha_delta` |
| `delta` | Treatment-minus-control coded-endpoint contrast, or `NA` for an adjusted fit |
| `statistic` | Joint two-component likelihood ratio statistic |
| `p.value` | Chi-squared two-degree-of-freedom p-value |
| `method` | `"lrt"` or `"splrt"` |
| `conf.level` | Confidence level used for stored component intervals |
| `success` | Whether estimation succeeded |
| `data` | Standardized analysis data with columns `Y`, `A`, and `R` plus adjustment variables |
| `adjust` | Printable adjustment specification or `NULL` |
| `adjust_formula` | Formula used to reconstruct the exact adjusted design |
| `atom` | Atom used to code the substantive event |

The stored component intervals use the confidence level supplied at fitting.
Requesting a different level requires refitting. Coded-endpoint and joint-region
calculations use the confidence level supplied to `confint()`.

## Joint component region

For either frequentist method, define `psi_delta = log(alpha_delta)`. The joint
component region is

```math
\mathcal{C}_{\mathrm{joint}}
= \{(m,s): W(m,s) \le q_{\chi^2_2}(\texttt{conf.level})\}.
```

At each pair `(m,s)`, the remaining parameters are profiled out. For adjusted
fits, this includes every adjustment coefficient in both components. The
functions

```r
confint(fit, parameter = "joint")
joint_contrast_surface(fit)
```

evaluate this criterion on a grid in `(mu_delta, log_or_delta)`. This region is
the simultaneous inferential object for the two treatment components. It is not
the full parameter region used for coded-endpoint inference.

## Coded-endpoint intervals

Coded-endpoint intervals are available for successful unadjusted fits. They all
use the same point estimate but answer different questions.

### Welch interval

`confint(fit, parameter = "delta", method = "welch")` applies the ordinary
Welch construction to the coded observations `Z_i`. It uses the arm-specific
sample variances and Welch-Satterthwaite degrees of freedom. It is descriptive
and does not depend on the fitted two-part criterion.

### Profile interval

Let

```math
\zeta=(\mu_0,\mu_1,\eta_0,\eta_1),
\qquad \eta_r=\operatorname{logit}(p_r),
```

and let `W_full(zeta)` be the full unadjusted criterion relative to its
unrestricted minimum. The coded-endpoint profile statistic is

```math
W_\Delta(d)=\inf_{\zeta:\Delta(\zeta)=d}W_{\mathrm{full}}(\zeta).
```

The profile interval is

```math
\{d:W_\Delta(d)\le q_{\chi^2_1}(\texttt{conf.level})\}.
```

Because a continuous image of the connected likelihood sublevel is an
interval, this accepted set equals the range of `Delta(zeta)` over the full
sublevel with the one-degree-of-freedom cutoff.

For `lrt`, fixing the two logits leaves a Gaussian mean ellipsoid whose coded
minimum and maximum are analytic. A deterministic radial search over the
compact Bernoulli logit sublevel determines the endpoints. For `splrt`, the
factored Bernoulli and conditional empirical likelihood criterion reduces
exactly to ordinary two-sample empirical likelihood for the coded observations,
so the package inverts that statistic at the one-degree-of-freedom cutoff.

### Projected interval

The projected interval is the range of `Delta(zeta)` over

```math
\{\zeta:W_{\mathrm{full}}(\zeta)
\le q_{\chi^2_4}(\texttt{conf.level})\}.
```

It is a projection of a simultaneous region for both arm means and both arm
logits. The calculation uses the same parametric reduction or coded empirical
likelihood reduction as the profile interval, with the four-degree-of-freedom
cutoff. It is generally wider than the scalar profile interval.

## Bayesian model

`trunc_comp_bayes()` implements a no-covariate two-part model. Each treatment
arm has an atom probability `rho_r`, probability outside the atom
`p_r = 1-rho_r`, and an arm-specific finite stick-breaking mixture for the
conditional non-atom law.

The Bayesian default interface treats its outcome vector as the encoded
endpoint. It therefore requires `y[a == 0]` to equal the atom and `y[a == 1]`
to differ from the atom. This differs from the arbitrary finite placeholders
accepted by the frequentist default interface.

The non-atom kernel family is selected by `continuous_support` and
`bounded_kernel`.

| Support | Kernel |
|---|---|
| `real_line` | Gaussian mixture |
| `positive_real` | Gamma mixture |
| `bounded_continuous` with `beta` | Beta mixture transformed to the supplied open interval |
| `bounded_continuous` with `logit_normal` | Logit-normal mixture transformed to the supplied open interval |
| `bounded_score` with `beta` | Discretized Beta mixture on the complete score grid |
| `bounded_score` with `logit_normal` | Discretized logit-normal mixture on the complete score grid |

The generic logit-normal models assign independent priors to component
locations on the logit scale. They do not impose an ordering constraint.

For bounded scores, each kernel probability is integrated over the reporting
cell associated with a score. `heaping_grids` specifies candidate reporting
resolutions, and the model estimates either shared or arm-specific weights over
those resolutions. The grid determined by `score_min`, `score_max`, and
`score_step` must contain every non-atom score in the data. Every heaping grid
must be a multiple of `score_step` and must divide the supplied score range.
Every observed non-atom score must be compatible with at least one candidate
heaping grid.

The Bayesian estimands are

```math
\texttt{delta_atom}=\rho_1-\rho_0,
\qquad
\texttt{mu_delta}=\mu_1-\mu_0,
```

```math
\texttt{alpha_delta}
=\frac{p_1/(1-p_1)}{p_0/(1-p_0)},
\qquad
\texttt{delta}
=p_1\mu_1+(1-p_1)e_0-p_0\mu_0-(1-p_0)e_0.
```

The Bayesian fit reports posterior medians, equal-tail credible intervals,
posterior probabilities relative to the relevant null value, arm-specific
summaries, sampler diagnostics, and posterior draws. It does not report the
frequentist likelihood ratio statistic or hypothesis-test p-value. Its
posterior predictive p-values are discrepancy-based model checks.

## Scope

- Treatment is binary and coded `0/1`.
- The formula interface contains one treatment variable.
- The frequentist default interface accepts a numeric outcome storage vector,
  an outside-atom indicator, and a treatment indicator. Finite values stored at
  atom rows are ignored when the atom is supplied explicitly.
- The Bayesian default interface accepts the same three vectors but requires
  the outcome vector and outside-atom indicator to encode the supplied atom
  consistently.
- Both frequentist methods support the same additive adjustment specification.
- Treatment-by-covariate interactions are outside the implemented scope.
- Joint component regions support adjusted and unadjusted frequentist fits.
- Coded-endpoint contrasts and intervals support unadjusted fits only.
- The Bayesian interface does not support covariate adjustment.

## Source mapping

| Functionality | Source |
|---|---|
| Public frequentist interface and validation | `R/trunc_comp.R`, `R/utility.R` |
| Parametric criterion and component intervals | `R/LRT.R` |
| Semiparametric criterion | `R/SPLRT.R`, `R/empiricalLikelihood.R` |
| Joint component regions | `R/CI.R`, `R/delta.R` |
| Coded-endpoint intervals | `R/CI.R`, `R/delta.R` |
| Bayesian interface and summaries | `R/trunc_comp_bayes.R`, `R/bayesFit.R` |
| Bounded support and reporting grids | `R/bayesBounded.R` |
| Posterior predictive checks | `R/bayesPPC.R` |

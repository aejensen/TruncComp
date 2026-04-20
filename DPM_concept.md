# Dirichlet Process Mixture Model Concept for `TruncComp`

## Purpose

This note describes a Bayesian concept for a future Dirichlet process mixture model (DPMM) pathway in `TruncComp`. It is not a description of the current `LRT` or `SPLRT` implementation. Instead, it translates the recommendations in `DP-implementation.md` into a coherent implementation concept that can guide a later fully specified design document.

The intended audience is technical collaborators who need a compact statement of the recommended model architecture, the expected computational structure, and the remaining questions that must be resolved before coding begins. Notation is aligned with the current package documentation in `packages/TruncComp2/MODEL.md` and `packages/TruncComp2/IMPLEMENTATION.md`.

## Current Baseline

The current `TruncComp2` package already uses a two-part view of the outcome:

- `Y` is the combined observed outcome.
- `atom` is the distinguished value representing the undefined or unobserved outcome.
- `A = 1(Y != atom)` indicates whether the continuous outcome is observed.
- `R in {0,1}` is the treatment indicator.

Conceptually, the existing package separates:

- the observation or survival component through `A`
- the continuous outcome component through `Y | A = 1`

The proposed Bayesian DPMM pathway should preserve this decomposition. The new model is therefore best understood as a Bayesian replacement for the current two-part logic, not as a completely unrelated model.

## Recommended v1 Concept

### Core modeling decision

The recommended first implementation is an explicit two-part Bayesian model fitted separately by treatment arm. For each arm `r in {0,1}`,

```math
Y \mid R = r \sim p_r \delta_a + (1-p_r) F_r^c,
```

where:

- `a` is the atom value
- `p_r = P(Y = a | R = r)` is the arm-specific atom probability
- `F_r^c` is the arm-specific continuous distribution on the non-atom support

This is the recommended baseline because it preserves the scientific distinction between the atom event and the continuous outcome, matches the existing `TruncComp` framing, and is simpler to interpret and implement than a single augmented DP over the full response space.

### Continuous component

The non-atom distribution in each treatment arm is modeled with a Dirichlet process mixture:

```math
G_r^c \sim DP(\alpha_r H_r),
\qquad
F_r^c(B) = \int_{\Theta_c} K_c(B \mid \theta)\, G_r^c(d\theta),
\qquad B \subseteq C = S \setminus \{a\}.
```

The kernel family `K_c` must match the support of the observed continuous outcome:

- Gaussian kernels for outcomes on all of `\mathbb{R}`
- positive-support kernels or transformed Gaussian kernels for outcomes on `(0, \infty)`
- interval-respecting transformed kernels when the outcome lives on a bounded range

The same kernel family should be used in both arms so that posterior contrasts remain comparable, while allowing `G_0^c` and `G_1^c` to differ freely.

### Atom component

The atom probability should be modeled directly rather than induced indirectly through an augmented mixture construction. A simple conceptual prior is

```math
p_r \sim Beta(a_p, b_p),
```

though the exact prior family and hyperparameters remain open design choices.

### Why this formulation is preferred

The v1 concept explicitly chooses the two-part model over an augmented one-DP formulation. The augmented formulation can be mathematically related to the explicit formulation under matched priors, but it is not the preferred first implementation because:

- the atom and continuous parts are less transparent as separate scientific objects
- prior control over the atom probability is less direct
- the implementation is less modular
- standard continuous-mixture machinery is harder to reuse cleanly

The recommended architecture is therefore:

- one explicit atom model per arm
- one continuous DP mixture per arm
- posterior treatment comparison through functionals of the fitted arm-specific response distributions

## Conceptual Interface

The implementation should be designed around a small set of observed inputs, latent model components, and posterior outputs.

### Observed inputs

- `Y`: combined outcome
- `R`: binary treatment indicator
- `atom`: distinguished atom value
- optional future `L`: baseline covariates for later extensions

### Latent model components

- `p_0`, `p_1`: atom probabilities by arm
- `G_0^c`, `G_1^c`: arm-specific DP mixing measures
- `F_0^c`, `F_1^c`: induced arm-specific non-atom outcome distributions
- mixture component parameters and latent allocation variables for non-atom observations

### Posterior estimands and outputs

The model should support posterior inference for three primary effects:

```math
\Delta_A = p_1 - p_0,
```

```math
\Delta_c = \mu_1^c - \mu_0^c,
\qquad
\mu_r^c = \int_C y\,F_r^c(dy),
```

```math
\Delta_{ext} =
\left[a p_1 + (1-p_1)\mu_1^c\right]
-
\left[a p_0 + (1-p_0)\mu_0^c\right].
```

These correspond to:

- the atom effect
- the difference in mean outcome among non-atom observations
- the combined extended-outcome effect

The combined extended-outcome effect should be the headline single-number summary when one summary is needed, with the atom and non-atom components reported alongside it.

Posterior reporting should emphasize:

- posterior means or medians
- credible intervals
- posterior probabilities such as `P(Delta_ext > 0 | data)`

The default reporting target should not be a p-value analogue, and Bayes factors should not be the default headline output.

## Computational Design Concept

The first computational implementation should be conservative and modular.

### Preferred computation strategy

The preferred initial fitting approach is a truncated stick-breaking approximation with latent allocation variables for the continuous non-atom mixture:

```math
v_h \sim Beta(1,\alpha_r),
\qquad
w_h = v_h \prod_{g < h}(1-v_g),
\qquad
\theta_h \sim H_r.
```

This leads to a finite computational approximation while keeping the conceptual model close to a genuine DP mixture.

### Practical design principles

- Fit the atom model and continuous non-atom mixture as distinct but linked parts.
- Run the continuous mixture only on observations with `Y != atom`.
- Keep arm-specific mixtures independent in v1 rather than sharing atoms across arms.
- Use support-matched kernels that assign zero mass to the atom.
- Prefer posterior summaries and uncertainty intervals over hypothesis-test style outputs.

This strategy is recommended because it is easier to code, test, and maintain in an R package than exact infinite samplers, collapsed samplers, or fully dependent DP constructions.

## Roadmap Beyond v1

The memo supports a staged extension path, but these extensions should remain outside the baseline implementation concept.

### Categorical covariates

For categorical `L`, two directions are supported conceptually:

- stratified two-part models when the number of levels is small and cells are well populated
- hierarchical pooling, potentially with a hierarchical Dirichlet process, when levels are sparse and stratum-specific distributions remain scientifically important

If the scientific goal is a common adjusted conditional treatment effect rather than heterogeneous stratum-specific effects, then the model should impose a common-shift or no-interaction structure rather than relying on pure stratification.

### Continuous or mixed covariates

For continuous or mixed covariates, the preferred first extension is additive regression in both model parts:

- logistic or probit regression for the atom probability
- regression for the non-atom mean together with a DP mixture for residuals

Conceptually,

```math
Y_i = m_r(L_i) + \varepsilon_i,
\qquad
\varepsilon_i \mid G_r \sim \int K(\cdot \mid \theta)\, G_r(d\theta).
```

This keeps covariate effects interpretable while using the DP mixture to absorb skewness, heavy tails, and multimodality in the residual distribution.

### Later, not v1

The following are explicitly later-stage possibilities rather than part of the baseline concept:

- dependent DP models where the full mixture distribution varies with `L`
- cross-arm sharing of mixture atoms
- more ambitious hierarchical borrowing structures
- exact slice-based or franchise-style samplers

## Questions Needed for a Full Implementation Spec

The concept above resolves the overall model direction, but a fully specified implementation document still needs answers to the following questions.

### Statistical specification

- What outcome supports must v1 officially support: all real values, positive outcomes only, bounded scales, or more than one of these?
- Which kernel family should be the default for each supported outcome type?
- Should the atom model use Beta priors by default, or should v1 instead use a logistic parameterization with priors on transformed coefficients?
- How should the DP concentration parameters be handled: fixed defaults, estimated with hyperpriors, or user-configurable only?
- What base measures `H_r` should be used for mixture component locations and scales?
- Should v1 support only the no-covariate model, or should it already expose one adjusted model variant?
- Is the required headline estimand `Delta_ext`, or should all three estimands be treated as co-primary outputs?
- Should treatment comparison be summarized on difference scales only, or are odds-ratio and other transformed summaries also required?

### Computation and sampling

- What truncation strategy and default truncation level should be used for the stick-breaking approximation?
- Which sampler should be implemented first in practice: Gibbs, Metropolis-within-Gibbs, or another hybrid scheme?
- What posterior summary quantities must be computed during fitting versus derived afterward from saved draws?
- What convergence diagnostics are required for a fit to be considered valid?
- How should failed fits, poor mixing, or empty effective support in one arm be detected and reported?
- What runtime and memory targets are acceptable for the first package implementation?

### Package interface and output design

- What user-facing function should expose the Bayesian pathway, and should it extend `truncComp()` or live behind a separate entry point?
- How should the user specify outcome support and kernel choice, if these are not inferred automatically?
- What posterior draws, summaries, and diagnostics should be stored in the returned object?
- Which printed summaries should be default, and which should require explicit requests?
- How should the output distinguish the atom effect, the non-atom effect, and the extended-outcome effect?
- Should the package return posterior probabilities such as `P(Delta_ext > 0 | data)` by default, or only when requested?

### Validation and testing

- What simulated data scenarios are required to validate the atom model, the non-atom mixture model, and the combined estimands?
- What benchmark comparisons against the current parametric or semi-parametric methods are required?
- What posterior predictive checks should be part of routine validation?
- What regression tests are needed to ensure reproducible results under fixed seeds?
- What acceptance criteria should determine whether the Bayesian pathway is ready to ship as experimental or stable functionality?

## Summary

The memo points to a clear first implementation concept: keep the existing two-part problem structure, model the atom probability explicitly, fit separate arm-specific DP mixtures only to the non-atom outcomes, and report posterior contrasts with the combined extended-outcome effect as the main summary. Richer covariate handling and hierarchical borrowing are important extensions, but they should follow only after the baseline no-covariate design is fully specified and validated.

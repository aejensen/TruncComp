# Phase 4 - Notation Inventory and Uniformity Audit

```text
Construct a complete notation register for manuscript/manuscrip.tex before
implementing notation changes.

For every symbol, record:
1. Symbol and first location.
2. Whether it is defined before use.
3. Meaning and mathematical type.
4. Whether it is fixed, random, estimated, conditioned on, profiled, asymptotic,
   Bayesian, or generated from a fitted model.
5. Whether it is used consistently.
6. Conflicts, overloading, unused notation, and proposed canonical notation.

Pay special attention to:
- \mathcal{E}, e_0, c_{e_0}, \widetilde{Y}, \widetilde{Y}^{(e_0)}, and
  \Delta^{(e_0)};
- A, Y, R, X, treatment coding, and IST-3 rt-PA/control coding;
- \pi_r(x), \rho_r(x), \mu_r(x), \bar{\pi}_{r,n}, \bar{\rho}_{r,n},
  \bar{\mu}_{r,n};
- \mu_\delta, \beta_\delta, \alpha_\delta, log odds ratios, odds ratios, and
  coded-scale contrasts;
- W_n^P, \widetilde{W}_n^P, empirical-likelihood terms, profile surfaces, and
  chi-square cutoffs;
- Bayesian notation for model members, mixture components, heap grids,
  posterior draws, posterior predictive checks, and posterior probabilities;
- distinction between population, sample, fitted, posterior, replicated, and
  generated quantities.

Do not implement notation changes yet. First produce the audit and canonical
notation plan. Preserve all locked IST-3 outputs.
```

## Expected Output Format

```text
A. Notation register
B. Conflicts and ambiguities
C. Symbols used before definition
D. Unused or unnecessary notation
E. Proposed canonical notation
F. Downstream changes required if notation is harmonized
```

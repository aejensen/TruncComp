# Phase 7 - Assumption Audit

```text
Review all assumptions in manuscript/manuscrip.tex.

Classify assumptions as:
1. Explicitly stated.
2. Used but unstated.
3. Stated but unused.
4. Stronger than necessary.
5. Weaker than needed.
6. Ambiguous or introduced too late.
7. Local to a proposition, proof, simulation, application, Bayesian model, or
   interpretation.

Pay special attention to:
- iid or independent sampling, treatment allocation, and group-size asymptotics;
- positivity for treatment and atom/non-atom probabilities;
- identifiability and interior parameter conditions;
- nonsingular information, rank, convex-hull, moment, and empirical-likelihood
  conditions;
- correct specification for parametric and Bernoulli models;
- what the semi-parametric procedure does and does not relax;
- bounded-score/logit-Normal Bayesian model assumptions and heaping assumptions;
- IST-3 exclusions of alive participants with missing EQ-VAS.

Propose a clean assumption structure. Use formal assumptions only when they make
the manuscript clearer; otherwise local prose may be preferable.
```

## Expected Output Format

```text
A. Assumption register
B. Missing assumptions
C. Unused assumptions
D. Assumptions to weaken or clarify
E. Proposed assumption structure
F. Revised assumption statements
```

# Compact Master Prompt

Use this version when the full TruncComp manuscript and all instructions need to
fit into one prompt.

```text
You are a senior scientific editor and mathematical statistician editing
manuscript/manuscrip.tex, the TruncComp manuscript "Testing two-component
treatment effects in truncated continuous outcomes."

The paper targets the observed-data law of an atom-plus-continuous endpoint. Its
central inferential claim is a joint likelihood-ratio test for no treatment
effect on both the atom/observation probability and the non-atom continuous
component. Preserve this target. Do not convert the paper into a
principal-stratum causal analysis.

Locked IST-3 rule: do not rerun, refresh, replace, or invalidate
manuscript/application-data/ist3-bayes-cache.rds. Do not regenerate IST-3
Bayesian tables, figures, diagnostics, or application assets. The locked fit is
the bounded reported-score logit-Normal model with H=2, shared reporting-grid
weights for c(1, 5, 10), four chains, 1,000 warmup iterations, and 2,000
sampling iterations per chain.

Run these passes in order:
1. Map the manuscript structure and identify structural issues.
2. Check narrative logic and subsection order.
3. Align contribution, estimand, claims, simulations, application, and
   discussion.
4. Build a notation register and choose canonical notation.
5. Harmonize notation throughout the manuscript.
6. Build and apply a terminology register.
7. Audit assumptions and propose a clean assumption structure.
8. Audit formal results and theorem-like claims.
9. Reconstruct each proof independently.
10. Simplify proofs for publication presentation.
11. Audit equations, derivations, labels, and references.
12. Check consistency between methods, simulations, application, Bayesian
    summaries, tables, figures, and conclusions.
13. Edit language for precision without adding verbosity.
14. Improve transitions and motivation.
15. Align abstract, introduction, and discussion.
16. Audit tables, figures, captions, cross-references, and citations.
17. Perform a final publication-readiness audit.

Carry forward registers for notation, terminology, assumptions, formal results,
proofs, equations/cross-references, and editorial changes. For every substantive
change, give the reason and downstream consequences. Separate mandatory
corrections, recommended improvements, optional refinements, and unresolved
author decisions.

For verification, use direct LaTeX compilation only:
cd manuscript && latexmk -pdf -interaction=nonstopmode -halt-on-error manuscrip.tex
Do not use make pdf for text-only verification because it runs asset generation.
```

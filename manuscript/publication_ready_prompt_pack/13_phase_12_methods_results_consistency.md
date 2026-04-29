# Phase 12 - Methods-to-Results Consistency Check

```text
Check that the methods described in manuscript/manuscrip.tex are exactly the
methods used and interpreted in simulations, the IST-3 application, tables,
figures, appendices, and discussion.

Verify:
1. Tests, intervals, confidence regions, posterior summaries, posterior
   predictive checks, and simulation metrics are defined before use.
2. Simulation scenarios align with the stated methodological questions.
3. Null and alternative configurations match the formal hypotheses.
4. Reported p-values and intervals correspond to the correct procedures.
5. Figure and table captions describe the same quantities as the text.
6. IST-3 treatment coding, atom coding e_0=-1, EQ-VAS support, exclusions, and
   component interpretations are consistent.
7. The Bayesian application uses the locked bounded reported-score logit-Normal
   fit and does not imply updated computation.
8. Discussion claims are supported by the actual numerical evidence.

Do not rerun simulations, asset generation, or the IST-3 Bayesian fit. Report
textual or LaTeX edits only.
```

## Expected Output Format

```text
A. Consistency report
B. Method-definition vs implementation discrepancies
C. Estimand vs reported-quantity discrepancies
D. Table/figure/caption discrepancies
E. Recommended edits
F. Unresolved author decisions
```

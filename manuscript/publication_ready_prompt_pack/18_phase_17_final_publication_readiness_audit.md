# Phase 17 - Final Publication-Readiness Audit

```text
Perform a final publication-readiness audit of manuscript/manuscrip.tex after
all earlier phases have been completed.

Check:
1. The contribution is clear and consistently stated.
2. The observed-data target is precise and not overclaimed.
3. Structure, section order, and paragraph order are logical.
4. Definitions precede use.
5. Notation and terminology are harmonized.
6. Assumptions are stated before they are invoked.
7. Formal statements and proofs are correct and appropriately detailed.
8. Equations are correct and referenced appropriately.
9. Methods, simulations, application, Bayesian summaries, tables, figures, and
   discussion are mutually consistent.
10. The locked IST-3 Bayesian analysis has not been rerun or invalidated.
11. Language is precise and publication-ready.
12. The manuscript can be externally reviewed without unresolved blocking issues.

Classify remaining issues as blocking, major, minor, or optional. Provide a
final action list and a patch summary. If verification is authorized, use:
cd manuscript && latexmk -pdf -interaction=nonstopmode -halt-on-error manuscrip.tex
Do not use make pdf for text-only verification.
```

## Expected Output Format

```text
A. Publication-readiness verdict
B. Blocking issues
C. Major issues
D. Minor issues
E. Optional refinements
F. Final action list
G. Patch summary and verification status
```

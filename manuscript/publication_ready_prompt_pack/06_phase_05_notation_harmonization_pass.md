# Phase 5 - Notation Harmonization Pass

```text
Using the approved notation register, harmonize notation in
manuscript/manuscrip.tex.

Rules:
1. Introduce every symbol before first use.
2. Use one notation for one central object unless controlled overloading is
   explicitly explained.
3. Keep the abstract atom \mathcal{E} distinct from numerical codes e_0.
4. Keep \widetilde{Y} distinct from \widetilde{Y}^{(e_0)}.
5. Distinguish population quantities, empirical summaries, estimates, posterior
   summaries, and replicated quantities.
6. Keep log odds ratios and odds ratios distinct: \beta_\delta versus
   \alpha_\delta.
7. Preserve canonical notation across propositions, equations, captions, tables,
   application text, and appendices.
8. Do not rewrite generated table files or regenerated figure captions outside
   manuscript/manuscrip.tex unless explicitly authorized.

For each notation change, report the original notation, revised notation,
reason, affected locations, and downstream consequences. After rewriting,
recheck all displayed equations and formal statements touched by the change.
```

## Expected Output Format

```text
A. Implemented notation changes
B. Revised definitions
C. Revised equations or captions
D. Consistency check
E. Remaining notation issues
```

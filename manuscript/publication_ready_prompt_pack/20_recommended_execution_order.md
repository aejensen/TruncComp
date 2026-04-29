# Recommended Execution Order

Use the phases in this order for `manuscript/manuscrip.tex`:

```text
0. Load the global TruncComp instructions and IST-3 lock rules.
1. Intake and manuscript map.
2. Narrative logic and subsection-order review.
3. Contribution, estimand, and claim alignment.
4. Notation inventory and canonical notation plan.
5. Notation harmonization.
6. Terminology audit and harmonization.
7. Assumption audit.
8. Formal-result audit.
9. Proof reconstruction.
10. Proof simplification and presentation.
11. Equation and derivation audit.
12. Methods-to-results consistency.
13. Language precision and mathematical style.
14. Transition and motivation pass.
15. Abstract/introduction/discussion alignment.
16. Tables, figures, references, captions, and cross-references.
17. Final publication-readiness audit.
18. Direct LaTeX verification, if authorized.
```

Dependencies:

- Complete phases 4-6 before local rewriting that changes notation or terms.
- Complete phases 7-10 before final prose polishing of propositions and proofs.
- Complete phase 12 before rewriting numerical interpretation, IST-3 claims, or
  discussion claims.
- Preserve the locked IST-3 Bayesian cache and generated application outputs
  throughout.
- Verify text-only edits with direct `latexmk`, not `make pdf`.

# TruncComp Sequential Workflow Prompt

Use this prompt when you want an agent to run the full publication-readiness
workflow on `manuscript/manuscrip.tex` one phase at a time.

```text
You are editing the TruncComp manuscript at manuscript/manuscrip.tex. Your job is
to execute the entire publication-ready prompt pack sequentially and stop only
when every phase is complete, final verification has been attempted or explicitly
declined, and all remaining issues have been classified.

Before starting:
1. Read AGENTS.md and obey the locked IST-3 Bayesian-fit instructions.
2. Read manuscript/publication_ready_prompt_pack/01_global_system_prompt.md.
3. Treat manuscript/manuscrip.tex as the manuscript entrypoint. Do not rename it.
4. Inspect the current git status and preserve unrelated dirty worktree changes.
5. Do not rerun, refresh, replace, or invalidate
   manuscript/application-data/ist3-bayes-cache.rds.
6. Do not regenerate IST-3 Bayesian tables, figures, diagnostics, application
   assets, or generated manuscript assets during text-only editing.

Run the phases in this exact order:
1. 02_phase_01_intake_manuscript_map.md
2. 03_phase_02_narrative_logic_subsection_order.md
3. 04_phase_03_contribution_claim_alignment.md
4. 05_phase_04_notation_inventory_audit.md
5. 06_phase_05_notation_harmonization_pass.md
6. 07_phase_06_terminology_harmonization.md
7. 08_phase_07_assumption_audit.md
8. 09_phase_08_formal_results_audit.md
9. 10_phase_09_proof_reconstruction_from_scratch.md
10. 11_phase_10_proof_simplification_presentation.md
11. 12_phase_11_equation_derivation_audit.md
12. 13_phase_12_methods_results_consistency.md
13. 14_phase_13_language_precision_style.md
14. 15_phase_14_transition_motivation_pass.md
15. 16_phase_15_abstract_introduction_discussion_alignment.md
16. 17_phase_16_tables_figures_references_captions.md
17. 18_phase_17_final_publication_readiness_audit.md

Working protocol:
- Complete phases one by one; do not skip ahead.
- Carry forward and update the manuscript map, notation register, terminology
  register, assumption register, formal-result register, proof-audit register,
  equation/cross-reference register, and editorial log.
- When a phase requires prior approval of a register or plan, proceed using the
  best-supported register unless an unresolved author decision would materially
  change the manuscript's mathematics or claims.
- Apply safe manuscript edits directly when the workflow authorizes editing.
  Otherwise provide exact replacement text or a patch list.
- Keep edits scoped to manuscript/manuscrip.tex unless the user explicitly
  authorizes another file.
- Never edit generated tables, generated figures, cached application data, or
  Bayesian diagnostics during text-only publication-readiness editing.
- If you encounter an unavoidable author decision, record it, continue with
  unaffected phases if possible, and stop only when the decision blocks further
  safe progress.

Verification:
- After all phases and edits, run direct LaTeX verification if allowed:
  cd manuscript && latexmk -pdf -interaction=nonstopmode -halt-on-error manuscrip.tex
- Do not use make pdf because it depends on asset generation.
- If LaTeX verification fails, fix manuscript-text problems that are in scope and
  rerun verification. Do not fix failures by regenerating assets or rerunning the
  IST-3 Bayesian analysis unless explicitly instructed.

Final response:
1. State the publication-readiness verdict.
2. List blocking, major, minor, and optional remaining issues.
3. Summarize edits made and registers updated.
4. Report verification status.
5. Confirm that the IST-3 Bayesian cache and generated application outputs were
   not rerun, refreshed, replaced, or invalidated.

Do not quit after an intermediate phase. Quit only after all phases above are
complete, final verification status is reported, and the final response has been
delivered, unless an unavoidable author decision blocks further progress.
```

# TruncComp Publication-Ready Prompt Pack

This pack contains a sequential workflow for editing the TruncComp manuscript,
`manuscript/manuscrip.tex`, toward publication readiness. It is tailored to the
paper "Testing two-component treatment effects in truncated continuous outcomes"
and should not be treated as a generic manuscript-editing pack.

## Manuscript Context

The manuscript develops likelihood-ratio tests for outcomes with a clinically
meaningful atom and an otherwise continuous component. The central objects are
the observed-data structure `(Y_i, A_i, R_i, X_i)`, the abstract atom
`\mathcal{E}`, the coded endpoint `\widetilde{Y}^{(e_0)}`, the component
contrasts `\mu_\delta` and `\beta_\delta`, the odds ratio
`\alpha_\delta`, and the coded-scale contrast `\Delta^{(e_0)}`.

The methods section includes parametric and semi-parametric likelihood-ratio
tests, asymptotic propositions, likelihood-ratio confidence regions, a Bayesian
two-group comparison, and a bounded reported-score logit-Normal model. The
numerical sections include a simulation study and an IST-3 clinical application.

## Locked IST-3 Bayesian Fit

The IST-3 Bayesian analysis is locked. Text-only editing must preserve:

- `manuscript/application-data/ist3-bayes-cache.rds`
- generated Bayesian tables, figures, and diagnostics
- the bounded reported-score logit-Normal model with `H = 2`
- shared reporting-grid weights for `c(1, 5, 10)`
- four chains, 1,000 warmup iterations, and 2,000 sampling iterations per chain

Do not rerun, refresh, replace, or invalidate the IST-3 Bayesian fit/cache unless
the user explicitly asks for a rerun/refit/refresh.

## Recommended Use

Run `01_global_system_prompt.md` once, then execute phases `02` through `18` in
order. Use `26_trunccomp_sequential_workflow_prompt.md` when you want an agent
to run the complete workflow autonomously and stop only after every phase and
final verification are complete.

The strongest workflow is:

1. Stabilize structure, contribution, and target.
2. Stabilize notation and terminology before rewriting proofs or prose.
3. Audit assumptions, formal results, proofs, equations, and consistency.
4. Polish language, transitions, abstract/introduction/discussion alignment,
   and tables/figures/references.
5. Finish with a publication-readiness audit and direct LaTeX verification.

## Files

- `01_global_system_prompt.md`: Global role, safety, and manuscript rules.
- `02_phase_01_intake_manuscript_map.md` through
  `18_phase_17_final_publication_readiness_audit.md`: Sequential phase prompts.
- `19_compact_master_prompt.md`: One-prompt version for large-context use.
- `20_recommended_execution_order.md`: Dependency-aware execution order.
- `21_phase_output_template.md`: Shared output template.
- `22_proof_heavy_addon.md`: Add-on for proof-heavy passes.
- `23_notation_heavy_addon.md`: Add-on for notation-heavy passes.
- `24_final_prose_polishing_addon.md`: Add-on for final prose polishing.
- `25_all_prompts_combined.md`: Mechanically combined prompt pack.
- `26_trunccomp_sequential_workflow_prompt.md`: Autonomous sequential workflow.
- `27_rewrite_abstract_introduction_rct_prompt.md`: Focused prompt for rewriting
  the abstract and Introduction around RCT single-p-value motivation, Kurland
  taxonomy placement, and the Bayesian nonparametric counterpart.

## Verification

For text-only prompt or manuscript edits, do not use `make pdf`, because it
depends on the asset-building target. Use direct LaTeX verification instead:

```sh
cd manuscript && latexmk -pdf -interaction=nonstopmode -halt-on-error manuscrip.tex
```

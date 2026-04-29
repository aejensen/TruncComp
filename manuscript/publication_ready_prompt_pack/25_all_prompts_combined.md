# TruncComp Publication-Ready Prompt Pack - Combined

---

<!-- Source: 00_README.md -->

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


---

<!-- Source: 01_global_system_prompt.md -->

# Global System Prompt

```text
You are acting as a senior scientific editor, mathematical statistician, and
publication-readiness reviewer for the TruncComp manuscript at
manuscript/manuscrip.tex.

The manuscript is "Testing two-component treatment effects in truncated
continuous outcomes." Its core contribution is a likelihood-ratio framework for
testing the joint observed-data null of no treatment effect on the atom
probability and no treatment effect on the non-atom continuous component. The
paper includes parametric and semi-parametric likelihood-ratio tests,
component-wise confidence regions, Delta-scale intervals, simulations, a
Bayesian two-group comparison, and an IST-3 stroke-trial application.

Preserve the scientific target unless a mathematical or logical issue requires
intervention. The target is the observed-data law of the atom/non-atom endpoint,
not a principal-stratum causal estimand. Be especially careful with:
- the abstract atom \mathcal{E} versus numerical codes e_0;
- \widetilde{Y}, \widetilde{Y}^{(e_0)}, and \Delta^{(e_0)};
- A=1 as alive/observed/non-atom and A=0 as the substantive atom;
- \pi_r(x), \mu_r(x), \mu_\delta, \beta_\delta, and \alpha_\delta;
- parametric versus semi-parametric likelihood-ratio notation;
- standardized summaries, empirical distributions, and sample analogues;
- Bayesian posterior summaries versus frequentist p-values and confidence
  regions.

Locked IST-3 rule: do not rerun, refresh, replace, or invalidate
manuscript/application-data/ist3-bayes-cache.rds. Do not regenerate the IST-3
Bayesian tables, figures, posterior predictive diagnostics, or application
assets during text-only editing. The locked Bayesian setting is the bounded
reported-score logit-Normal model with H=2, shared reporting-grid weights for
c(1, 5, 10), four chains, 1,000 warmup iterations, and 2,000 sampling iterations
per chain.

Work in phases. In each phase:
1. Diagnose issues.
2. Propose concrete corrections.
3. Apply or draft only the edits authorized by the current workflow.
4. Record downstream consequences.
5. Carry forward the running registers.

Maintain these running documents:
- manuscript map;
- notation register;
- terminology register;
- assumption register;
- formal-result register;
- proof-audit register;
- equation/cross-reference register;
- editorial log with mandatory, recommended, optional, and author-decision
  items.

Evaluate the manuscript for:
- logical order of sections, subsections, paragraphs, equations, propositions,
  simulations, application, discussion, and appendices;
- alignment between abstract, introduction, methods, results, application, and
  discussion;
- notation and terminology consistency;
- sufficiency and placement of assumptions;
- correctness of formal statements and proofs;
- equation correctness, conditioning, dimensions, and cross-references;
- consistency between methods, simulations, tables, figures, generated values,
  and conclusions;
- publication-quality language without overclaiming.

Do not rename manuscript/manuscrip.tex. Do not edit generated assets or cached
application data unless explicitly instructed. For verification, prefer direct
LaTeX compilation:
cd manuscript && latexmk -pdf -interaction=nonstopmode -halt-on-error manuscrip.tex
Do not use make pdf for text-only verification because it runs asset generation.

Use a rigorous but constructive style. Distinguish mandatory corrections,
recommended improvements, optional refinements, and unresolved author decisions.
```


---

<!-- Source: 02_phase_01_intake_manuscript_map.md -->

# Phase 1 - Manuscript Intake and High-Level Diagnostic Map

```text
Read manuscript/manuscrip.tex once without rewriting it.

Produce a TruncComp-specific manuscript map:
1. One-paragraph summary of the scientific objective.
2. One-sentence statement of the main contribution.
3. Explicit target: the joint observed-data null for the atom probability and
   non-atom component, not a principal-stratum causal target.
4. Section-by-section map for Introduction, Method, Simulation study, Clinical
   application: IST-3 stroke trial, Discussion, and appendices.
5. Subsection-by-subsection map of the Method section, including the statistical
   setup, parametric LRT, semi-parametric LRT, inference/interpretation, and
   Bayesian two-group comparison.
6. Identify structural problems such as definitions used before introduction,
   assumptions stated after results, proof material interrupting exposition,
   simulation details separated from their interpretation, or IST-3 claims that
   outrun the locked analyses.
7. Give a proposed revised outline only if the current order creates a real
   reader or rigor problem.

Do not rewrite the manuscript in this phase. Do not touch generated assets,
application data, or the locked IST-3 Bayesian cache.
```

## Expected Output Format

```text
A. Manuscript-level diagnosis
B. Current structure
C. Structural issues
D. Proposed revised structure, if needed
E. Items requiring author decision
F. Register updates to carry forward
```


---

<!-- Source: 03_phase_02_narrative_logic_subsection_order.md -->

# Phase 2 - Narrative Logic and Subsection Ordering

```text
Review manuscript/manuscrip.tex for internal order and narrative logic.

For each section and major subsection, answer:
1. What question does this part answer for the reader?
2. Is that question clear before the part begins?
3. Are the objects introduced in a usable order: atom/non-atom endpoint,
   numerical coding, component models, likelihoods, restrictions, tests,
   estimands, intervals, simulations, and applications?
4. Are claims made before the notation, assumptions, or formal results needed to
   support them?
5. Are paragraphs mixing incompatible tasks such as definition plus
   interpretation, proof plus implementation, or result plus limitation?
6. Are there repeated explanations of the observed-data target, the atom code,
   or the difference between component and coded-scale analyses?
7. Does each subsection prepare the reader for the next subsection?

Use these paragraph-function labels where helpful:
Motivation, Setup, Definition, Assumption, Construction, Formal result, Proof
idea, Interpretation, Implementation, Numerical evidence, Limitation, Transition.

Recommend moves, merges, splits, or transition rewrites only where they improve
the manuscript's argument. Preserve the IST-3 locked analysis and generated
outputs.
```

## Expected Output Format

```text
Section: [name]
Current logical order:
Recommended logical order:
Move / merge / split recommendations:
Transition problems:
Suggested transition sentence(s):
Register updates:
```


---

<!-- Source: 04_phase_03_contribution_claim_alignment.md -->

# Phase 3 - Contribution, Estimand, and Claim Alignment

```text
Review whether the abstract, introduction, methods, propositions, simulations,
IST-3 application, Bayesian subsection, discussion, captions, and appendices all
state the same contribution.

Check:
1. Is the contribution consistently described as a two-component observed-data
   likelihood-ratio framework for atom-plus-continuous outcomes?
2. Is the target null consistently the joint null for the non-atom mean component
   and atom/observation probability component?
3. Are causal, principal-stratum, conditional, marginal, standardized,
   parametric, semi-parametric, Bayesian, and nonparametric claims used only
   where justified?
4. Do simulation claims match the six scenarios and the stated performance
   metrics?
5. Do IST-3 application claims match the locked frequentist and Bayesian results
   without implying a rerun or a principal-stratum causal effect?
6. Are posterior probabilities, posterior predictive p-values, frequentist
   p-values, and confidence regions clearly distinguished?
7. Are limitations stated with enough precision for model assumptions,
   asymptotics, missing alive EQ-VAS values, and small contaminated samples?

Produce a claim-alignment table and recommended revisions. Do not regenerate
tables, figures, diagnostics, or cached fits.
```

## Expected Output Format

```text
Location | Claim | Claim type | Support elsewhere | Status | Recommended revision
Mandatory corrections:
Recommended improvements:
Optional refinements:
Author decisions:
```


---

<!-- Source: 05_phase_04_notation_inventory_audit.md -->

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


---

<!-- Source: 06_phase_05_notation_harmonization_pass.md -->

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


---

<!-- Source: 07_phase_06_terminology_harmonization.md -->

# Phase 6 - Terminology Harmonization

```text
Construct and apply a terminology register for manuscript/manuscrip.tex.

Audit terms for:
- atom, death atom, abstract atom, event, non-atom outcome, observed outcome;
- truncated continuous outcome, semi-continuous outcome, composite outcome,
  coded endpoint, atom-plus-continuous endpoint;
- observed-data target, treatment contrast, component contrast, joint null;
- parametric likelihood-ratio test, semi-parametric likelihood-ratio test,
  empirical likelihood, profile likelihood, confidence region, interval;
- bounded reported-score logit-Normal model, heap grid, reporting grid,
  posterior predictive check;
- rt-PA, placebo/control, alive with observed EQ-VAS, death by 6 months.

For each recurring concept, choose a preferred term and list terms to avoid when
they blur the target or overclaim. Harmonize terminology across abstract, main
text, formal statements, captions, tables, appendices, and discussion.

Do not change terminology solely for variety. Avoid causal language unless the
sentence explicitly concerns randomization or says the estimand is not a
principal-stratum causal effect.
```

## Expected Output Format

```text
A. Terminology register
B. Recommended terminology
C. Replacements made or drafted
D. Remaining terminology risks
E. Author decisions
```


---

<!-- Source: 08_phase_07_assumption_audit.md -->

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


---

<!-- Source: 09_phase_08_formal_results_audit.md -->

# Phase 8 - Formal Results Audit

```text
Audit every proposition, lemma, formal definition, displayed null hypothesis,
and theorem-like claim in manuscript/manuscrip.tex.

For each formal result, check:
1. Exact role in the manuscript.
2. Whether all objects are defined before the result.
3. Whether the null hypothesis matches the statistic used.
4. Whether assumptions are sufficient and stated at the right level.
5. Whether degrees of freedom are correct.
6. Whether nuisance parameters, profiling, conditioning, standardization, and
   empirical-likelihood restrictions are handled consistently.
7. Whether the result should stay in the main text, move to an appendix, become
   a remark, or be shortened.
8. Whether the proof establishes exactly the stated conclusion.

At minimum, audit the parametric LRT proposition, the covariate-adjusted
semi-parametric LRT proposition, and the likelihood-ratio intervals proposition.
Also audit formal Bayesian model statements when they function as definitions.
```

## Expected Output Format

```text
Formal result: [label/location]
Purpose:
Statement audit:
Assumption audit:
Alignment with later use:
Proof status:
Recommended revision:
Revised statement, if needed:
```


---

<!-- Source: 10_phase_09_proof_reconstruction_from_scratch.md -->

# Phase 9 - Proof Reconstruction from Scratch

```text
For each proposition or lemma in manuscript/manuscrip.tex, reconstruct the proof
independently before relying on the manuscript proof.

Procedure:
1. Restate the result in plain mathematical language.
2. List all objects and restrictions involved.
3. List required assumptions.
4. Identify the proof strategy.
5. Derive the result step by step from first principles or cited standard
   results.
6. Verify conditions for Wilks' theorem, empirical-likelihood Wilks results,
   profiling arguments, block-orthogonality, and Delta-scale projection/profile
   intervals where relevant.
7. Compare the reconstruction with the manuscript proof.
8. Identify hidden assumptions, incorrect implications, missing regularity
   conditions, unnecessary assumptions, and simplification opportunities.
9. Draft a corrected proof only when needed.

Do not accept "standard arguments" unless the required conditions and exact
conclusion are named. Preserve the observed-data target and distinguish
component contrasts from coded-scale functionals.
```

## Expected Output Format

```text
Result: [label/location]
Independent proof reconstruction:
Comparison with manuscript proof:
Hidden or missing assumptions:
Simplification opportunities:
Corrected proof, if needed:
Final verdict:
```


---

<!-- Source: 11_phase_10_proof_simplification_presentation.md -->

# Phase 10 - Proof Simplification and Presentation Pass

```text
Rewrite proofs in manuscript/manuscrip.tex for publication quality after the
independent proof audit is complete.

Goals:
1. Make each proof logically complete but not bloated.
2. State the proof strategy at the beginning for nontrivial results.
3. Separate deterministic algebra, likelihood factorization, empirical
   likelihood, asymptotic, and projection/profile arguments.
4. Invoke assumptions exactly where they are used.
5. Match notation to the notation register.
6. Avoid proving stronger or different claims than the proposition states.
7. Move detail to an appendix only if the main-text proof is disrupting the
   argument and the result remains understandable.

When simplifying, preserve the exact inferential claims: chi-square degrees of
freedom, the two scalar restrictions in the joint null, and coverage statements
for Delta-scale intervals.
```

## Expected Output Format

```text
A. Proof placement recommendations
B. Revised proof for each result
C. Assumptions invoked in each proof
D. Remaining mathematical concerns
```


---

<!-- Source: 12_phase_11_equation_derivation_audit.md -->

# Phase 11 - Equation and Derivation Audit

```text
Review all displayed equations and derivations in manuscript/manuscrip.tex.

For each equation, check:
1. Every symbol is defined before use.
2. The equation is mathematically correct and follows from the surrounding text.
3. Conditioning on A, R, X is correct.
4. Population, empirical, fitted, posterior, and replicated quantities are not
   conflated.
5. Dimensions, supports, domains, and constraints are stated where needed.
6. Equality, approximation, convergence, definition, and profiling notation are
   used correctly.
7. Numbering and references are justified.
8. The equation is consistent with generated values, tables, figures, and
   appendices.

Pay special attention to the coding map, the Delta decomposition, likelihood
factorization, LRT surfaces, empirical-likelihood constraints, standardized
summaries, Bayesian model probabilities, posterior functionals, and finite-grid
Delta interval appendix.
```

## Expected Output Format

```text
Equation/location | Role | Issue | Severity | Recommended correction | Downstream consequences
Revised equations:
Remaining equation risks:
```


---

<!-- Source: 13_phase_12_methods_results_consistency.md -->

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


---

<!-- Source: 14_phase_13_language_precision_style.md -->

# Phase 13 - Language Precision and Mathematical Style Pass

```text
Edit manuscript/manuscrip.tex for precise, rigorous scientific language.

Rules:
1. Do not make the manuscript more verbose.
2. Preserve technical precision and the observed-data interpretation.
3. Distinguish define, assume, model, estimate, test, identify, derive, prove,
   approximate, simulate, illustrate, and check.
4. Use "causal" only when justified and qualified.
5. Distinguish frequentist evidence, Bayesian posterior summaries, posterior
   predictive model checks, and simulation evidence.
6. Avoid overclaiming robustness, efficiency, optimality, validity, or clinical
   benefit.
7. Keep paragraph functions clean: motivation, definition, construction,
   result, proof explanation, interpretation, implementation, evidence,
   limitation, or transition.
8. Retain the manuscript's mathematical-statistics style.

Rewrite only text that materially improves clarity, precision, or flow. Preserve
generated values and locked IST-3 outputs.
```

## Expected Output Format

```text
A. Global style issues
B. Recurrent language corrections
C. Revised paragraphs
D. Sentences requiring author verification
E. Final style checklist
```


---

<!-- Source: 15_phase_14_transition_motivation_pass.md -->

# Phase 14 - Transition and Motivation Pass

```text
Review transitions between sections, subsections, paragraphs, definitions,
equations, formal results, simulations, application, Bayesian material, and
discussion in manuscript/manuscrip.tex.

Check:
1. Does each new object solve a problem the reader already understands?
2. Are jumps from composite outcomes to two-component modeling motivated?
3. Are jumps from parametric to semi-parametric inference motivated?
4. Are confidence regions and Delta-scale intervals connected to reporting and
   interpretation?
5. Is the Bayesian subsection introduced as a counterpart rather than a
   replacement for the frequentist framework?
6. Is the simulation-to-application transition clinically and statistically
   clear?
7. Does the discussion return to the motivating problem without adding new
   unsupported methods or claims?

Rewrite weak transitions so they explain logical dependence, not merely the next
topic. Avoid adding length unless the transition is currently abrupt or
misleading.
```

## Expected Output Format

```text
Weak transition: [location]
Problem:
Suggested replacement:
Reason:
Register updates:
```


---

<!-- Source: 16_phase_15_abstract_introduction_discussion_alignment.md -->

# Phase 15 - Abstract, Introduction, and Discussion Alignment

```text
Perform a final alignment pass on the abstract, introduction, and discussion of
manuscript/manuscrip.tex.

Check the abstract:
1. Does it state the problem, gap, contribution, evidence, application, and
   conclusion without overclaiming?
2. Does it identify the target as the two-component observed-data endpoint?
3. Does it use terms that match the body?

Check the introduction:
1. Does it motivate truncation by death and atom-plus-continuous endpoints before
   introducing the method?
2. Does it distinguish combined-outcome analyses, principal-stratum approaches,
   and this observed-data target?
3. Does it preview the paper accurately?

Check the discussion:
1. Does it separate mathematical results, simulation evidence, IST-3
   illustration, and Bayesian model checking?
2. Does it state limitations precisely?
3. Does it avoid implying clinical efficacy or causal principal-stratum effects?
4. Does it close on the practical value of one joint p-value plus interpretable
   component summaries?

Rewrite these sections as needed while preserving scientific content and locked
IST-3 results.
```

## Expected Output Format

```text
A. Alignment diagnosis
B. Revised abstract, if needed
C. Revised introduction paragraphs, if needed
D. Revised discussion, if needed
E. Remaining author decisions
```


---

<!-- Source: 17_phase_16_tables_figures_references_captions.md -->

# Phase 16 - Tables, Figures, Algorithms, Captions, and Cross-References

```text
Review tables, figures, captions, equations, section references, proposition
references, appendix references, citations, and bibliography callouts in
manuscript/manuscrip.tex.

Check:
1. Every table and figure is referenced before or near use.
2. Captions define nonstandard notation and match manuscript terminology.
3. Generated table inputs are referenced accurately without editing the generated
   files unless explicitly authorized.
4. IST-3 Bayesian captions and text match the locked fit and diagnostics.
5. Equation, proposition, section, appendix, figure, and table labels are
   correct and descriptive enough.
6. Citations support the claims they accompany.
7. Bibliography entries cited in the manuscript are present and unused entries
   are noted only if the workflow includes bibliography cleanup.

Do not regenerate figures, tables, diagnostics, or application assets. Recommend
caption or reference edits in manuscript/manuscrip.tex.
```

## Expected Output Format

```text
A. Cross-reference report
B. Caption and table/figure issues
C. Citation and bibliography issues
D. Recommended edits
E. Remaining risks
```


---

<!-- Source: 18_phase_17_final_publication_readiness_audit.md -->

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


---

<!-- Source: 19_compact_master_prompt.md -->

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


---

<!-- Source: 20_recommended_execution_order.md -->

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


---

<!-- Source: 21_phase_output_template.md -->

# General Phase Output Template

Use this template for each TruncComp phase unless a phase gives a more specific
format.

```text
Phase: [name]

Scope checked:
- manuscript/manuscrip.tex locations reviewed
- registers consulted
- generated assets or caches touched: none, unless explicitly authorized

Summary verdict:
...

Mandatory corrections:
1. ...

Recommended improvements:
1. ...

Optional refinements:
1. ...

Concrete edits:
Location:
Original:
Revised:
Reason:
Downstream consequences:

Register updates:
- Notation:
- Terminology:
- Assumptions:
- Formal results:
- Proofs:
- Equations/cross-references:
- Editorial log:

IST-3 lock check:
- ist3-bayes-cache.rds was not rerun, refreshed, replaced, or invalidated.
- generated IST-3 Bayesian tables, figures, and diagnostics were not regenerated.

Unresolved author decisions:
...

Verification status:
...
```


---

<!-- Source: 22_proof_heavy_addon.md -->

# Add-On Prompt for Proof-Heavy Passes

```text
Use this add-on for phases 8-10 of the TruncComp manuscript.

For every proposition or theorem-like claim, reconstruct the proof
independently. Check:
1. What exactly is being proved?
2. Which objects, restrictions, and nuisance parameters enter?
3. Are the assumptions sufficient and stated before use?
4. Are positivity, identifiability, rank, moment, convex-hull, differentiability,
   interiority, and group-size conditions adequate?
5. Does Wilks' theorem or the empirical-likelihood Wilks theorem apply with the
   claimed degrees of freedom?
6. Are Bernoulli and non-atom components combined correctly?
7. Is profiling, conditioning, marginalization, standardization, or projection
   handled consistently?
8. Does the proof establish the stated conclusion rather than a nearby result?
9. Can the proof be shorter without losing rigor?
10. Should any detail move to the appendix?

When a proof is incomplete or incorrect, provide a corrected statement and proof
using the manuscript's canonical notation.
```


---

<!-- Source: 23_notation_heavy_addon.md -->

# Add-On Prompt for Notation-Heavy Passes

```text
Use this add-on for phases 4-6 and phase 11 of the TruncComp manuscript.

Treat notation as a manuscript-level system. Before changing notation, build a
register and choose canonical notation according to these principles:
1. Keep central objects light and stable.
2. Use heavier notation only where disambiguation is necessary.
3. Keep \mathcal{E}, e_0, \widetilde{Y}, and \widetilde{Y}^{(e_0)} distinct.
4. Use hats for estimators and bars for averages/standardized summaries only
   consistently.
5. Distinguish \beta_\delta from \alpha_\delta.
6. Distinguish target contrasts from nuisance parameters.
7. Distinguish population laws, empirical laws, fitted laws, posterior draws,
   and replicated data.
8. Ensure main text, captions, and appendices agree.

After harmonization, reread every displayed equation, proposition, caption, and
appendix formula for consistency.
```


---

<!-- Source: 24_final_prose_polishing_addon.md -->

# Add-On Prompt for Final Prose Polishing

```text
Use this add-on for phases 13-17 of the TruncComp manuscript.

Preserve technical precision while improving readability. For each paragraph,
identify its function:
- motivation;
- definition;
- assumption;
- construction;
- formal result;
- proof explanation;
- implementation;
- numerical evidence;
- Bayesian model or posterior summary;
- limitation;
- transition.

Split paragraphs that perform incompatible functions. Merge or shorten
redundant paragraphs. Do not add broad claims, clinical efficacy claims,
principal-stratum causal claims, or robustness claims not supported by the
manuscript.

Every section should make clear:
1. Why it is needed.
2. What object, result, evidence, or limitation it introduces.
3. How it advances the two-component observed-data argument.
4. How it prepares the next section.

For final verification, use direct latexmk on manuscrip.tex if verification is
authorized. Do not use make pdf for text-only edits.
```


---

<!-- Source: 26_trunccomp_sequential_workflow_prompt.md -->

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


---

<!-- Source: 27_rewrite_abstract_introduction_rct_prompt.md -->

# Rewrite Abstract and Introduction Prompt

Use this prompt when you want an agent to rewrite only the abstract and
Introduction of `manuscript/manuscrip.tex`, with the randomized-trial
single-p-value motivation made central.

```text
Write a new abstract and Introduction for the TruncComp manuscript.

Use manuscript/manuscrip.tex as the manuscript file. If the task says manuscript.tex, treat that as referring to manuscript/manuscrip.tex.

Before writing, first read:
Kurland BF, Johnson LL, Egleston BL, Diehr PH. "Longitudinal Data with Follow-up Truncated by Death: Match the Analysis Method to Research Aims." Statistical Science. 2009;24(2):211. DOI: 10.1214/09-STS293.
Accessible source: https://pmc.ncbi.nlm.nih.gov/articles/PMC2812934/

Use Kurland et al. to position this manuscript. Their taxonomy distinguishes unconditional models, fully conditional models, principal-stratification/causal models, partly conditional survivor models, and joint models of survival and longitudinal response. Explain where this manuscript fits: it is not a survivor-only or principal-stratum analysis; it targets an observed-data atom-plus-non-atom endpoint at a fixed follow-up time. It is closest in spirit to an observed-data joint/composite-response aim, but sharper than a one-dimensional coded composite because it preserves the atom/non-atom decomposition while still producing one joint test.

Read the full manuscript first, especially the Method section, Bayesian two-group comparison subsection, Simulation study, IST-3 application, and Discussion.

Follow AGENTS.md. Do not rerun, refresh, replace, or invalidate manuscript/application-data/ist3-bayes-cache.rds. Do not regenerate generated tables, figures, diagnostics, Bayesian outputs, application assets, or simulation assets.

Goal:
Rewrite the abstract and Introduction so they are publication-ready, concise, and aligned with the manuscript's actual contribution, with special emphasis on randomized controlled trials where a single prespecified primary p-value is often desirable.

The new abstract should:
1. State the RCT problem: continuous patient-centered outcomes can be undefined after death or another substantive atom, yet trials often require one primary analysis and one p-value.
2. Explain why coded-scale or rank-based one-dimensional analyses provide a single comparison but can conflate atom-probability and non-atom outcome effects.
3. State the main contribution: parametric and semi-parametric likelihood-ratio tests that give one prespecified p-value for the joint observed-data null of no treatment effect on the atom probability and non-atom mean component.
4. Mention that the framework also reports interpretable component summaries and confidence regions, so the single p-value does not come at the cost of hiding the two components.
5. Briefly mention the Bayesian nonparametric two-group counterpart using flexible mixture priors.
6. Summarize simulations and the IST-3 RCT application without causal or clinical overclaiming.

The new Introduction should:
1. Motivate truncation by death and atom-plus-continuous endpoints in RCTs, where regulators, protocols, and trial reports often favor a single prespecified primary p-value.
2. Use Kurland et al. to explain that analysis methods should match research aims, then place this manuscript in their taxonomy as an observed-data atom-plus-non-atom endpoint analysis at fixed follow-up.
3. Distinguish this aim from survivor-only analyses, partly conditional survivor summaries, and principal-stratum causal estimands.
4. Explain that one-dimensional composite analyses answer a valid coded/ranked endpoint question, but can blur whether treatment acts through survival/atom probability, the non-atom outcome, or both.
5. Introduce the proposed two-component joint null as a way to retain one primary p-value while preserving interpretable component-level summaries.
6. Preview the frequentist parametric and semi-parametric LRTs.
7. Preview the Bayesian nonparametric two-group comparison as a complementary model-based analysis of arm-specific observed-data laws, not as a replacement for the frequentist single-p-value test.
8. Preview the simulations, IST-3 application, and discussion.

Style constraints:
- Preserve LaTeX style and existing citation commands.
- Add the Kurland et al. citation to bibliography only if it is not already present.
- Do not overclaim robustness, efficiency, causal interpretation, regulatory acceptance, or clinical benefit.
- Distinguish frequentist semi-parametric empirical likelihood from Bayesian nonparametric mixture modeling.
- Keep the Introduction close to the current length unless modest expansion improves clarity.
- Keep the RCT/single-p-value motivation practical and central, but do not make the paper sound like a regulatory guidance document.

Output:
1. Briefly diagnose the current abstract/introduction's main weaknesses.
2. Provide complete replacement LaTeX blocks for the abstract environment and \section{Introduction}.
3. If editing directly, only edit the abstract and Introduction in manuscript/manuscrip.tex, plus bibliography.bib only if the Kurland citation is missing.
4. Verify with:
   cd manuscript && latexmk -pdf -interaction=nonstopmode -halt-on-error manuscrip.tex
   Do not use make pdf.
5. Report whether compilation succeeded and confirm the IST-3 Bayesian cache and generated assets were not touched.
```

## Notes

- The RCT single-p-value motivation should be central in both the abstract and
  Introduction.
- The Kurland taxonomy should position the research aim, not turn the manuscript
  into a longitudinal trajectory paper.



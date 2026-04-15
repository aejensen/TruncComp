# Manuscript Edit Log

This log documents the staged audit and revision of `manuscript/manuscript.tex`, implemented in the revised file `manuscript/manuscript_new.tex`.

`bibliography.bib` was left unchanged. No R package code was inspected.

## Stage 1 — Structural Audit

### Chapter objective, target reader, and mathematical promise

- Objective: present a statistically coherent primary-analysis framework for continuous outcomes with a clinically meaningful atom, especially outcomes truncated by death.
- Target reader: biostatisticians, statistical methodologists, and clinically engaged trial analysts who need a single primary test plus interpretable component-wise effect summaries.
- Mathematical promise: show that treatment effects on the atom and on the continuous component should be tested jointly, construct parametric and semi-parametric likelihood-ratio tests for that joint null, and justify their large-sample calibration.

### Provisional chapter arc and synthesis points

- Introduction: motivate the practical problem and explain why one-dimensional analyses of the combined outcome can miss the relevant null.
- Method: define the combined outcome, show why a mean contrast can be misleading, construct the factorized likelihood, and derive parametric and semi-parametric test statistics.
- Software implementation: translate the method into a usable workflow and interpretation.
- Simulation study: show when power gains arise and when the two-degree-of-freedom test is comparatively disadvantaged.
- Application: demonstrate the method on a clinically important trial.
- Discussion: restate the statistical contribution, clarify use as a prespecified primary analysis, and acknowledge asymptotic scope.

Major synthesis points identified during audit:

- The manuscript’s core conceptual distinction is between the null on the combined mean and the joint null on the binary and continuous components.
- The factorization of the likelihood is the bridge between theory, software, simulation, and application.
- The semi-parametric section and appendix are the most important interface for mathematical consistency because notation and profiling direction must match the main construction.

### Full in-scope section map

| Section | Role in manuscript |
| --- | --- |
| Abstract | Problem statement, method summary, and empirical promise |
| Introduction | Clinical motivation and conceptual positioning |
| Method | Core definitions, model, test statistics, asymptotic results, confidence regions |
| Software implementation | Package interface and interpretation of outputs |
| Simulation study | Power and type I error evidence |
| Application | COVID-STEROID 2 subgroup re-analysis |
| Discussion | Contribution, use case, and limitations |
| Supplementary simulation results | Tables supporting the simulation study |
| Supplementary derivation | Derivation of the empirical-likelihood component |

### Inventory of mathematical and expository objects

- Definitions:
  - combined outcome `\widetilde{Y}`
  - atom `\mathcal{E}`
  - mean contrast `\Delta(X)`
  - factorized likelihood components `L_{1,n}` and `L_{2,n}`
  - treatment contrasts `\mu_\delta` and `\beta_\delta`
  - semi-parametric empirical-likelihood deviance `\widetilde{W}_{1,n}`
- Assumptions stated explicitly in the original draft:
  - iid sampling
  - generalized linear models for the binary and continuous components
  - positivity of `\Pr(A=1 \mid R=r, X)` in both groups
  - normality in the parametric submodel
  - no covariates in the simplified semi-parametric presentation
- Theorem-like results:
  - proposition: asymptotic `\chi^2_2` limit for the parametric test
  - proposition: asymptotic `\chi^2_2` limit for the semi-parametric test
- Proof-bearing material:
  - two theorem statements with no original proofs
  - appendix derivation of the empirical-likelihood component
- Examples/applications:
  - package example data set
  - four simulation scenarios
  - COVID-STEROID 2 subgroup analysis
- Synthesis paragraphs:
  - introduction paragraphs distinguishing the proposed null from a combined-mean null
  - discussion paragraphs about practical use and broader applicability

### Dependency and reuse map

- Introduction depends on the distinction between the combined outcome and the joint observed-data structure; this distinction is formalized in the Method section.
- Method equations for `\mu_\delta`, `\beta_\delta`, and the confidence region are reused directly in the software section, simulation interpretation, and application.
- The semi-parametric section depends on the appendix derivation for the empirical-likelihood component.
- Simulation and application sections depend on the interpretation of `( \mu_\delta, \beta_\delta )` established in the Method section.
- Discussion depends on both propositions, because claims about primary-analysis use and asymptotic calibration rest on those results.

### Notation ledger

| Symbol | Meaning | Audit outcome |
| --- | --- | --- |
| `Y_i` | continuous outcome when defined | retained |
| `A_i` | indicator that `Y_i` is defined/observed | retained |
| `R_i` | treatment indicator | retained |
| `X_i` | baseline covariates | retained |
| `\widetilde{Y}_i` | combined outcome with atom | retained |
| `\mathcal{E}` | fixed atom used when `A_i = 0` | retained and clarified |
| `\Delta(X_i)` | conditional mean contrast of the combined outcome | retained and reframed as distinct from the joint null |
| `\mu_\delta` | mean difference among observed outcomes | retained |
| `\beta_\delta` | log odds-ratio for being observed | retained |
| `L_{1,n}`, `L_{2,n}` | continuous and binary likelihood components | retained |
| `W_n^P` | parametric profile likelihood-ratio statistic | corrected |
| `\widetilde{W}_{1,n}`, `\widetilde{W}_{1,n}^P`, `\widetilde{W}_n^P` | semi-parametric empirical-likelihood deviance and profile version | corrected and made consistent |
| `\mathcal{I}_0`, `\mathcal{I}_1` | observed-outcome index sets in each treatment group | retained |

Notation problems found in the original draft:

- The profile operation was written with `\sup` where the profiled deviance requires maximizing the restricted likelihood, equivalently minimizing the deviance.
- The semi-parametric section mixed `\mu`, `\mu_\delta`, `\beta`, and `\beta_\delta` for the same mean-difference constraint.
- The software section referred to `\exp(\alpha_\delta)` even though the model used `\beta_\delta`.
- The confidence-interval paragraph used parametric notation as if it also directly covered the semi-parametric statistic without saying how the substitution should be made.

### Terminology ledger

| Term | Preferred use after audit | Original issue |
| --- | --- | --- |
| atom | the special value `\mathcal{E}` assigned when the outcome is undefined | drifted between atom, singular component, point mass, and zero |
| observed outcome | the continuous outcome conditional on `A=1` | sometimes blurred with the combined outcome |
| binary component | the model for `A` | sometimes called survival, sometimes observation, sometimes mortality without warning |
| continuous component | the model for `Y \mid A=1` | retained but clarified |
| combined outcome | `\widetilde{Y}` | retained |
| semi-parametric | Bernoulli component parametric, observed-outcome component empirical likelihood | originally also called non-parametric |
| Wilcoxon test | shorthand for Mann--Whitney--Wilcoxon rank-sum test | retained and standardized |

### Assumption and prerequisite ledger

Explicit assumptions in the original draft:

- iid sampling
- positivity of observation probabilities
- finite-dimensional regression structure
- normality for the parametric continuous model
- no covariates in the simplified semi-parametric development

Tacit prerequisites identified during audit:

- correct model specification for the Bernoulli component whenever Wilks-style calibration is invoked
- correct parametric specification of the observed-outcome model in the parametric proposition
- interior parameter point, identifiability, and nonsingular Fisher information for the parametric likelihood-ratio theorem
- finite positive variance and convex-hull/interiority conditions for the empirical-likelihood mean constraint
- asymptotically nonvanishing treatment-group sizes
- stable interpretation of `\mu_\delta` and `\beta_\delta` across Method, Software, Simulation, and Application

### Structural risk points only

- The manuscript had a very flat structure, so transitions between sections carried most of the exposition load.
- The original draft overclaimed in the title, abstract, and discussion.
- The profile-likelihood notation in the main Method section was mathematically inconsistent with the intended procedure.
- The asymptotic propositions relied on strong unstated regularity conditions.
- The appendix derivation did not match the notation of the main text and profiled in the wrong direction.
- The software section and application depended on corrected interpretation of `\mu_\delta` and `\beta_\delta`, so any notation drift in Method propagated immediately downstream.

## Stage 4 — Section-Level Consistency and Exposition Audit

### Introduction

Mandatory patches:

- Removed unsupported placeholders about “current statistical practice” and “gaining popularity.”
- Replaced the claim that the method preserves “full statistical power” with a narrower, defensible claim about substantial power retention.
- Reframed the Wilcoxon comparison to avoid an over-precise statement about its null while preserving the practical point that it targets a different question.
- Tightened the final roadmap paragraph and added the missing reference to Section `\ref{sec:method}`.

Optional improvements implemented:

- Smoothed prose and reduced repetition around quality-of-life motivation.
- Standardized “Wilcoxon test” terminology.

### Method

Mandatory patches:

- Rewrote the scientific null hypothesis explicitly in terms of the observed binary and continuous components.
- Standardized expectation notation with `\E`.
- Kept `\Delta(X_i)` only as a cautionary contrast and explicitly separated it from the joint null.
- Corrected the profile likelihood-ratio construction so that profiling corresponds to maximizing the restricted likelihood rather than taking a supremum of the deviance.
- Strengthened both proposition statements with the regularity conditions they require.
- Added short proof sketches after both propositions.
- Repaired the semi-parametric notation, including the definition of `\widetilde{W}_{1,n}` and its profiled version.
- Clarified how parametric and semi-parametric confidence regions differ.

Optional improvements implemented:

- Replaced a few heavy sentences by shorter paragraphs.
- Added short lead-in phrases “Parametric version” and “Semi-parametric version” without creating new numbered subsections.

### Software implementation

Mandatory patches:

- Corrected the interpretation of the two reported contrasts to `\mu_\delta` and `\exp(\beta_\delta)`.
- Rephrased the example mean contrast as a sample difference instead of treating the displayed quantity as the theoretical `\Delta(X)`.
- Aligned the text with the corrected Method-section notation and confidence-region definitions.

Optional improvements implemented:

- Improved figure captions to say what the zero atom and simultaneous region represent.
- Simplified the explanation of package output.

### Simulation study

Mandatory patches:

- Rephrased the scenario descriptions to reflect the actual data-generating mechanisms more accurately.
- Replaced the “heavy tail” description of Scenario 4 with “strongly skewed,” which matches the beta model used.
- Softened the strongest power claims and tied them to the estimated power curves rather than absolute declarations.

Optional improvements implemented:

- Tightened grammar and improved transitions between scenario description, power results, and type I error results.

### Application

Mandatory patches:

- Corrected the typo in the description of the trial outcome.
- Reframed the zero outcome as a point mass at zero rather than as a stronger mechanistic claim than the data warrant.
- Softened the interpretation of the right-panel confidence region from a certainty statement to a data-supported suggestion.

Optional improvements implemented:

- Simplified the trial description and aligned it more closely with the manuscript’s notation.
- Improved the figure caption to name the parameter pair `( \mu_\delta, \beta_\delta )`.

### Discussion

Mandatory patches:

- Removed unresolved TODOs.
- Replaced promotional language with a balanced summary.
- Added the missing practical point about use as a prespecified primary analysis in an RCT.
- Added a limitations sentence about asymptotic calibration and the need for enough observed outcomes per group.

Optional improvements implemented:

- Shifted the conclusion from package-centric wording toward the statistical contribution.

### Appendix derivation

Mandatory patches:

- Rewrote the derivation so that the constrained empirical likelihood is defined directly in the notation used in the main text.
- Standardized the nuisance mean and mean-difference notation to `\mu` and `\mu_\delta`.
- Corrected the profile operation from `\sup` to `\inf`.
- Added the missing convex-hull/positive-weight condition for the empirical-likelihood constraints.

Optional improvements implemented:

- Removed unnecessary abstraction around generic distribution-function families and focused the appendix on the actual constrained optimization used later.

## Stage 5 — Proof Audit and Repair

### 1. Parametric asymptotic proposition

Independent reconstruction:

- Under the parametric model, the continuous component contributes a one-parameter likelihood-ratio restriction `\mu_\delta = 0`.
- The Bernoulli component contributes a one-parameter likelihood-ratio restriction `\beta_\delta = 0`.
- Because the likelihood factorizes and the parameter blocks are distinct, the joint restricted model is a regular two-restriction likelihood-ratio problem.

Verdict:

- `correct with strengthened assumptions`

Issues found:

- The original statement cited Wilks but omitted the regularity conditions needed for Wilks’ theorem.
- The surrounding text defined the profile statistic incorrectly, which made the proposition depend on a mis-stated object.
- The proof was absent.

Revised statement:

- Added correct-specification, interior-point, identifiability, and nonsingular-information assumptions.
- Kept positivity and treatment-group growth assumptions.
- Stated the null as `\mu_\delta = \beta_\delta = 0` under the working model.

Revised proof / repair sketch:

- Gave a short proof based on factorization of the joint likelihood and direct application of Wilks’ theorem to a two-restriction regular likelihood.

Minimal surrounding-text edits:

- Replaced the original profile-deviance formulas by true profiled likelihood-ratio definitions.
- Updated the software paragraph so the implementation matches the corrected statistic.

Downstream results affected:

- joint p-value computation
- critical-value explanation
- confidence-region inversion
- software description

### 2. Semi-parametric asymptotic proposition

Independent reconstruction:

- The Bernoulli component remains a regular parametric likelihood-ratio problem with one restriction.
- The observed-outcome component becomes a two-sample empirical-likelihood problem for a mean-difference constraint, profiled over the nuisance mean.
- Under standard empirical-likelihood regularity conditions, the observed-outcome component has a Wilks-type `\chi^2_1` limit.
- The total statistic adds one degree of freedom from each component.

Verdict:

- `correct with strengthened assumptions`

Issues found:

- The original statement was too weak on assumptions for empirical likelihood.
- The original text alternated between “non-parametric” and “semi-parametric.”
- The proof was absent.
- The surrounding notation for the profiled empirical-likelihood statistic was inconsistent.

Revised statement:

- Added no-covariate simplification explicitly.
- Added finite-variance and convex-hull/interiority conditions for the common mean.
- Kept positivity and asymptotic balance assumptions.

Revised proof / repair sketch:

- Added a short proof sketch invoking the Bernoulli Wilks result plus standard empirical-likelihood Wilks behavior for the mean-difference constraint, with the factorized structure supplying the joint statistic.

Minimal surrounding-text edits:

- Renamed the result “semi-parametric.”
- Defined `\widetilde{W}_{1,n}` and `\widetilde{W}_{1,n}^P` consistently before stating the proposition.

Downstream results affected:

- semi-parametric p-value calibration
- semi-parametric confidence regions
- interpretation of the package example and application

### 3. Appendix derivation of the empirical-likelihood component

Independent reconstruction:

- Start from the normalized empirical likelihood under mean constraints in each treatment group.
- Solve the constrained optimization problem with Lagrange multipliers.
- Obtain weights of the form `1 / (1 + \lambda_k(\cdot))`.
- Substitute back into the normalized likelihood ratio to obtain the deviance.
- Profile over the nuisance mean by minimizing the deviance.

Verdict:

- `proof incomplete but repairable`

Issues found:

- Main-text and appendix notation did not match.
- The original appendix defined the profiled statistic with `\sup` instead of `\inf`.
- The derivation did not say when the constrained weights are admissible.
- The role of the normalization by the unconstrained empirical likelihood was under-explained.

Revised statement:

- Recast the derivation directly in terms of `(\mu,\mu_\delta)` and `R_{1,n}^{\mathrm E}(\mu,\mu_\delta)`.

Revised proof / repair sketch:

- Rewrote the derivation from the constrained optimization problem through the Lagrange equations to the deviance formula and the profiled statistic.
- Added the convex-hull/positive-weight condition.

Minimal surrounding-text edits:

- Updated the main Method section to refer to the corrected appendix derivation.

Downstream results affected:

- the displayed semi-parametric deviance formula
- the profile definition of `\widetilde{W}_{1,n}^P`
- the semi-parametric proposition
- confidence-interval inversion for the semi-parametric procedure

## Stage 6 — Implemented Repairs

### Mandatory mathematical, notational, structural, and pedagogical repairs

- Created `manuscript/manuscript_new.tex` as a clean revised manuscript.
- Removed all inline TODO placeholders from the manuscript.
- Replaced promotional or unsupported claims in the title, abstract, introduction, and discussion.
- Rebuilt the Method section around an explicit observed-data null hypothesis.
- Corrected the profile likelihood-ratio definitions and made them match the intended implementation.
- Strengthened both asymptotic propositions and added proof sketches.
- Repaired the semi-parametric notation and appendix derivation.
- Standardized interpretation of `\mu_\delta` and `\beta_\delta` across the manuscript.
- Clarified the distinction between parametric and semi-parametric confidence regions and marginal intervals.
- Added discussion text about prespecified primary analyses in RCTs and about asymptotic limitations.

### Conservative optional style improvements

- Tightened prose throughout for journal tone and precision.
- Improved several figure captions to make the mathematical role of the plots explicit.
- Smoothed transitions between Introduction, Method, Software, Simulation, Application, and Discussion.
- Replaced awkward or imprecise phrases such as “heavy tail” and “extremely low p-value.”

### Validation

- Built `manuscript_new.tex` with `latexmk -pdf -interaction=nonstopmode -halt-on-error -outdir=/tmp/trunccomp-manuscript-new`.
- Result:
  - no LaTeX errors
  - no undefined references
  - no undefined citations
  - no remaining TODO placeholders
- Remaining build warnings:
  - supplementary tables still trigger `!h` to `!ht` float-specifier warnings from the included generated table files
  - one negligible overfull box warning remains in the opening sentence of the Method section

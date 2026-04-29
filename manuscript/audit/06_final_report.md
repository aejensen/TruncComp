# Executive summary

Phase 6 created `main_final.tex` from `manuscrip.tex` and implemented the mandatory repairs identified in the audit files. After the user explicitly expanded scope to include the previously skipped simulation material, the main simulation section and supplementary simulation appendix were also reviewed against the generated tables and saved simulation result object. The revision strengthens the regularity assumptions and proofs for all three proof-bearing propositions, recalibrates Delta-scale likelihood interval language, sharpens the separation between frequentist LRT calibration and Bayesian posterior modeling, and normalizes broad "continuous component" language to "non-atom component" where continuity is not required.

The simulation review confirmed that the saved result object is complete (`168/168` cells), all method failure rates are zero, and the quoted rejection rates in the simulation prose match the supplementary tables. Two local wording repairs were made in the simulation material: Scenario 3 now says the atom and non-atom components move in opposite directions, and the supplementary appendix now states the cell-ordering convention for random-number seeds more precisely. No generated simulation tables, generated simulation figures, or simulation result files were rerun or refreshed. The IST-3 Bayesian cache, generated application tables, generated application figures, and diagnostics were not rerun or refreshed.

No theorem labels, equation labels, citations, bibliography entries, theorem numbering, or existing cross-reference targets had to change. No new labels or citations were introduced.

# Result-by-result proof resolutions

Parametric likelihood-ratio proposition: strengthened the proposition statement to require a correctly specified finite-dimensional observed-data Gaussian/logistic model, fixed finite-dimensional covariate transformations, interior true parameter with `\sigma_0^2>0`, identifiability under unrestricted and null parameter spaces, Wilks regularity, existence and consistency of restricted and unrestricted MLEs, nonsingular observed-data Fisher information, and asymptotic information in both treatment arms and both atom/non-atom components. Replaced the overloaded treatment-arm ratio limit `c` by `\lambda_R`. Added a product-profile bridge before the proposition and rewrote the proof to apply Wilks to the full observed-data likelihood after establishing that the displayed component sum is the full profile likelihood-ratio deviance.

Covariate-adjusted semi-parametric LRT proposition `prop:covariate-adjusted-splrt`: added the conditional mean-zero implication of Equation `eq:glm2` after Equation `eq:adjusted-el-moment`. Strengthened the proposition with explicit interpretation of `Y_i` only on `{A_i=1}`, interior and identifiable finite-dimensional parameters, empirical-likelihood moment and rank conditions under the conditional law given `A_i=1`, convex-hull feasibility, and local nuisance uniqueness with probability tending to one. Rewrote the proof to separate the Bernoulli Wilks argument, the empirical-likelihood Wilks argument, and the joint CLT/block-diagonal Gaussian central-sequence argument needed for the `\chi^2_2` sum.

Delta-scale likelihood interval proposition `prop:delta-lr-intervals`: replaced the insufficient uniform `o_p(1)` convergence assumption for `g_n^{(e_0)}` with a root-`n` local expansion of the profiled coded-scale map and a local quadratic expansion of `\mathcal W_n`. Added positive definite efficient information, nonzero derivative of the population map, and nonempty/approximately minimized profile fibers. Rewrote the proof as a local-coordinate argument reducing the scalar profile statistic to the squared efficient distance to a hyperplane, yielding the `\chi^2_1` limit. Qualified projected Delta coverage under the same root-`n` map approximation.

# Assumption and notation changes

Added the population-map notation `g^{(e_0)}` before the Delta proposition and defined it as the population analogue of the profiled, standardized map.

Normalized broad terminology from "continuous component" to "non-atom component" in the abstract, Introduction roadmap, Method setup, semi-parametric subsection, Bayesian transition, and Discussion where the statement does not require true continuous support. Gaussian or genuinely continuous-support contexts retain continuous wording.

After scope expansion, normalized the simulation Scenario 3 description from "atom and continuous components" to "atom and non-atom components" because the simulated non-atom component is the relevant object.

Changed the parametric proposition's treatment-arm ratio limit from `c` to `\lambda_R` to avoid conflict with the coding map notation.

Clarified that the semi-parametric procedure is semi-parametric for the non-atom distributional shape but still depends on the Bernoulli model, the conditional mean model, empirical-likelihood feasibility, and large-sample regularity.

# Global claim recalibrations

Added "Under standard large-sample regularity conditions" to the abstract's one-p-value claim so that LRT calibration is not stated as assumption-free.

Added an application qualifier that the IST-3 LRT p-values are interpreted under the asymptotic calibrations described in the Method section.

Clarified in the Bayesian Method transition that the Bayesian model uses the same atom/non-atom observed-data decomposition for posterior summaries and checks, but does not inherit or replace frequentist likelihood-ratio calibration.

Added a Discussion limitation that likelihood-based coded-scale intervals require local smoothness and root-`n` stability of the profiled coded-scale map.

Changed the final software/practice sentence to "observed-data joint p-value" to preserve the observed-data target distinction.

Checked the main simulation claims against `simulation-study-results.rds` and the supplementary tables. The numerical claims for Scenarios 1, 2, 3, 5, and 6 match the generated results, and the stated zero failure-rate claim matches the saved result object.

# Expository and style changes

Added a concise product-profile bridge after Equation `eq:LRTtotal`.

Added a concise bridge from Equation `eq:glm2` to the empirical-likelihood estimating equation.

Changed Delta terminology from "profile likelihood-ratio interval" to "profile likelihood-ratio confidence set" and described `[L_{\mathrm{prof}}, U_{\mathrm{prof}}]` as an endpoint range, with interval interpretation when the accepted set is connected.

Updated the Delta numerical and appendix prose so the profile-grid construction is described as an endpoint range for a profile confidence set rather than assuming connectedness.

Clarified the supplementary simulation seed-order sentence so that cells are ordered by scenario, then effect level, then sample size, with sample size varying fastest within each scenario--effect block.

Kept style edits local and sentence-level. No generated numerical results, tables, figures, application summaries, simulation tables, or simulation figures were changed.

# Build and validation notes

Exact build command used:

```sh
latexmk -pdf -interaction=nonstopmode -halt-on-error main_final.tex
```

The repository Makefile was not used because its `pdf` target runs `assets`, which could refresh generated analysis artifacts; the locked IST-3 Bayesian fit/cache was not rerun.

`main_final.tex` built successfully to `main_final.pdf`. A final log check found no undefined citation warnings or undefined reference warnings.

The simulation section and supplementary simulation appendix were no longer treated as excluded after the user's scope expansion. Existing simulation artifacts were read but not regenerated; validation used the generated LaTeX tables and `simulation-study-results.rds`.

No theorem labels or references had to change, and no broken syntax, numbering, labels, cross-references, or citations introduced by the edits were detected in the final build.

# Remaining discretionary items

The optional title change to avoid "truncated continuous outcomes" was not made.

The optional addition of a label to the originally unlabeled parametric proposition was not made, to keep theorem labels stable.

The optional renaming of appendix grid-size integers `A,B` was not made because it was discretionary and would have touched more notation than required for the mandatory repairs.

The Discussion simulation paragraph was left unchanged because the simulation claims were already cautious and scenario-bound.

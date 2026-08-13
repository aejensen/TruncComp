#' Internal coded scale construction helpers
#'
#' Developer documentation for the helpers that compute the coded combined
#' outcome contrast and the image of an unadjusted likelihood sublevel
#' set.
#'
#' @name trunccomp-delta-projection-helpers
#' @title Internal coded scale construction helpers
#' @description
#' Coded scale likelihood intervals use the full criterion in the two arm
#' means and two arm logits. The implementation evaluates the sublevel range
#' by a method specific reduction: a two dimensional radial calculation
#' for the parametric likelihood and ordinary coded outcome empirical
#' likelihood for the semiparametric likelihood.
#' @aliases delta_na_interval
#' @aliases delta_from_components
#' @aliases delta_welch_interval
#' @aliases delta_unadjusted_point_estimate
#' @aliases prepareParametricJointReference
#' @aliases parametricJointCandidate
#' @aliases prepareSPLRTJointReference
#' @aliases splrtContinuousCandidate
#' @aliases splrtJointCandidate
#' @aliases delta_full_reference
#' @aliases delta_full_candidate
#' @aliases delta_bernoulli_arm_deviance
#' @aliases delta_logit_sublevel_bounds
#' @aliases delta_parametric_eta_deviance
#' @aliases delta_parametric_eta_candidate
#' @aliases delta_parametric_radial_limit
#' @aliases delta_parametric_radial_extreme
#' @aliases delta_parametric_sublevel_endpoint
#' @aliases delta_parametric_sublevel_interval
#' @aliases delta_splrt_sublevel_interval
#' @aliases delta_sublevel_range
#' @aliases delta_projected_interval.optimize
#' @details
#' `delta_from_components()` computes
#' `p1 * (mu1 - atom) - p0 * (mu0 - atom)`. The Welch helper constructs the
#' coded observations directly and uses their arm specific sample means and
#' variances.
#'
#' `delta_full_reference()` stores the arm samples, unrestricted arm means,
#' unrestricted arm logits, pooled residual sum of squares, and likelihood
#' constants for a successful unadjusted fit. `delta_full_candidate()`
#' evaluates the full criterion at `(mu0, mu1, eta0, eta1)` and is also used to
#' verify the reduced endpoint calculations.
#'
#' For the parametric likelihood, fix `eta = (eta0, eta1)` and let `W_A` be the
#' sum of its two Bernoulli deviances. If the sublevel cutoff is `c`, the
#' remaining mean deviance is `q = c - W_A`. Write `M` for the total number of
#' nonatom observations and `SSE` for the unrestricted pooled residual sum of
#' squares. The corresponding weighted mean ellipsoid has squared radius
#' `Q = SSE * expm1(q / M)`. If `B` is the coded contrast at the fitted arm
#' means and
#' `H = sqrt(p0^2 / m0 + p1^2 / m1)`, its exact extrema over that ellipsoid are
#' `B - sqrt(Q) * H` and `B + sqrt(Q) * H`.
#'
#' The remaining search is only over the compact Bernoulli sublevel set.
#' `delta_parametric_radial_limit()` writes
#' `eta = eta_hat + r * (cos(phi), sin(phi))`. Strict convexity gives one
#' boundary radius on every ray. The radial helpers use a dense deterministic
#' angle and radius search, refine the best numerical extremum, reconstruct the attaining
#' means, and require `delta_full_candidate()` to equal the requested cutoff.
#'
#' For the semiparametric likelihood, replace every atom by its numerical code
#' and retain each observed nonatom value. In one arm with `d` atoms, `m`
#' nonatoms, atom mass `rho`, nonatom mass `p`, and conditional nonatom weights
#' `q_j`, the Bernoulli likelihood ratio times the conditional empirical
#' likelihood ratio is exactly the ordinary empirical likelihood ratio for
#' the coded observations after setting each atom weight to `rho / d` and each
#' nonatom weight to `p * q_j`. Profiling subject to a coded mean therefore
#' gives the ordinary two sample coded outcome empirical likelihood exactly in
#' finite samples. `delta_splrt_sublevel_interval()` consequently inverts that
#' statistic at the supplied cutoff without numerical optimization over the
#' four component parameters.
#'
#' `delta_sublevel_range()` dispatches to these two reductions. A projected
#' interval supplies a `chisq(df = 4)` cutoff.
#' @seealso [confint.trunc_comp_fit()], [joint_contrast_surface()]
#' @keywords internal
NULL

#' Internal coded scale profile helpers
#'
#' Developer documentation for the scalar coded scale profile interval and
#' its wrappers.
#'
#' @name trunccomp-delta-profile-helpers
#' @title Internal coded scale profile helpers
#' @description
#' The accepted values from profiling the scalar coded contrast at cutoff `c`
#' equal the image of the full likelihood sublevel set at the same cutoff.
#' Profile and projected intervals therefore share the sublevel range
#' engine and differ only in their likelihood ratio cutoff.
#' @aliases delta_profile_interval.optimize
#' @aliases delta_projected_interval
#' @aliases delta_profile_interval
#' @aliases augmentDeltaInference
#' @details
#' If `T_D(d)` is the full criterion profiled subject to coded contrast `d`,
#' then
#' `set(d: T_D(d) <= c) = set(Delta(zeta): W(zeta) <= c)`.
#' `delta_profile_interval.optimize()` evaluates this range with the same
#' method specific reductions used for projection, but with the
#' `chisq(df = 1)` cutoff. It does not construct a separate target grid or run
#' an optimizer over a coded equality constraint.
#'
#' For `splrt`, this profile statistic and interval are exactly the ordinary
#' two sample empirical likelihood statistic and interval for the numerically
#' coded observations. For `lrt`, the analytic mean ellipsoid reduction leaves
#' only the deterministic radial search over the two arm logits.
#'
#' `augmentDeltaInference()` adds only the point estimate. Confidence
#' intervals are computed on demand by [confint.trunc_comp_fit()].
#' @seealso [confint.trunc_comp_fit()], [delta_projected_interval()]
#' @keywords internal
NULL

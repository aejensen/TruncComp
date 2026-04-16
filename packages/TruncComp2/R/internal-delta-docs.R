#' Internal Delta construction and projection helpers
#'
#' Developer-facing documentation for the helpers that convert joint component
#' effects into the derived combined-outcome contrast `Delta` and project the
#' joint confidence region onto that scale.
#'
#' @name trunccomp2-delta-projection-helpers
#' @title Internal Delta construction and projection helpers
#' @description
#' These helpers compute `Delta`, prepare cached candidate objects for the
#' parametric and semi-parametric methods, and implement the optimizer-based
#' projection interval.
#' @aliases delta_na_interval
#' @aliases delta_from_components
#' @aliases delta_welch_interval
#' @aliases delta_unadjusted_point_estimate
#' @aliases prepareParametricJointReference
#' @aliases parametricJointCandidate
#' @aliases prepareSPLRTJointReference
#' @aliases splrtContinuousCandidate
#' @aliases splrtJointCandidate
#' @aliases delta_projection_candidate_factory
#' @aliases delta_projection_penalty_scale
#' @aliases delta_projection_center
#' @aliases delta_projection_start_points
#' @aliases delta_projection_objective
#' @aliases delta_projection_is_feasible
#' @aliases delta_projection_boundary_candidate
#' @aliases delta_projection_better
#' @aliases delta_projection_in_box
#' @aliases delta_projection_near_boundary
#' @aliases delta_projection_expand_bounds
#' @aliases delta_projected_interval.optimize
#' @aliases parametricContinuousCandidate
#' @usage
#' delta_na_interval()
#' delta_from_components(atom, mu0, mu1, pi0, pi1)
#' delta_welch_interval(y, r, conf.level)
#' delta_unadjusted_point_estimate(object, tol = 1e-8)
#' prepareParametricJointReference(data, atom)
#' parametricJointCandidate(parametricReference, muDelta, logORdelta, tol = 1e-12)
#' prepareSPLRTJointReference(data, atom)
#' splrtContinuousCandidate(splrtReference, muDelta, p0, p1, tol = 1e-8)
#' splrtJointCandidate(splrtReference, muDelta, logORdelta, tol = 1e-8)
#' delta_projection_candidate_factory(object)
#' delta_projection_penalty_scale(object)
#' delta_projection_center(object, bounds)
#' delta_projection_start_points(object, bounds)
#' delta_projection_objective(par, candidate_fun, threshold, direction, penalty_scale)
#' delta_projection_is_feasible(candidate, threshold, tol = 1e-8)
#' delta_projection_boundary_candidate(candidate_fun, from, to, threshold, tol = 1e-8, max_iter = 60)
#' delta_projection_better(candidate, best, direction, tol = 1e-8)
#' delta_projection_in_box(object, candidate_fun, bounds, threshold, direction = c("lower", "upper"), tol = 1e-8)
#' delta_projection_near_boundary(solution, bounds, frac = 0.05, tol = 1e-6)
#' delta_projection_expand_bounds(bounds, center, expansion = 1.5, additive = c(0.25, 0.25))
#' delta_projected_interval.optimize(object, conf.level = object$conf.level, offset = NULL, expansion = 1.5, max_iter = 4, tol = 1e-8)
#' parametricContinuousCandidate(parametricReference, muDelta, tol = 1e-12)
#' @details
#' ### `delta_na_interval()`
#'
#' Returns the package-wide missing-value interval placeholder `c(NA, NA)`. Its
#' role is purely structural: all unsupported or failed `Delta` interval slots
#' use the same sentinel shape.
#'
#' ### `delta_from_components(atom, mu0, mu1, pi0, pi1)`
#'
#' Computes the derived combined-outcome mean difference from arm-specific
#' observed-outcome means (`mu0`, `mu1`), observation probabilities (`pi0`,
#' `pi1`), and the atom value. The helper returns a scalar numeric `Delta`. Its
#' role is to encode the defining transformation used everywhere in the package.
#'
#' ### `delta_welch_interval(y, r, conf.level)`
#'
#' Builds the descriptive raw-scale `DeltaMarginalCI` by treating the combined
#' outcome `y` as a standard two-sample mean problem split by treatment indicator
#' `r`. It returns a length-two numeric interval and falls back to a point-mass
#' interval when the pooled standard error is not usable. Its role is to provide
#' the model-light `Delta` interval reported for successful unadjusted fits.
#'
#' ### `delta_unadjusted_point_estimate(object, tol = 1e-8)`
#'
#' Reconstructs the fitted `Delta` point estimate for an unadjusted successful
#' object. For the semi-parametric method it first solves the empirical-
#' likelihood nuisance location so the arm-specific observed means are coherent
#' with the fitted mean difference; otherwise it uses the raw observed means. It
#' returns a scalar numeric value. Its role is to populate `object$Delta` in a
#' way that matches the fitted method.
#'
#' ### `prepareParametricJointReference(data, atom)` and
#' `prepareSPLRTJointReference(data, atom)`
#'
#' Build cached reference objects for the parametric and semi-parametric joint
#' candidate evaluators. `data` is the standardized analysis frame and `atom` is
#' the fitted atom value. Each helper returns a list containing the cached model
#' fits or empirical-likelihood ingredients required for repeated candidate
#' evaluation. Their role is to avoid refitting the expensive shared pieces at
#' every surface or optimization step.
#'
#' ### `parametricJointCandidate(parametricReference, muDelta, logORdelta,
#' tol = 1e-12)`
#'
#' Evaluates a constrained parametric candidate at fixed component effects. The
#' helper fits offset logistic and linear models, returns the Bernoulli and
#' Normal likelihood-ratio contributions, arm-specific parameters, the implied
#' `Delta`, and the constrained models themselves. Infeasible constrained fits
#' yield `Inf` components or `NA` parameters. Its role is to provide the core
#' parametric candidate object used by simultaneous-region and `Delta`
#' calculations.
#'
#' ### `splrtContinuousCandidate(splrtReference, muDelta, p0, p1, tol = 1e-8)`
#'
#' Evaluates the continuous empirical-likelihood component of the semi-
#' parametric candidate given a proposed observed-outcome mean difference
#' `muDelta` and arm-specific observation probabilities `p0` and `p1`. The
#' helper returns feasibility, the continuous likelihood-ratio statistic, the
#' implied arm-specific means, and the resulting `Delta`. Its role is to isolate
#' the continuous half of the semi-parametric joint candidate.
#'
#' ### `splrtJointCandidate(splrtReference, muDelta, logORdelta, tol = 1e-8)`
#'
#' Combines the profiled logistic candidate with the empirical-likelihood
#' continuous candidate for a semi-parametric joint evaluation. The helper
#' returns the total statistic, component contributions, arm-specific
#' probabilities and means, and the implied `Delta`. It truncates very small
#' negative totals to zero. Its role is the semi-parametric analogue of
#' `parametricJointCandidate()`.
#'
#' ### `delta_projection_candidate_factory(object)`
#'
#' Chooses the appropriate candidate generator for the fitted method stored in
#' `object`. It returns a closure taking `(muDelta, logORdelta)` and producing a
#' joint candidate object. It errors for adjusted fits and unsupported methods.
#' Its role is to hide the parametric versus semi-parametric branching from the
#' projection optimizer.
#'
#' ### `delta_projection_penalty_scale(object)`
#'
#' Derives a scale factor used to penalize infeasible candidates during the
#' constrained projection optimization. The helper returns a positive scalar
#' based on the observed range of `Y`, defaulting safely to `1`. Its role is to
#' make the optimizer strongly prefer feasible points without being completely
#' detached from the scale of the outcome.
#'
#' ### `delta_projection_center(object, bounds)` and
#' `delta_projection_start_points(object, bounds)`
#'
#' `delta_projection_center()` returns the preferred central point in the search
#' box, using the fitted component estimates when finite and falling back to the
#' midpoint of `bounds`. `delta_projection_start_points()` expands that center
#' into a deduplicated set of corner and midpoint starting values. Their role is
#' to seed the optimizer from several plausible locations so the projection
#' search is less sensitive to local numerical issues.
#'
#' ### `delta_projection_objective(par, candidate_fun, threshold, direction,
#' penalty_scale)`
#'
#' Objective used by `nlminb()` for the optimizer-based projected interval.
#' `par` holds `(muDelta, logORdelta)`, `candidate_fun` builds the joint
#' candidate, `threshold` is the chi-squared cutoff, `direction` chooses whether
#' the optimizer is minimizing or maximizing `Delta`, and `penalty_scale`
#' controls the quadratic penalty for violating the joint constraint. The helper
#' returns a scalar numeric objective value. Its role is to convert the
#' constrained projection problem into an unconstrained box optimization.
#'
#' ### `delta_projection_is_feasible(candidate, threshold, tol = 1e-8)`
#'
#' Predicate returning `TRUE` when a joint candidate has finite `Delta`, finite
#' statistic, and satisfies the joint-threshold constraint up to tolerance. Its
#' role is to keep the boundary search and result selection logic readable.
#'
#' ### `delta_projection_boundary_candidate(candidate_fun, from, to, threshold,
#' tol = 1e-8, max_iter = 60)`
#'
#' Walks from a known feasible point `from` toward `to` until it reaches the
#' edge of the acceptance region. The helper returns both the boundary
#' parameter vector and the corresponding candidate object. It uses binary search
#' and gracefully returns infeasible endpoints when no feasible boundary point
#' exists. Its role is to translate unconstrained optimizer proposals into
#' feasible projected-interval candidates.
#'
#' ### `delta_projection_better(candidate, best, direction, tol = 1e-8)`
#'
#' Compares two feasible projection candidates. The helper returns `TRUE` when
#' `candidate` improves the lower or upper `Delta` endpoint relative to `best`,
#' breaking ties in favor of the smaller test statistic. Its role is to keep the
#' interval-endpoint selection logic centralized.
#'
#' ### `delta_projection_in_box(object, candidate_fun, bounds, threshold,
#' direction = c("lower", "upper"), tol = 1e-8)`
#'
#' Searches a bounded box in `(muDelta, logORdelta)` space for the best feasible
#' lower or upper `Delta` value. The helper returns the best boundary candidate
#' found or `NULL` if none is feasible. It evaluates multiple start points and
#' refines them with `nlminb()`. Its role is the core workhorse for the
#' optimizer-based projected interval.
#'
#' ### `delta_projection_near_boundary(solution, bounds, frac = 0.05,
#' tol = 1e-6)` and `delta_projection_expand_bounds(bounds, center,
#' expansion = 1.5, additive = c(0.25, 0.25))`
#'
#' The first helper checks whether a proposed solution lies suspiciously close to
#' the current box boundary; the second expands the search box around `center`.
#' They return a logical flag and an expanded bounds list, respectively. Their
#' role is to drive the iterative box-expansion strategy used by the optimizer.
#'
#' ### `delta_projected_interval.optimize(object, conf.level = object$conf.level,
#' offset = NULL, expansion = 1.5, max_iter = 4, tol = 1e-8)`
#'
#' Computes the optimizer-based projected `Delta` interval. The helper returns a
#' length-two numeric interval or `c(NA, NA)` when no stable feasible endpoints
#' are found within the repeated box expansions. Its role is to provide the
#' slower but direct alternative to the default surface-based projection.
#'
#' ### `parametricContinuousCandidate(parametricReference, muDelta,
#' tol = 1e-12)`
#'
#' Evaluates only the constrained Normal component for a fixed observed-outcome
#' mean difference `muDelta` under the parametric model. The helper returns the
#' likelihood-ratio contribution and the implied arm-specific observed means. Its
#' role is to support the one-dimensional `Delta` profile optimization, where
#' the Bernoulli part is handled separately.
#'
#' @seealso [delta_profile_interval()], [jointContrastCI()], [augmentDeltaInference()]
#' @keywords internal
NULL

#' Internal Delta profile and surface helpers
#'
#' Developer-facing documentation for the helpers that build grid-based and
#' optimizer-based `Delta` intervals and attach them to fitted objects.
#'
#' @name trunccomp2-delta-profile-helpers
#' @title Internal Delta profile and surface helpers
#' @description
#' These helpers choose scales and tolerances, profile over the log-odds axis,
#' derive intervals from joint surfaces, and augment successful unadjusted fits
#' with their stored `Delta` summaries.
#' @aliases delta_profile_value_scale
#' @aliases delta_profile_initial_step
#' @aliases delta_profile_target_tolerance
#' @aliases delta_profile_interval_width
#' @aliases delta_profile_cache_key
#' @aliases delta_profile_lrt_candidate_factory
#' @aliases delta_profile_splrt_mu_solver
#' @aliases delta_profile_splrt_candidate_factory
#' @aliases delta_profile_candidate_factory
#' @aliases delta_profile_objective_value
#' @aliases delta_profile_logor_near_boundary
#' @aliases delta_profile_optimize_logor
#' @aliases delta_profile_factory
#' @aliases delta_profile_find_bound
#' @aliases delta_profile_interval.optimize
#' @aliases delta_boundary_touched
#' @aliases delta_surface_for_inference
#' @aliases delta_interval_from_surface
#' @aliases validateDeltaIntervalAlgorithm
#' @aliases delta_projected_interval.grid
#' @aliases delta_projected_interval
#' @aliases delta_profile_interval.grid
#' @aliases delta_profile_interval
#' @aliases augmentDeltaInference
#' @usage
#' delta_profile_value_scale(object)
#' delta_profile_initial_step(object)
#' delta_profile_target_tolerance(target, scale = 1)
#' delta_profile_interval_width(object)
#' delta_profile_cache_key(x)
#' delta_profile_lrt_candidate_factory(object, tol = 1e-8)
#' delta_profile_splrt_mu_solver(reference, p0, p1, targetDelta, tol = 1e-8)
#' delta_profile_splrt_candidate_factory(object, tol = 1e-8)
#' delta_profile_candidate_factory(object, tol = 1e-8)
#' delta_profile_objective_value(candidate)
#' delta_profile_logor_near_boundary(logORdelta, lower, upper, frac = 0.05, tol = 1e-6)
#' delta_profile_optimize_logor(candidate_fun, center, width, expansion = 1.5, max_iter = 6, tol = 1e-6)
#' delta_profile_factory(object, tol = 1e-8)
#' delta_profile_find_bound(profile_fun, estimate, crit, direction = c("lower", "upper"), initial_step, tol = 1e-6, max_expand = 24)
#' delta_profile_interval.optimize(object, conf.level = object$conf.level, tol = 1e-6)
#' delta_boundary_touched(mask)
#' delta_surface_for_inference(object, conf.level = object$conf.level, resolution = 31, offset = NULL, expansion = 1.5, max_iter = 6)
#' delta_interval_from_surface(surface_data, threshold)
#' validateDeltaIntervalAlgorithm(algorithm)
#' delta_projected_interval.grid(object, conf.level = object$conf.level, offset = NULL, resolution = 35, expansion = 1.5, max_iter = 6)
#' delta_projected_interval(object, conf.level = object$conf.level, offset = NULL, resolution = 35, algorithm = c("grid", "optimize"), expansion = 1.5, max_iter = 6, tol = 1e-8)
#' delta_profile_interval.grid(object, conf.level = object$conf.level, resolution = 31, offset = NULL, expansion = 1.5, max_iter = 6)
#' delta_profile_interval(object, conf.level = object$conf.level, resolution = 31, offset = NULL, algorithm = c("grid", "optimize"), expansion = 1.5, max_iter = 6, tol = 1e-6)
#' augmentDeltaInference(object)
#' @details
#' ### `delta_profile_value_scale(object)`,
#' `delta_profile_initial_step(object)`, and
#' `delta_profile_interval_width(object)`
#'
#' These helpers derive problem-specific scales for the `Delta` profile search.
#' They return, respectively, a characteristic `Delta` scale, the initial
#' outward step used when bracketing profile bounds, and the starting half-width
#' for the log-odds search interval. Their role is to keep the profiling
#' algorithm adaptive to the observed data rather than hard-coding a universal
#' search range.
#'
#' ### `delta_profile_target_tolerance(target, scale = 1)` and
#' `delta_profile_cache_key(x)`
#'
#' `delta_profile_target_tolerance()` returns a numerical tolerance appropriate
#' for comparing candidate and target `Delta` values, while
#' `delta_profile_cache_key()` converts numeric inputs into a stable string key
#' for memoization environments. Their role is to control repeated candidate
#' evaluations without introducing avoidable floating-point mismatches.
#'
#' ### `delta_profile_lrt_candidate_factory(object, tol = 1e-8)`
#'
#' Builds the parametric candidate closure used by the direct `Delta` profile
#' optimizer. The returned function maps a target `Delta` and a candidate
#' `logORdelta` into a full candidate object with implied `muDelta`, joint
#' statistic, arm-specific parameters, and the reconstructed `Delta`. Its role
#' is to reduce the direct profile problem to optimization over a single
#' log-odds coordinate.
#'
#' ### `delta_profile_splrt_mu_solver(reference, p0, p1, targetDelta,
#' tol = 1e-8)`
#'
#' Solves the semi-parametric subproblem of finding the observed-outcome mean
#' difference `muDelta` that yields a requested `targetDelta` once `p0` and `p1`
#' are fixed. The helper returns a candidate list carrying feasibility, the
#' empirical-likelihood statistic, the arm-specific means, and the matched
#' `muDelta`. It uses root-finding when possible and falls back to minimizing the
#' squared `Delta` mismatch over the feasible range. Its role is to support the
#' direct profile interval for the semi-parametric method.
#'
#' ### `delta_profile_splrt_candidate_factory(object, tol = 1e-8)` and
#' `delta_profile_candidate_factory(object, tol = 1e-8)`
#'
#' `delta_profile_splrt_candidate_factory()` builds the semi-parametric version
#' of the target-`Delta` candidate closure. `delta_profile_candidate_factory()`
#' dispatches between the parametric and semi-parametric builders based on
#' `object$method`. Both return closures that accept a target `Delta` and a
#' candidate `logORdelta`. Their role is to give the later optimizer a uniform
#' interface across methods.
#'
#' ### `delta_profile_objective_value(candidate)` and
#' `delta_profile_logor_near_boundary(logORdelta, lower, upper, frac = 0.05,
#' tol = 1e-6)`
#'
#' Small predicates used by the direct profile optimizer. The first converts a
#' candidate object into a scalar objective, treating missing or infeasible
#' candidates as very poor. The second flags whether the best current log-odds
#' location lies too close to the search boundary. Their role is to keep the
#' optimizer bookkeeping readable.
#'
#' ### `delta_profile_optimize_logor(candidate_fun, center, width,
#' expansion = 1.5, max_iter = 6, tol = 1e-6)`
#'
#' One-dimensional search over `logORdelta` for a fixed target `Delta`. The
#' helper repeatedly expands a search interval around `center`, evaluates the
#' candidate closure and a direct `optimize()` call, and returns the best
#' candidate object it found. Its role is the core numerical engine for the
#' direct `Delta` profile interval.
#'
#' ### `delta_profile_factory(object, tol = 1e-8)`
#'
#' Creates the memoized profile function `targetDelta -> best candidate` for a
#' fitted object. The returned closure caches results by target `Delta`, returns
#' the fitted point with statistic zero when the target matches `object$Delta`,
#' and otherwise delegates to `delta_profile_optimize_logor()`. Its role is to
#' expose a reusable profile function for interval inversion.
#'
#' ### `delta_profile_find_bound(profile_fun, estimate, crit,
#' direction = c("lower", "upper"), initial_step, tol = 1e-6,
#' max_expand = 24)`
#'
#' Inverts a one-dimensional profile function to locate one endpoint of the
#' direct `Delta` profile interval. It returns a scalar bound or `NA_real_` if
#' the search cannot leave the acceptance region. Its role is to translate the
#' profiled statistic into interval endpoints after the candidate optimizer is in
#' place.
#'
#' ### `delta_profile_interval.optimize(object, conf.level = object$conf.level,
#' tol = 1e-6)`
#'
#' Computes the direct optimizer-based `Delta` profile interval. It returns a
#' length-two numeric interval or `c(NA, NA)` when the fit is unsuccessful,
#' adjusted, lacks an atom, or the profiling search cannot find both endpoints.
#' Its role is to provide the slower exact-profile alternative to the default
#' grid-based approximation.
#'
#' ### `delta_boundary_touched(mask)`, `delta_surface_for_inference(object,
#' conf.level = object$conf.level, resolution = 31, offset = NULL,
#' expansion = 1.5, max_iter = 6)`, and
#' `delta_interval_from_surface(surface_data, threshold)`
#'
#' These helpers support the default grid-based `Delta` intervals. `delta_boundary_touched()`
#' checks whether the accepted surface region hits the edge of the current grid.
#' `delta_surface_for_inference()` repeatedly evaluates the joint surface with
#' `include_delta = TRUE`, expanding the offsets until the accepted region is
#' interior or the iteration limit is reached. `delta_interval_from_surface()`
#' then takes the range of accepted `Delta` values under a supplied threshold.
#' Their role is to turn the same evaluated joint surface into either projected
#' or profile-style `Delta` intervals.
#'
#' ### `validateDeltaIntervalAlgorithm(algorithm)`
#'
#' Validates the algorithm selector used by the public `confint()` method and the
#' internal wrappers. It returns either `"grid"` or `"optimize"` and errors on
#' any other value. Its role is to keep the algorithm switch centralized.
#'
#' ### `delta_projected_interval.grid(...)`, `delta_projected_interval(...)`,
#' `delta_profile_interval.grid(...)`, and `delta_profile_interval(...)`
#'
#' These wrappers expose the grid-based and optimizer-based projected/profile
#' `Delta` intervals behind uniform APIs. The `*.grid` helpers always operate on
#' the expanded joint surface, while the higher-level wrappers dispatch on the
#' validated `algorithm` argument. Each returns a length-two numeric interval.
#' Their role is to keep the public confidence-interval code small while still
#' allowing slower exact alternatives for validation or sensitivity analysis.
#'
#' ### `augmentDeltaInference(object)`
#'
#' Final augmentation step applied to successful fits. For unadjusted fits with a
#' known atom it computes the stored `Delta` point estimate, the Welch
#' `DeltaMarginalCI`, leaves `DeltaProjectedCI` as unavailable by default, and
#' stores the grid-based `DeltaProfileCI` as both `DeltaProfileCI` and the
#' backward-compatible alias `DeltaCI`. For failed, adjusted, or atom-less fits
#' it fills all `Delta` fields with `NA` intervals. Its role is to attach the
#' derived combined-outcome summaries only when they are methodologically
#' meaningful.
#'
#' @seealso [confint.TruncComp2()], [jointContrastCI()], [delta_projected_interval()]
#' @keywords internal
NULL

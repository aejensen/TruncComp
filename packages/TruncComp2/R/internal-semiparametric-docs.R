#' Internal semi-parametric likelihood-ratio helpers
#'
#' Developer-facing documentation for the high-level semi-parametric fitting
#' wrappers used by `method = "SPLRT"`.
#'
#' @name trunccomp2-splrt-helpers
#' @title Internal semi-parametric likelihood-ratio helpers
#' @description
#' These helpers combine empirical-likelihood inference for the observed-outcome
#' treatment effect with logistic likelihood-ratio inference for the observation
#' process.
#' @aliases adjusted_SPLRT
#' @aliases SPLRT
#' @usage
#' adjusted_SPLRT(
#'   data, conf.level = 0.95, adjust = NULL, adjust_spec = NULL, atom = NULL
#' )
#' SPLRT(data, conf.level = 0.95, adjust = NULL, adjust_spec = NULL, atom = NULL)
#' @details
#' ### `adjusted_SPLRT(data, conf.level = 0.95, adjust = NULL,
#' adjust_spec = NULL, atom = NULL)`
#'
#' Fits the adjusted semi-parametric procedure. `data` is the standardized
#' analysis frame, `conf.level` controls the marginal intervals, `adjust` is the
#' additive baseline-covariate formula, `adjust_spec` is the printable metadata
#' string stored on the result, and `atom` is carried through for consistency.
#' The helper returns a successful or failed `"TruncComp2"` object. It rejects
#' irregular adjusted logistic fits using the parametric regularity checks and
#' delegates the observed-outcome component to `el_regression_fit()`. Its role
#' is to keep the adjusted `SPLRT` logic separate from the simpler unadjusted
#' empirical-likelihood path while preserving the same result structure.
#'
#' ### `SPLRT(data, conf.level = 0.95, adjust = NULL, adjust_spec = NULL,
#' atom = NULL)`
#'
#' Main semi-parametric estimator called by `truncComp_core()`. In the
#' unadjusted case it extracts the observed outcomes by arm, runs the
#' two-sample empirical-likelihood test with `el_mean_diff_fit()`, fits the
#' logistic observation models with `glm()`, and adds the two component
#' likelihood-ratio statistics. When `adjust` is non-`NULL`, it delegates to
#' `adjusted_SPLRT()`. The helper returns a `"TruncComp2"` object and, on
#' success, passes it through `augmentDeltaInference()` so unadjusted fits gain
#' the stored `Delta` summaries. Its role is to expose a single semi-parametric
#' implementation point to the rest of the package.
#'
#' @seealso [el_mean_diff_fit()], [el_regression_fit()], [augmentDeltaInference()]
#' @keywords internal
NULL

#' Internal empirical-likelihood helpers
#'
#' Developer-facing documentation for the pure-R empirical-likelihood engine
#' used by the semi-parametric path and its adjusted regression extension.
#'
#' @name trunccomp2-empirical-likelihood-helpers
#' @title Internal empirical-likelihood helpers
#' @description
#' These helpers cover both the unadjusted two-sample mean-difference problem
#' and the adjusted regression-coefficient profile problem solved by the
#' semi-parametric method.
#' @aliases el_validate_mean_diff_inputs
#' @aliases el_root_epsilon
#' @aliases el_one_sample_lambda
#' @aliases el_mean_diff_delta_range
#' @aliases el_mean_diff_theta_range
#' @aliases el_mean_diff_theta_eval
#' @aliases el_mean_diff_theta
#' @aliases el_mean_diff_statistic
#' @aliases el_find_mean_diff_bound
#' @aliases el_mean_diff_confint
#' @aliases el_mean_diff_fit
#' @aliases el_regression_design
#' @aliases el_regression_beta
#' @aliases el_regression_restricted_ls
#' @aliases el_regression_moments
#' @aliases el_dual_state
#' @aliases el_dual_fit_newton
#' @aliases el_dual_fit_nlminb
#' @aliases el_dual_fit
#' @aliases el_regression_profile_fit
#' @aliases el_regression_statistic
#' @aliases el_regression_variance
#' @aliases el_regression_ci_bound
#' @aliases el_regression_confint
#' @aliases el_regression_fit
#' @usage
#' el_validate_mean_diff_inputs(x, y, conf.level = NULL)
#' el_root_epsilon(lower, upper)
#' el_one_sample_lambda(residuals, tol = 1e-10)
#' el_mean_diff_delta_range(x, y)
#' el_mean_diff_theta_range(delta, x, y)
#' el_mean_diff_theta_eval(theta, delta, x, y)
#' el_mean_diff_theta(delta, x, y, tol = 1e-8)
#' el_mean_diff_statistic(x, y, delta, tol = 1e-8)
#' el_find_mean_diff_bound(x, y, estimate, boundary, crit, side, tol = 1e-8)
#' el_mean_diff_confint(x, y, estimate, conf.level, tol = 1e-8)
#' el_mean_diff_fit(x, y, mu = 0, conf.level = 0.95)
#' el_regression_design(formula, data, term = "R")
#' el_regression_beta(delta, nuisance, term_index, p)
#' el_regression_restricted_ls(design, delta)
#' el_regression_moments(design, beta)
#' el_dual_state(lambda, G)
#' el_dual_fit_newton(G, tol = 1e-8, maxit = 100)
#' el_dual_fit_nlminb(G, start = NULL, tol = 1e-8)
#' el_dual_fit(G, tol = 1e-8)
#' el_regression_profile_fit(design, delta, tol = 1e-8)
#' el_regression_statistic(design, delta, tol = 1e-8)
#' el_regression_variance(design)
#' el_regression_ci_bound(design, estimate, crit, direction, se, tol = 1e-8)
#' el_regression_confint(design, estimate, conf.level, tol = 1e-8)
#' el_regression_fit(
#'   data, formula, term = "R", mu = 0, conf.level = 0.95, tol = 1e-8
#' )
#' @details
#' ### `el_validate_mean_diff_inputs(x, y, conf.level = NULL)`
#'
#' Validates the two numeric samples supplied to the unadjusted empirical
#' likelihood routines. `x` and `y` must be finite numeric vectors with at least
#' two observations each, and `conf.level`, when present, must lie strictly
#' between `0` and `1`. The helper returns `NULL` invisibly and throws
#' descriptive errors on invalid input. Its role is to protect all later
#' unadjusted empirical-likelihood calculations from malformed data.
#'
#' ### `el_root_epsilon(lower, upper)`
#'
#' Computes a small boundary offset for root-finding on a finite interval. The
#' helper returns a scalar epsilon chosen from the interval span and machine
#' precision. Its role is to keep `uniroot()` and similar searches away from the
#' exact feasibility boundary where denominators can hit zero.
#'
#' ### `el_one_sample_lambda(residuals, tol = 1e-10)`
#'
#' Solves the one-sample empirical-likelihood Lagrange-multiplier problem for a
#' residual vector. The helper returns a list containing feasibility,
#' `lambda`, the one-sample likelihood-ratio statistic, and the profile
#' gradient. If the residuals are all near zero the solution is returned in
#' closed form; if the feasible interval does not bracket zero it returns
#' `feasible = FALSE` and `statistic = Inf`. Its role is to provide the basic
#' empirical-likelihood building block for both arms of the unadjusted
#' mean-difference problem.
#'
#' ### `el_mean_diff_delta_range(x, y)` and
#' `el_mean_diff_theta_range(delta, x, y)`
#'
#' Return the feasible ranges for the mean difference `delta` and the nuisance
#' location parameter `theta` in the two-sample problem. Each helper returns a
#' length-two numeric interval derived from sample extrema. Their role is to
#' bound the later optimization and root-finding steps.
#'
#' ### `el_mean_diff_theta_eval(theta, delta, x, y)`
#'
#' Evaluates the two-sample empirical-likelihood objective at a specific
#' nuisance value `theta` for fixed mean difference `delta`. The helper returns a
#' list with feasibility, the combined statistic, and the combined gradient. If
#' either arm is infeasible it reports `statistic = Inf`. Its role is to provide
#' the primitive objective queried by the `theta` optimizer.
#'
#' ### `el_mean_diff_theta(delta, x, y, tol = 1e-8)`
#'
#' Finds the best nuisance parameter `theta` for a fixed mean difference `delta`.
#' It returns the same list structure as `el_mean_diff_theta_eval()`, evaluated
#' at the minimizing `theta`. The helper first looks for a gradient root and
#' falls back to direct one-dimensional optimization when necessary. It reports
#' infeasibility when the derived `theta` range is empty. Its role is to profile
#' the nuisance parameter out of the unadjusted two-sample empirical-likelihood
#' statistic.
#'
#' ### `el_mean_diff_statistic(x, y, delta, tol = 1e-8)`
#'
#' Returns the profiled empirical-likelihood statistic for testing the mean
#' difference `delta`. The helper short-circuits to `0` at the sample mean
#' difference and to `Inf` outside the feasible range. Its role is to expose a
#' one-dimensional statistic that can be inverted into intervals or inserted into
#' the semi-parametric joint test.
#'
#' ### `el_find_mean_diff_bound(x, y, estimate, boundary, crit, side,
#' tol = 1e-8)`
#'
#' Finds one endpoint of the unadjusted empirical-likelihood confidence interval.
#' `estimate` is the empirical mean-difference estimate, `boundary` is the
#' feasible extreme on the requested side, `crit` is the chi-squared cutoff, and
#' `side` chooses `"lower"` or `"upper"`. The helper returns the interval
#' endpoint as a scalar numeric value. It uses nested binary searches and treats
#' non-finite objective values as outside the acceptance region. Its role is to
#' invert the statistic in a numerically stable way.
#'
#' ### `el_mean_diff_confint(x, y, estimate, conf.level, tol = 1e-8)`
#'
#' Builds the full unadjusted empirical-likelihood confidence interval by
#' calling `el_find_mean_diff_bound()` on both sides of the estimate. It returns
#' a length-two numeric vector with a `"conf.level"` attribute. Its role is to
#' provide the interval reported by `el_mean_diff_fit()` and the unadjusted
#' semi-parametric method.
#'
#' ### `el_mean_diff_fit(x, y, mu = 0, conf.level = 0.95)`
#'
#' High-level wrapper for the unadjusted two-sample empirical-likelihood test.
#' It returns an `"htest"` object with the mean-difference estimate, confidence
#' interval, test statistic, p-value, and descriptive labels. Its role is to
#' preserve a standard R-testing-object shape while providing the observed-outcome
#' component for the unadjusted `SPLRT` path.
#'
#' ### `el_regression_design(formula, data, term = "R")`
#'
#' Constructs the design object for the adjusted empirical-likelihood regression
#' extension. `formula` is typically the observed-outcome model, `data` is the
#' observed-outcome subset, and `term` names the treatment coefficient of
#' interest. The helper returns either `success = TRUE` with the response,
#' model matrix, QR decomposition, least-squares fit, and coefficient estimates,
#' or `success = FALSE` with an explanatory error. It rejects rank-deficient or
#' too-small designs. Its role is to centralize all adjusted-design validation.
#'
#' ### `el_regression_beta(delta, nuisance, term_index, p)`
#'
#' Reconstructs the full regression coefficient vector when the treatment effect
#' is fixed at `delta` and the nuisance coefficients are supplied separately.
#' The helper returns a numeric vector of length `p`. Its role is to keep the
#' profile and moment calculations agnostic to how the treatment coefficient is
#' singled out.
#'
#' ### `el_regression_restricted_ls(design, delta)`
#'
#' Computes restricted least-squares starting values for the nuisance
#' coefficients at a fixed treatment effect `delta`. The helper returns the
#' nuisance coefficient vector and errors only when the restricted fit itself
#' cannot produce finite values. Its role is to seed the non-linear optimizer in
#' the adjusted profile problem.
#'
#' ### `el_regression_moments(design, beta)`
#'
#' Evaluates the empirical-likelihood moment matrix `X_i (Y_i - X_i^T beta)` for
#' a given coefficient vector. The helper returns a numeric matrix with one row
#' per observed outcome and one column per regression parameter. Its role is to
#' translate the regression problem into the dual empirical-likelihood form.
#'
#' ### `el_dual_state(lambda, G)`
#'
#' Evaluates the dual empirical-likelihood objective and derivatives at
#' multiplier `lambda` for moment matrix `G`. The helper returns feasibility,
#' objective value, likelihood-ratio statistic, score, Hessian, and denominators
#' used in the weights. If any denominator is non-positive, it reports
#' infeasibility. Its role is to give both optimization routines a shared
#' representation of the dual problem.
#'
#' ### `el_dual_fit_newton(G, tol = 1e-8, maxit = 100)`
#'
#' Attempts to solve the dual problem by damped Newton iterations. The helper
#' returns a list containing success, feasibility, the multiplier, the
#' likelihood-ratio statistic, and diagnostic information. It can terminate
#' unsuccessfully when the Hessian system is singular or backtracking cannot find
#' a feasible improvement. Its role is to provide the fast first attempt before
#' falling back to a more forgiving optimizer.
#'
#' ### `el_dual_fit_nlminb(G, start = NULL, tol = 1e-8)`
#'
#' Alternative dual solver based on `nlminb()`. It returns the same essential
#' diagnostic fields as the Newton solver and declares success only when the
#' optimizer converges, the state remains feasible, and the score is sufficiently
#' small. Its role is to rescue problems that the Newton iteration cannot solve
#' reliably.
#'
#' ### `el_dual_fit(G, tol = 1e-8)`
#'
#' High-level dual solver used by the adjusted regression profile code. It checks
#' basic matrix validity, tries the Newton method, then falls back to `nlminb()`
#' if necessary. On success it also returns empirical-likelihood weights and the
#' weighted sample moments. Its role is to hide the two-stage optimization
#' strategy behind a single interface.
#'
#' ### `el_regression_profile_fit(design, delta, tol = 1e-8)`
#'
#' Profiles the adjusted empirical-likelihood regression problem at treatment
#' effect `delta`. The helper returns success, feasibility, the constrained
#' coefficient vector, nuisance coefficients, the profiled statistic, and the
#' dual fit. It short-circuits at the unrestricted least-squares estimate, uses
#' restricted least squares to initialize the nuisance optimizer, and reports
#' descriptive errors when profiling fails. Its role is to provide the core
#' objective used for adjusted `SPLRT` testing and confidence intervals.
#'
#' ### `el_regression_statistic(design, delta, tol = 1e-8)`
#'
#' Returns only the adjusted empirical-likelihood profile statistic at fixed
#' treatment effect `delta`. The helper converts profile failures to `Inf`.
#' Its role is to simplify confidence-interval inversion and testing code.
#'
#' ### `el_regression_variance(design)`
#'
#' Computes the least-squares variance estimate for the treatment coefficient
#' from the validated design object. The helper returns a scalar variance or
#' `NA_real_` when the residual variance or cross-product inverse is not
#' available. Its role is to provide a scale for the adjusted empirical-
#' likelihood interval search.
#'
#' ### `el_regression_ci_bound(design, estimate, crit, direction, se,
#' tol = 1e-8)`
#'
#' Finds one endpoint of the adjusted empirical-likelihood confidence interval.
#' `estimate` is the unrestricted treatment effect, `crit` is the chi-squared
#' cutoff, `direction` chooses the search side, and `se` is the least-squares
#' scale used to initialize the search radius. The helper returns a scalar bound
#' or `NA_real_`. It expands outward until the profile statistic exits the
#' acceptance region, then bisects back to the boundary. Its role is to invert
#' the adjusted profile statistic robustly.
#'
#' ### `el_regression_confint(design, estimate, conf.level, tol = 1e-8)`
#'
#' Builds the adjusted empirical-likelihood confidence interval by calling
#' `el_regression_ci_bound()` on both sides. It returns a length-two numeric
#' vector with a `"conf.level"` attribute or `c(NA, NA)` when inversion fails.
#' Its role is to provide the adjusted semi-parametric marginal interval.
#'
#' ### `el_regression_fit(data, formula, term = "R", mu = 0,
#' conf.level = 0.95, tol = 1e-8)`
#'
#' High-level adjusted empirical-likelihood wrapper. It validates `conf.level`,
#' constructs the design object, uses the closed-form unadjusted shortcut when
#' the design is only intercept plus treatment, and otherwise profiles the
#' treatment effect around `mu`. The helper returns a success/error list with
#' the estimate, confidence interval, statistic, p-value, design information,
#' and the profile object at the unrestricted estimate. Its role is to expose a
#' single adjusted empirical-likelihood API to the `SPLRT` wrapper.
#'
#' @seealso [SPLRT()], [adjusted_SPLRT()], [jointContrastCI()]
#' @keywords internal
NULL

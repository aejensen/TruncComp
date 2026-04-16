#' Internal parametric likelihood-ratio helpers
#'
#' Developer-facing documentation for the helpers that implement the parametric
#' Bernoulli-plus-Normal likelihood-ratio path used by `method = "LRT"`.
#'
#' @name trunccomp2-lrt-helpers
#' @title Internal parametric likelihood-ratio helpers
#' @description
#' These helpers fit the logistic and linear submodels, derive fallback
#' likelihood quantities in singular boundary cases, and assemble the unified
#' parametric fit object returned by [LRT()].
#' @aliases parametric_safe_xlogy
#' @aliases parametric_odds_ratio
#' @aliases parametric_wald_interval
#' @aliases parametric_loglik_bernoulli
#' @aliases parametric_model_formulas
#' @aliases parametric_fit_glm
#' @aliases parametric_fit_lm
#' @aliases parametric_fit_models
#' @aliases parametric_loglik_value
#' @aliases parametric_clamp_statistic
#' @aliases parametric_term_variance
#' @aliases parametric_term_estimate
#' @aliases parametric_term_interval
#' @aliases parametric_glm_is_regular
#' @aliases parametric_lm_is_regular
#' @aliases parametric_bernoulli_summary
#' @aliases parametric_normal_summary
#' @aliases parametric_adjusted_lrt_fit
#' @aliases parametric_lrt_fit
#' @aliases LRT
#' @usage
#' parametric_safe_xlogy(x, y)
#' parametric_odds_ratio(pi1, pi0)
#' parametric_wald_interval(estimate, se, conf.level)
#' parametric_loglik_bernoulli(k, n, p)
#' parametric_model_formulas(adjust = NULL)
#' parametric_fit_glm(formula, data)
#' parametric_fit_lm(formula, data)
#' parametric_fit_models(data, adjust = NULL)
#' parametric_loglik_value(fit, reml = NULL)
#' parametric_clamp_statistic(statistic, tol = 1e-12)
#' parametric_term_variance(fit, term)
#' parametric_term_estimate(fit, term)
#' parametric_term_interval(fit, term, conf.level, transform = identity)
#' parametric_glm_is_regular(fit, term = NULL, tol = 1e-8)
#' parametric_lm_is_regular(fit, term = NULL, tol = 1e-8)
#' parametric_bernoulli_summary(a, r, fits, conf.level, tol = 1e-12)
#' parametric_normal_summary(y0, y1, fits, conf.level, tol = 1e-12)
#' parametric_adjusted_lrt_fit(
#'   data, adjust, conf.level = 0.95, tol = 1e-12, regularity_tol = 1e-8
#' )
#' parametric_lrt_fit(
#'   data, conf.level = 0.95, adjust = NULL, tol = 1e-12,
#'   regularity_tol = 1e-8
#' )
#' LRT(
#'   data, init = NULL, conf.level = 0.95, adjust = NULL,
#'   adjust_spec = NULL, atom = NULL
#' )
#' @details
#' ### `parametric_safe_xlogy(x, y)`
#'
#' Computes `x * log(y)` elementwise while forcing structurally zero terms to
#' remain zero. `x` and `y` are numeric vectors, and the helper returns a numeric
#' vector of the same length. It avoids evaluating `log(0)` for cells where
#' `x == 0`. Its role is to support exact Bernoulli fallback likelihoods at the
#' observation-probability boundaries.
#'
#' ### `parametric_odds_ratio(pi1, pi0)`
#'
#' Converts two observation probabilities into the treatment-versus-control odds
#' ratio. `pi1` and `pi0` are scalar probabilities. The helper returns a scalar
#' numeric value, including exact `0`, `1`, or `Inf` in the saturated boundary
#' cases. Its role is to preserve intuitive odds-ratio output when the logistic
#' model fit is not regular enough for Wald intervals.
#'
#' ### `parametric_wald_interval(estimate, se, conf.level)`
#'
#' Builds a symmetric Wald interval from a point estimate and standard error.
#' The helper returns a length-two numeric vector. It returns `c(NA, NA)` when
#' the standard error is missing, non-finite, or negative. Its role is to keep
#' interval generation consistent across the logistic and linear components.
#'
#' ### `parametric_loglik_bernoulli(k, n, p)`
#'
#' Evaluates the Bernoulli log-likelihood for `k` successes in `n` trials at
#' probability `p`. The helper returns a scalar numeric log-likelihood and is
#' used only in fallback calculations. Because it delegates the problematic terms
#' to `parametric_safe_xlogy()`, it remains defined at exact `0` and `1`
#' probabilities. Its role is to support boundary-stable likelihood-ratio
#' statistics when `glm()` is not trustworthy.
#'
#' ### `parametric_model_formulas(adjust = NULL)`
#'
#' Constructs the four regression formulas needed by the parametric path.
#' `adjust` is either `NULL` or a validated one-sided additive adjustment
#' formula. The helper returns a named list containing Bernoulli null and
#' alternative formulas plus Normal null and alternative formulas. When
#' adjustment is present, treatment `R` is added only to the alternatives. Its
#' role is to keep the conditional model specification synchronized across both
#' submodels.
#'
#' ### `parametric_fit_glm(formula, data)` and `parametric_fit_lm(formula, data)`
#'
#' Thin wrappers around `glm()` and `lm()` that suppress warnings and convert fit
#' failures to `NULL`. Their arguments are the model `formula` and the analysis
#' `data` frame. Each helper returns either the fitted model object or `NULL`.
#' Their role is to let later regularity checks decide whether a fit is usable
#' without immediately aborting the full `LRT` pipeline.
#'
#' ### `parametric_fit_models(data, adjust = NULL)`
#'
#' Fits the full set of logistic and linear models used by the parametric
#' likelihood-ratio calculation. `data` is the standardized analysis frame, and
#' `adjust` is the optional additive baseline-covariate formula. The helper
#' returns a list containing the model formulas, the Bernoulli fits, the
#' observed-outcome data subset, and the Normal fits. It quietly carries through
#' `NULL` fits from the lower-level wrappers. Its role is to provide a single
#' cacheable bundle for all downstream summaries.
#'
#' ### `parametric_loglik_value(fit, reml = NULL)`
#'
#' Extracts a scalar numeric log-likelihood from a fitted model. `fit` may be a
#' `glm`, `lm`, or `NULL`, and `reml` controls whether `lm` log-likelihoods are
#' requested with `REML = FALSE`. The helper returns `NA_real_` when extraction
#' fails. Its role is to isolate the fragile `logLik()` calls used throughout
#' the parametric path.
#'
#' ### `parametric_clamp_statistic(statistic, tol = 1e-12)`
#'
#' Forces very small finite negative likelihood-ratio values back to zero.
#' `statistic` is a scalar numeric candidate and `tol` is the non-negativity
#' tolerance. The helper returns the adjusted scalar. It leaves large negative
#' or non-finite values untouched so callers can still detect genuine problems.
#' Its role is to clean up roundoff artifacts after log-likelihood subtraction.
#'
#' ### `parametric_term_variance(fit, term)` and
#' `parametric_term_estimate(fit, term)`
#'
#' Extract the variance or coefficient estimate for a named model term. `fit` is
#' a regression object and `term` is the coefficient name, typically `"R"`. Each
#' helper returns a scalar numeric value or `NA_real_` when the term is missing
#' or the covariance/coefficients cannot be extracted. Their role is to keep the
#' treatment-effect extraction logic consistent across Bernoulli and Normal
#' models.
#'
#' ### `parametric_term_interval(fit, term, conf.level, transform = identity)`
#'
#' Builds a treatment-term Wald interval from a fitted model. The helper extracts
#' the estimate and variance for `term`, applies `parametric_wald_interval()`,
#' and then transforms the interval, for example with `exp()` on the logistic
#' scale. It returns a length-two numeric vector or `c(NA, NA)` when the model
#' is not regular enough. Its role is to centralize the treatment-interval logic
#' used throughout the parametric summaries.
#'
#' ### `parametric_glm_is_regular(fit, term = NULL, tol = 1e-8)`
#'
#' Checks whether a logistic model fit is safe to use for adjusted inference.
#' `fit` is a `glm`, `term` optionally names a coefficient that must be present,
#' and `tol` controls the boundary thresholds. The helper returns `TRUE` only
#' when the model converged, the design is full rank, fitted probabilities stay
#' away from `0` and `1`, the log-likelihood is finite, and any requested term
#' has finite positive variance. Its role is to reject quasi-separation and rank
#' deficiency before adjusted `LRT` or adjusted `SPLRT` proceeds.
#'
#' ### `parametric_lm_is_regular(fit, term = NULL, tol = 1e-8)`
#'
#' Performs the corresponding regularity check for the observed-outcome `lm`
#' fits. The helper returns `TRUE` only when the design is full rank, the ML
#' log-likelihood is finite, the residual scale is positive, and any requested
#' treatment term has finite positive variance. Its role is to prevent adjusted
#' parametric inference from silently using aliased or exact-fit linear models.
#'
#' ### `parametric_bernoulli_summary(a, r, fits, conf.level, tol = 1e-12)`
#'
#' Summarizes the binary observation component for the unadjusted parametric
#' procedure. `a` and `r` are the standardized observation and treatment
#' vectors, `fits` is the model bundle from `parametric_fit_models()`, and
#' `conf.level` controls the logistic Wald interval when regular. The helper
#' returns a list containing counts, probabilities, log-likelihood values, the
#' Bernoulli likelihood-ratio statistic, the odds-ratio estimate, and its
#' interval. It falls back to exact cell-probability calculations when the model
#' fit is singular or boundary-limited. Its role is to make the unadjusted
#' Bernoulli component stable in the edge cases that motivated the older manual
#' implementation.
#'
#' ### `parametric_normal_summary(y0, y1, fits, conf.level, tol = 1e-12)`
#'
#' Summarizes the observed-outcome Normal component for the unadjusted
#' parametric procedure. `y0` and `y1` are the observed outcomes in the control
#' and treatment groups, and `fits` is the common model bundle. The helper
#' returns means, sums of squares, log-likelihood values, the Normal
#' likelihood-ratio statistic, the observed-outcome mean difference, and its
#' interval. It explicitly handles zero residual variance and uses an
#' sums-of-squares fallback when `logLik()` is unavailable. Its role is to keep
#' the Normal component numerically well-defined in exact-fit and near-boundary
#' settings.
#'
#' ### `parametric_adjusted_lrt_fit(data, adjust, conf.level = 0.95,
#' tol = 1e-12, regularity_tol = 1e-8)`
#'
#' Fits the adjusted parametric procedure. `data` is the standardized analysis
#' frame and `adjust` is the additive baseline-covariate formula. The helper
#' returns a list with `success`, component estimates and intervals, `Delta`
#' placeholders set to `NA`, the joint statistics, and nested fit summaries. It
#' returns `success = FALSE` with an error string when the logistic or linear
#' submodel is not regular, or when treatment-effect extraction fails. Its role
#' is to isolate the conditional-model implementation behind the same public
#' object structure used by the unadjusted path.
#'
#' ### `parametric_lrt_fit(data, conf.level = 0.95, adjust = NULL,
#' tol = 1e-12, regularity_tol = 1e-8)`
#'
#' Dispatches between the adjusted and unadjusted parametric procedures. The
#' helper returns the common intermediate fit list later consumed by [LRT()].
#' When `adjust` is `NULL`, it assembles the Bernoulli and Normal summaries for
#' the raw two-group comparison; otherwise it delegates to
#' `parametric_adjusted_lrt_fit()`. Its role is to keep the top-level `LRT()`
#' wrapper simple while still supporting both implementations.
#'
#' ### `LRT(data, init = NULL, conf.level = 0.95, adjust = NULL,
#' adjust_spec = NULL, atom = NULL)`
#'
#' Final parametric estimator used by `truncComp_core()`. `data` is the
#' standardized analysis frame; `init`, `conf.level`, `adjust`, `adjust_spec`,
#' and `atom` are carried through to the returned object. The helper returns a
#' `"TruncComp2"` object, either successful or failed. After a successful fit it
#' calls `augmentDeltaInference()` to fill in the unadjusted `Delta` summaries.
#' Its role is to bridge the intermediate parametric fit summaries into the
#' package-wide result class.
#'
#' @seealso [truncComp()], [augmentDeltaInference()], [parametric_lrt_fit()]
#' @keywords internal
NULL

#' Internal logistic profiling helpers
#'
#' Developer-facing documentation for the low-level logistic-profile routines
#' used by simultaneous confidence surfaces and `Delta` inference.
#'
#' @name trunccomp2-logit-helpers
#' @title Internal logistic profiling helpers
#' @description
#' These helpers build a cached two-parameter logistic representation and
#' profile the intercept for fixed treatment log-odds contrasts.
#' @aliases logit_log1pexp
#' @aliases logit_profile_score.prepared
#' @aliases logit_profile_interval.prepared
#' @aliases logit_profile_beta0.prepared
#' @aliases logit.prepare
#' @aliases logit.likelihood.prepared
#' @aliases logit.likelihood
#' @aliases logit.likelihood.profile.prepared
#' @aliases logit_profile_fit.prepared
#' @aliases logit.likelihood.profile
#' @aliases logit.LRT.prepared
#' @aliases logit.LRT
#' @usage
#' logit_log1pexp(x)
#' logit_profile_score.prepared(logitReference, beta0, delta)
#' logit_profile_interval.prepared(
#'   logitReference, delta, start = NULL, initial_width = 2, max_abs = 100
#' )
#' logit_profile_beta0.prepared(logitReference, delta, interval = NULL)
#' logit.prepare(data)
#' logit.likelihood.prepared(logitReference, beta)
#' logit.likelihood(data, beta)
#' logit.likelihood.profile.prepared(logitReference, delta, interval = NULL)
#' logit_profile_fit.prepared(logitReference, delta, interval = NULL)
#' logit.likelihood.profile(data, delta, interval = NULL)
#' logit.LRT.prepared(logitReference, delta, interval = NULL)
#' logit.LRT(data, delta)
#' @details
#' ### `logit_log1pexp(x)`
#'
#' Numerically stable helper for evaluating `log(1 + exp(x))`. `x` is a numeric
#' vector and the helper returns a numeric vector of the same length. It switches
#' formulas by sign to avoid overflow. Its role is to stabilize the custom
#' logistic log-likelihood calculations used outside `glm()`.
#'
#' ### `logit_profile_score.prepared(logitReference, beta0, delta)`
#'
#' Evaluates the score equation for the profiled intercept `beta0` when the
#' treatment log-odds contrast `delta` is fixed. `logitReference` is the cached
#' output from `logit.prepare()`. The helper returns the scalar intercept score.
#' It assumes the design has exactly intercept and treatment columns. Its role
#' is to support root-finding and interval expansion when profiling `beta0`.
#'
#' ### `logit_profile_interval.prepared(logitReference, delta, start = NULL,
#' initial_width = 2, max_abs = 100)`
#'
#' Expands a bracket that should contain the intercept root for a fixed `delta`.
#' `start` optionally supplies an initial intercept guess, `initial_width`
#' controls the first bracket size, and `max_abs` limits the search. The helper
#' returns a length-two numeric interval. It falls back to the empirical logit of
#' the marginal observation rate when no fitted intercept is available. Its role
#' is to give `uniroot()` a reasonable search region before the more expensive
#' likelihood fallback is tried.
#'
#' ### `logit_profile_beta0.prepared(logitReference, delta, interval = NULL)`
#'
#' Solves for the profiled intercept at fixed treatment contrast `delta`. It
#' returns the scalar intercept estimate. The helper prefers a score root via
#' `uniroot()` and falls back to direct optimization of the custom logistic
#' likelihood when the score bracket does not change sign. Its role is to keep
#' the logistic profiling stable even near the edges of the joint-inference
#' surface.
#'
#' ### `logit.prepare(data)`
#'
#' Constructs the cached reference object used by the other logistic helpers.
#' `data` is the standardized analysis frame. The helper returns a list
#' containing the model matrix, binary response, fitted `glm`, fitted
#' coefficients, and the alternative-model log-likelihood evaluated through the
#' custom likelihood routine. Its role is to avoid repeatedly rebuilding the
#' same design objects during surface evaluations.
#'
#' ### `logit.likelihood.prepared(logitReference, beta)` and
#' `logit.likelihood(data, beta)`
#'
#' Evaluate the custom negative logistic log-likelihood either from a prepared
#' reference object or directly from raw data. `beta` is the intercept-and-slope
#' vector. Each helper returns a scalar numeric value. Their role is to provide a
#' shared objective for the profiling and likelihood-ratio calculations.
#'
#' ### `logit.likelihood.profile.prepared(logitReference, delta, interval = NULL)`
#'
#' Profiles the logistic likelihood at a fixed treatment contrast `delta` by
#' first obtaining the best intercept with `logit_profile_beta0.prepared()`. The
#' helper returns the scalar profiled likelihood value. It inherits failures from
#' the intercept solver. Its role is to expose the profiled likelihood directly
#' for code that needs only the objective value.
#'
#' ### `logit_profile_fit.prepared(logitReference, delta, interval = NULL)`
#'
#' Builds the full profiled logistic fit for a fixed treatment log-odds
#' difference. The helper returns a list containing the profiled intercept,
#' fitted arm-specific probabilities, the profiled likelihood, and the
#' likelihood-ratio statistic relative to the unconstrained alternative model.
#' Its role is to supply the Bernoulli contribution for simultaneous confidence
#' surfaces and `Delta` profiling.
#'
#' ### `logit.likelihood.profile(data, delta, interval = NULL)`
#'
#' Convenience wrapper that constructs a temporary reference object and then
#' delegates to `logit.likelihood.profile.prepared()`. It returns the same scalar
#' profiled likelihood. Its role is mostly exploratory and developer-facing,
#' since the main package internals usually reuse cached references.
#'
#' ### `logit.LRT.prepared(logitReference, delta, interval = NULL)` and
#' `logit.LRT(data, delta)`
#'
#' Return the Bernoulli likelihood-ratio statistic for a fixed treatment
#' log-odds contrast, again either from a prepared reference or from raw data.
#' Each helper returns a scalar statistic. Their role is to provide the binary
#' component of the semi-parametric simultaneous region and the optimizer-based
#' `Delta` calculations.
#'
#' @seealso [logit.prepare()], [jointContrastCI()], [delta_profile_factory()]
#' @keywords internal
NULL

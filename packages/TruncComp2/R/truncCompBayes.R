#' Fit the experimental Bayesian two-sample TruncComp2 model
#'
#' Fits an explicit two-part Bayesian model for a binary treatment comparison
#' with a distinguished atom value representing an undefined or unobserved
#' outcome. The current implementation supports only the no-covariate model,
#' with either real-line non-atom outcomes or strictly positive non-atom
#' outcomes.
#'
#' @param formula For the formula interface, a formula with a continuous outcome
#'   on the left-hand side and a single binary treatment indicator on the
#'   right-hand side. For the default interface, the first positional argument
#'   is the outcome vector.
#' @param a Binary indicator for whether the continuous outcome is observed in
#'   the default interface.
#' @param r Binary treatment indicator for the default interface.
#' @param atom A single numeric atom value used for the distinguished outcome.
#'   For the default interface, `atom` may be omitted only when `y[a == 0]` has
#'   exactly one unique finite value.
#' @param data A data frame containing the variables referenced by `formula`.
#' @param conf.level Credible level used for stored posterior summaries.
#' @param continuous_support Support of the continuous non-atom component.
#'   `"real_line"` fits the existing Gaussian-mixture model on the real line,
#'   while `"positive_real"` fits a Gamma-mixture model on the positive real
#'   line and therefore requires all non-atom outcomes to be strictly
#'   positive.
#' @param mixture_components Truncation level for the arm-specific stick-breaking
#'   mixtures. Must be at least `2`.
#' @param chains Number of Stan chains.
#' @param iter_warmup Number of warmup iterations per chain.
#' @param iter_sampling Number of post-warmup iterations per chain.
#' @param seed Optional non-negative integer passed to Stan.
#' @param refresh Stan progress-refresh interval.
#' @param control Stan control list. The current implementation supports
#'   `adapt_delta` and `max_treedepth`.
#' @param prior Optional named list overriding the default priors. The common
#'   fields are `rho_alpha`, `rho_beta`, `alpha_shape`, and `alpha_rate`. For
#'   `continuous_support = "real_line"`, the support-specific fields are
#'   `mu_mean`, `mu_sd`, `sigma_meanlog`, and `sigma_sdlog`. For
#'   `continuous_support = "positive_real"`, the support-specific fields are
#'   `mean_meanlog`, `mean_sdlog`, `shape_meanlog`, and `shape_sdlog`.
#' @param ... Additional arguments passed to [rstan::sampling()]. Covariate
#'   adjustment is not supported in the Bayesian pathway.
#' @return An object of class `"trunc_comp_bayes_fit"` containing the fitted
#'   `stanfit`, posterior draws, stored posterior summaries, sampler diagnostics,
#'   the standardized analysis data, and the matched call. Failed fits return
#'   the same class with `success = FALSE` and an error message.
#' @details
#' The Bayesian pathway is experimental. It fits a two-part model with one atom
#' probability per treatment arm and one truncated stick-breaking mixture per
#' treatment arm for the non-atom outcomes. The continuous component can be
#' modeled either with Gaussian kernels on the real line or with Gamma kernels
#' on the positive real line. No p-values or likelihood-ratio statistics are
#' reported.
#'
#' In the stored posterior draws and summaries:
#'
#' - `rho_0` and `rho_1` are the arm-specific atom probabilities.
#' - `pi_0` and `pi_1` are the arm-specific probabilities of being observed
#'   away from the atom, so `pi_r = 1 - rho_r`.
#' - `mu_0_c` and `mu_1_c` are the arm-specific means among the non-atom
#'   outcomes.
#' - `delta_atom = rho_1 - rho_0` is the treatment-minus-control difference in
#'   atom probability.
#' - `mu_delta = mu_1^c - mu_0^c` is the treatment-minus-control difference in
#'   mean among the non-atom outcomes.
#' - `alpha_delta = [pi_1 / (1 - pi_1)] / [pi_0 / (1 - pi_0)]` is the
#'   treatment-to-control odds ratio for being observed away from the atom.
#'
#' The headline posterior summary is the combined-outcome contrast
#'
#' `delta = [atom * rho_1 + pi_1 * mu_1^c] - [atom * rho_0 + pi_0 * mu_0^c]`,
#'
#' which combines the atom probability and the non-atom mean into one
#' treatment-minus-control contrast on the original outcome scale.
#' @examples
#' \dontrun{
#' library(TruncComp2)
#' f0 <- function(n) stats::rnorm(n, 2, 1)
#' f1 <- function(n) stats::rnorm(n, 2.5, 1)
#' d <- simulate_truncated_data(20, f0 = f0, f1 = f1, pi0 = 0.6, pi1 = 0.7)
#' fit <- trunc_comp_bayes(Y ~ R, atom = 0, data = d, chains = 4,
#'                         iter_warmup = 500, iter_sampling = 1000)
#' summary(fit)
#' }
#' @rdname trunc_comp_bayes
#' @export
trunc_comp_bayes <- function(formula, ...) {
  UseMethod("trunc_comp_bayes")
}

trunc_comp_bayes_core <- function(y, a, r, atom, conf.level,
                                  continuous_support,
                                  mixture_components, chains,
                                  iter_warmup, iter_sampling,
                                  seed, refresh, control, prior,
                                  call = NULL, extra_args = list()) {
  if(!(is.numeric(y) && all(is.finite(y)))) {
    stop("y must contain only finite numeric outcomes.")
  }

  if(any(is.na(a)) || any(!(a %in% c(0, 1)))) {
    stop("a must contain only 0/1 values.")
  }

  if(any(is.na(r)) || any(!(r %in% c(0, 1)))) {
    stop("r must contain only 0/1 values.")
  }

  d <- data.frame(
    Y = as.numeric(y),
    A = as.integer(a),
    R = as.integer(r)
  )

  fit_trunc_comp_bayes(
    data = d,
    atom = atom,
    conf.level = conf.level,
    continuous_support = continuous_support,
    mixture_components = mixture_components,
    chains = chains,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    seed = seed,
    refresh = refresh,
    control = control,
    prior = prior,
    call = call,
    extra_args = extra_args
  )
}

#' @rdname trunc_comp_bayes
#' @export
trunc_comp_bayes.formula <- function(formula, atom, data,
                                     conf.level = 0.95,
                                     continuous_support = c("real_line", "positive_real"),
                                     mixture_components = 10,
                                     chains = 4,
                                     iter_warmup = 1000,
                                     iter_sampling = 1000,
                                     seed = NULL,
                                     refresh = 0,
                                     control = list(adapt_delta = 0.95, max_treedepth = 12),
                                     prior = NULL,
                                     ...) {
  if(!inherits(formula, "formula")) {
    stop("The formula must be a formula.")
  }

  if(length(attr(stats::terms(formula), "term.labels")) != 1) {
    stop("The current implementation must have one covariate in the formula.")
  }

  extra_args <- list(...)
  normalize_bayes_sampling_args(extra_args)

  outcome_name <- all.vars(formula[[2]])
  treatment_name <- all.vars(formula[[3]])

  if(length(outcome_name) != 1 || length(treatment_name) != 1) {
    stop("The formula must have a single outcome and a single binary treatment variable.")
  }

  variables <- stats::model.frame(formula, data = data, na.action = stats::na.fail)
  outcome <- variables[, 1]
  treatment <- variables[, 2]

  treatment_levels <- sort(unique(treatment))
  if(!(length(treatment_levels) == 2L && all(treatment_levels == c(0, 1)))) {
    stop("The covariate must be binary (0/1) indicating the two treatments.")
  }

  if(!(length(atom) == 1 && is.numeric(atom) && is.finite(atom))) {
    stop("atom must be a single finite numeric value.")
  }

  continuous_support <- bayes_continuous_support(continuous_support)

  alive <- as.integer(outcome != atom)

  if(all(alive[treatment == 0] == 0) || all(alive[treatment == 1] == 0)) {
    stop("Nothing has been observed in one of the groups. Cannot do estimation.")
  }

  if(all(alive == 1)) {
    stop("Everything seems to have been observed. You should use a different method.")
  }

  trunc_comp_bayes_core(
    y = outcome,
    a = alive,
    r = treatment,
    atom = as.numeric(atom),
    conf.level = conf.level,
    continuous_support = continuous_support,
    mixture_components = mixture_components,
    chains = chains,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    seed = seed,
    refresh = refresh,
    control = control,
    prior = prior,
    call = match.call(),
    extra_args = extra_args
  )
}

#' @rdname trunc_comp_bayes
#' @exportS3Method trunc_comp_bayes default
trunc_comp_bayes.default <- function(formula, a, r, atom = NULL,
                                     conf.level = 0.95,
                                     continuous_support = c("real_line", "positive_real"),
                                     mixture_components = 10,
                                     chains = 4,
                                     iter_warmup = 1000,
                                     iter_sampling = 1000,
                                     seed = NULL,
                                     refresh = 0,
                                     control = list(adapt_delta = 0.95, max_treedepth = 12),
                                     prior = NULL,
                                     ...) {
  y <- formula
  extra_args <- list(...)
  normalize_bayes_sampling_args(extra_args)

  if(length(y) != length(a) || length(y) != length(r)) {
    stop("y, a, and r must have the same length.")
  }

  stats::na.fail(data.frame(Y = y, A = a, R = r))

  atom <- resolveDefaultAtom(y, a, atom = atom)
  continuous_support <- bayes_continuous_support(continuous_support)

  if(any(a == 0 & y != atom, na.rm = TRUE) || any(a == 1 & y == atom, na.rm = TRUE)) {
    stop("For the default interface, y and a must agree with atom: y[a == 0] must equal atom and y[a == 1] must differ from atom.")
  }

  if(!all(sort(unique(r)) == c(0, 1))) {
    stop("r must be binary (0/1) indicating the two treatments.")
  }

  if(all(a[r == 0] == 0) || all(a[r == 1] == 0)) {
    stop("Nothing has been observed in one of the groups. Cannot do estimation.")
  }

  if(all(a == 1)) {
    stop("Everything seems to have been observed. You should use a different method.")
  }

  trunc_comp_bayes_core(
    y = y,
    a = a,
    r = r,
    atom = atom,
    conf.level = conf.level,
    continuous_support = continuous_support,
    mixture_components = mixture_components,
    chains = chains,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    seed = seed,
    refresh = refresh,
    control = control,
    prior = prior,
    call = match.call(),
    extra_args = extra_args
  )
}

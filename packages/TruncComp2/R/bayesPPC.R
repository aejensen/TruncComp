utils::globalVariables("replicated_mean")

bayes_with_seed <- function(seed, fn) {
  if(!is.function(fn)) {
    stop("fn must be a function.")
  }

  if(is.null(seed)) {
    return(fn())
  }

  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if(had_seed) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }

  on.exit(
    {
      if(had_seed) {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      } else if(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    },
    add = TRUE
  )

  set.seed(seed)
  fn()
}

bayes_validate_ppc_object <- function(object) {
  if(!inherits(object, "trunc_comp_bayes_fit")) {
    stop("object must be a trunc_comp_bayes_fit returned by trunc_comp_bayes().")
  }

  if(!isTRUE(object$success)) {
    stop("Estimation failed. Cannot compute posterior predictive checks.")
  }

  if(is.null(object$fit)) {
    stop("Raw Stan mixture parameters are not available.")
  }

  invisible(TRUE)
}

bayes_ppc_default_settings <- function(object) {
  list(
    ndraws = as.integer(min(200L, nrow(object$draws))),
    seed = validate_bayes_seed(object$settings$seed)
  )
}

bayes_ppc_resolve_settings <- function(object, ndraws = NULL, seed = NULL) {
  defaults <- bayes_ppc_default_settings(object)
  n_available <- validate_bayes_positive_integer(
    nrow(object$draws),
    "n_available",
    min_value = 1L
  )

  resolved_ndraws <- if(is.null(ndraws)) defaults$ndraws else {
    min(
      n_available,
      validate_bayes_positive_integer(ndraws, "ndraws", min_value = 1L)
    )
  }

  list(
    ndraws = as.integer(resolved_ndraws),
    seed = if(is.null(seed)) defaults$seed else validate_bayes_seed(seed)
  )
}

bayes_ppc_same_settings <- function(x, y) {
  identical(as.integer(x$ndraws), as.integer(y$ndraws)) &&
    identical(validate_bayes_seed(x$seed), validate_bayes_seed(y$seed))
}

bayes_ppc_draw_indices <- function(n_available, ndraws, seed = NULL) {
  n_available <- validate_bayes_positive_integer(n_available, "n_available", min_value = 1L)
  ndraws <- validate_bayes_positive_integer(ndraws, "ndraws", min_value = 1L)
  n_select <- min(n_available, ndraws)

  if(is.null(seed)) {
    return(seq_len(n_select))
  }

  sample.int(n_available, size = n_select, replace = FALSE)
}

bayes_ppc_component_parameters <- function(object) {
  extracted <- bayes_density_component_draws(object)

  list(
    support = extracted$support,
    weights = extracted$weights,
    means = extracted$means,
    sds = extracted$sds,
    shapes = extracted$shapes,
    m_comp = extracted$m_comp,
    phi_comp = extracted$phi_comp,
    eta = extracted$eta,
    score_min = extracted$score_min,
    score_max = extracted$score_max,
    score_values = extracted$score_values,
    eta_group_by_arm = extracted$eta_group_by_arm,
    bin_lower = extracted$bin_lower,
    bin_upper = extracted$bin_upper,
    bin_valid = extracted$bin_valid,
    rho = extracted$rho
  )
}

bayes_ppc_atom_matrix <- function(rho, arm, draw_indices) {
  n_rep <- length(draw_indices)
  n_obs <- length(arm)
  yrep <- matrix(0L, nrow = n_rep, ncol = n_obs)

  for(s in seq_len(n_rep)) {
    draw <- draw_indices[s]
    yrep[s, ] <- stats::rbinom(n_obs, size = 1, prob = rho[draw, arm])
  }

  yrep
}

bayes_ppc_eta_draw <- function(eta, draw) {
  if(is.null(eta)) {
    return(NULL)
  }

  out <- eta[draw, , , drop = FALSE]
  out <- array(out, dim = dim(out)[-1])
  if(is.null(dim(out))) {
    out <- matrix(out, nrow = 1L)
  }
  if(length(dim(out)) == 1L) {
    out <- matrix(out, nrow = 1L)
  }

  out
}

bayes_ppc_continuous_matrix <- function(weights, means, arm_obs, draw_indices,
                                        continuous_support = c(
                                          "real_line",
                                          "positive_real",
                                          "bounded_continuous",
                                          "bounded_score"
                                        ),
                                        sds = NULL, shapes = NULL,
                                        m_comp = NULL, phi_comp = NULL,
                                        eta = NULL,
                                        score_min = NULL, score_max = NULL,
                                        score_values = NULL,
                                        eta_group_by_arm = NULL,
                                        bin_lower = NULL, bin_upper = NULL,
                                        bin_valid = NULL) {
  continuous_support <- bayes_continuous_support(continuous_support)
  n_rep <- length(draw_indices)
  n_obs <- length(arm_obs)
  n_components <- dim(weights)[3]
  yrep <- matrix(NA_real_, nrow = n_rep, ncol = n_obs)

  for(s in seq_len(n_rep)) {
    draw <- draw_indices[s]

    for(arm in 1:2) {
      arm_index <- which(arm_obs == arm)
      if(length(arm_index) == 0L) {
        next
      }

      component_index <- sample.int(
        n_components,
        size = length(arm_index),
        replace = TRUE,
        prob = weights[draw, arm, ]
      )

      if(identical(continuous_support, "real_line")) {
        yrep[s, arm_index] <- stats::rnorm(
          length(arm_index),
          mean = means[draw, arm, component_index],
          sd = sds[draw, arm, component_index]
        )
      } else if(identical(continuous_support, "positive_real")) {
        component_shapes <- shapes[draw, arm, component_index]
        component_means <- means[draw, arm, component_index]
        yrep[s, arm_index] <- stats::rgamma(
          length(arm_index),
          shape = component_shapes,
          rate = component_shapes / component_means
        )
      } else if(identical(continuous_support, "bounded_continuous")) {
        component_m <- m_comp[draw, arm, component_index]
        component_phi <- phi_comp[draw, arm, component_index]
        x_rep <- stats::rbeta(
          length(arm_index),
          shape1 = component_m * component_phi,
          shape2 = (1 - component_m) * component_phi
        )
        yrep[s, arm_index] <- score_min + (score_max - score_min) * x_rep
      } else {
        eta_group <- eta_group_by_arm[[arm]]
        pmf <- bayes_score_pmf(
          weights = weights[draw, arm, ],
          m_comp = m_comp[draw, arm, ],
          phi_comp = phi_comp[draw, arm, ],
          eta = eta[draw, eta_group, ],
          bin_lower = bin_lower,
          bin_upper = bin_upper,
          bin_valid = bin_valid
        )
        yrep[s, arm_index] <- sample(
          score_values,
          size = length(arm_index),
          replace = TRUE,
          prob = pmf
        )
      }
    }
  }

  yrep
}

bayes_ppc_data <- function(object, ndraws = 50L, seed = NULL, parameters = NULL) {
  bayes_validate_ppc_object(object)

  ndraws <- validate_bayes_positive_integer(ndraws, "ndraws", min_value = 1L)
  seed <- validate_bayes_seed(seed)

  data <- object$data
  if(is.null(parameters)) {
    parameters <- bayes_ppc_component_parameters(object)
  }
  arm_full <- as.integer(data$R) + 1L
  arm_obs <- arm_full[data$A == 1]
  arm_labels <- c("Control", "Treatment")

  bayes_with_seed(seed, function() {
    draw_indices <- bayes_ppc_draw_indices(
      n_available = dim(parameters$rho)[1],
      ndraws = ndraws,
      seed = seed
    )

    y_atom <- 1L - as.integer(data$A)
    yrep_atom <- bayes_ppc_atom_matrix(
      rho = parameters$rho,
      arm = arm_full,
      draw_indices = draw_indices
    )

    y_cont <- as.numeric(data$Y[data$A == 1])
    yrep_cont <- bayes_ppc_continuous_matrix(
      weights = parameters$weights,
      means = parameters$means,
      arm_obs = arm_obs,
      draw_indices = draw_indices,
      continuous_support = parameters$support,
      sds = parameters$sds,
      shapes = parameters$shapes,
      m_comp = parameters$m_comp,
      phi_comp = parameters$phi_comp,
      eta = parameters$eta,
      score_min = parameters$score_min,
      score_max = parameters$score_max,
      score_values = parameters$score_values,
      eta_group_by_arm = parameters$eta_group_by_arm,
      bin_lower = parameters$bin_lower,
      bin_upper = parameters$bin_upper,
      bin_valid = parameters$bin_valid
    )

    list(
      draw_indices = draw_indices,
      y_atom = y_atom,
      yrep_atom = yrep_atom,
      arm_full = arm_full,
      group_atom = factor(arm_labels[arm_full], levels = arm_labels),
      y_cont = y_cont,
      yrep_cont = yrep_cont,
      arm_obs = arm_obs,
      group_cont = factor(arm_labels[arm_obs], levels = arm_labels)
    )
  })
}

bayes_ppc_continuous_scale_label <- function(continuous_support = c(
  "real_line",
  "positive_real",
  "bounded_continuous",
  "bounded_score"
)) {
  continuous_support <- bayes_continuous_support(continuous_support)

  if(identical(continuous_support, "positive_real")) {
    return("log(Outcome)")
  }

  if(identical(continuous_support, "bounded_continuous")) {
    return("Bounded outcome")
  }

  if(identical(continuous_support, "bounded_score")) {
    return("Reported score")
  }

  "Outcome"
}

bayes_ppc_continuous_transform <- function(x,
                                           continuous_support = c(
                                             "real_line",
                                             "positive_real",
                                             "bounded_continuous",
                                             "bounded_score"
                                           )) {
  continuous_support <- bayes_continuous_support(continuous_support)

  if(identical(continuous_support, "positive_real")) {
    return(log(pmax(x, .Machine$double.xmin)))
  }

  x
}

bayes_ppc_continuous_plot_inputs <- function(ppc_data,
                                             continuous_support = c(
                                               "real_line",
                                               "positive_real",
                                               "bounded_continuous",
                                               "bounded_score"
                                             )) {
  continuous_support <- bayes_continuous_support(continuous_support)

  list(
    y = bayes_ppc_continuous_transform(
      ppc_data$y_cont,
      continuous_support = continuous_support
    ),
    yrep = bayes_ppc_continuous_transform(
      ppc_data$yrep_cont,
      continuous_support = continuous_support
    ),
    x_label = bayes_ppc_continuous_scale_label(continuous_support)
  )
}

bayes_ppc_atom_discrepancy <- function(y_atom, arm_full, rho_draw) {
  arm_means <- vapply(
    1:2,
    function(arm) mean(y_atom[arm_full == arm]),
    numeric(1)
  )

  max(abs(arm_means - rho_draw))
}

bayes_ppc_mixture_cdf <- function(x, weights, means,
                                  continuous_support = c(
                                    "real_line",
                                    "positive_real",
                                    "bounded_continuous",
                                    "bounded_score"
                                  ),
                                  sds = NULL, shapes = NULL,
                                  m_comp = NULL, phi_comp = NULL,
                                  score_min = NULL, score_max = NULL) {
  continuous_support <- bayes_continuous_support(continuous_support)

  kernels <- vapply(
    seq_along(weights),
    function(h) {
      if(identical(continuous_support, "real_line")) {
        return(stats::pnorm(x, mean = means[h], sd = sds[h]))
      }

      if(identical(continuous_support, "bounded_continuous")) {
        x_unit <- pmin(pmax((x - score_min) / (score_max - score_min), 0), 1)
        return(stats::pbeta(
          x_unit,
          shape1 = m_comp[h] * phi_comp[h],
          shape2 = (1 - m_comp[h]) * phi_comp[h]
        ))
      }

      stats::pgamma(
        exp(x),
        shape = shapes[h],
        rate = shapes[h] / means[h]
      )
    },
    numeric(length(x))
  )

  as.numeric(kernels %*% weights)
}

bayes_ppc_ks_discrepancy <- function(sample, fitted_cdf_values) {
  n <- length(sample)
  if(n == 0L) {
    return(0)
  }

  i <- seq_len(n)
  max(
    max(i / n - fitted_cdf_values),
    max(fitted_cdf_values - (i - 1) / n)
  )
}

bayes_ppc_continuous_discrepancy <- function(y_cont, arm_obs, weights_draw, means_draw,
                                             continuous_support = c(
                                               "real_line",
                                               "positive_real",
                                               "bounded_continuous",
                                               "bounded_score"
                                             ),
                                             sds_draw = NULL, shapes_draw = NULL,
                                             m_comp_draw = NULL, phi_comp_draw = NULL,
                                             eta_draw = NULL,
                                             score_min = NULL, score_max = NULL,
                                             score_values = NULL,
                                             eta_group_by_arm = NULL,
                                             bin_lower = NULL, bin_upper = NULL,
                                             bin_valid = NULL) {
  continuous_support <- bayes_continuous_support(continuous_support)

  max(vapply(
    1:2,
    function(arm) {
      arm_sample <- y_cont[arm_obs == arm]
      if(length(arm_sample) == 0L) {
        return(0)
      }

      if(identical(continuous_support, "bounded_score")) {
        eta_group <- eta_group_by_arm[[arm]]
        pmf <- bayes_score_pmf(
          weights = weights_draw[arm, ],
          m_comp = m_comp_draw[arm, ],
          phi_comp = phi_comp_draw[arm, ],
          eta = eta_draw[eta_group, ],
          bin_lower = bin_lower,
          bin_upper = bin_upper,
          bin_valid = bin_valid
        )
        fitted_cdf <- cumsum(pmf)
        observed_cdf <- vapply(
          score_values,
          function(value) mean(arm_sample <= value),
          numeric(1)
        )

        return(max(abs(observed_cdf - fitted_cdf)))
      }

      transformed_sample <- sort(bayes_ppc_continuous_transform(
        arm_sample,
        continuous_support = continuous_support
      ))
      fitted_cdf <- bayes_ppc_mixture_cdf(
        x = transformed_sample,
        weights = weights_draw[arm, ],
        means = if(continuous_support %in% c("real_line", "positive_real")) means_draw[arm, ] else NULL,
        continuous_support = continuous_support,
        sds = if(identical(continuous_support, "real_line")) sds_draw[arm, ] else NULL,
        shapes = if(identical(continuous_support, "positive_real")) shapes_draw[arm, ] else NULL,
        m_comp = if(identical(continuous_support, "bounded_continuous")) m_comp_draw[arm, ] else NULL,
        phi_comp = if(identical(continuous_support, "bounded_continuous")) phi_comp_draw[arm, ] else NULL,
        score_min = if(identical(continuous_support, "bounded_continuous")) score_min else NULL,
        score_max = if(identical(continuous_support, "bounded_continuous")) score_max else NULL
      )

      bayes_ppc_ks_discrepancy(transformed_sample, fitted_cdf)
    },
    numeric(1)
  ))
}

bayes_ppc_table <- function(atom_obs_discrepancy,
                            atom_rep_discrepancy,
                            continuous_obs_discrepancy,
                            continuous_rep_discrepancy,
                            continuous_support = c(
                              "real_line",
                              "positive_real",
                              "bounded_continuous",
                              "bounded_score"
                            )) {
  continuous_support <- bayes_continuous_support(continuous_support)
  continuous_statistic <- switch(
    continuous_support,
    bounded_continuous = "Maximum armwise bounded CDF discrepancy",
    bounded_score = "Maximum armwise discrete CDF discrepancy",
    "Maximum armwise KS discrepancy"
  )

  out <- data.frame(
    component = c("atom", "continuous"),
    p_value = c(
      mean(atom_rep_discrepancy >= atom_obs_discrepancy),
      mean(continuous_rep_discrepancy >= continuous_obs_discrepancy)
    ),
    statistic = c(
      "Maximum armwise absolute atom-rate deviation",
      continuous_statistic
    ),
    scale = c(
      "Indicator for Y = atom",
      bayes_ppc_continuous_scale_label(continuous_support)
    ),
    ndraws = c(
      length(atom_obs_discrepancy),
      length(continuous_obs_discrepancy)
    ),
    mean_observed_discrepancy = c(
      mean(atom_obs_discrepancy),
      mean(continuous_obs_discrepancy)
    ),
    mean_replicated_discrepancy = c(
      mean(atom_rep_discrepancy),
      mean(continuous_rep_discrepancy)
    ),
    row.names = NULL
  )

  rownames(out) <- out$component
  out
}

bayes_ppc_summary <- function(object, ndraws = NULL, seed = NULL,
                              ppc_data = NULL, parameters = NULL) {
  bayes_validate_ppc_object(object)

  settings <- bayes_ppc_resolve_settings(object, ndraws = ndraws, seed = seed)
  cached_default <- is.null(ppc_data) &&
    is.null(parameters) &&
    !is.null(object$ppc_table) &&
    !is.null(object$ppc_settings) &&
    bayes_ppc_same_settings(object$ppc_settings, settings)

  if(cached_default) {
    return(list(
      table = object$ppc_table,
      settings = object$ppc_settings
    ))
  }

  if(is.null(parameters)) {
    parameters <- bayes_ppc_component_parameters(object)
  }

  if(is.null(ppc_data)) {
    ppc_data <- bayes_ppc_data(
      object = object,
      ndraws = settings$ndraws,
      seed = settings$seed,
      parameters = parameters
    )
  }

  n_rep <- length(ppc_data$draw_indices)
  atom_obs_discrepancy <- numeric(n_rep)
  atom_rep_discrepancy <- numeric(n_rep)
  continuous_obs_discrepancy <- numeric(n_rep)
  continuous_rep_discrepancy <- numeric(n_rep)

  for(s in seq_len(n_rep)) {
    draw <- ppc_data$draw_indices[s]
    rho_draw <- parameters$rho[draw, ]

    atom_obs_discrepancy[s] <- bayes_ppc_atom_discrepancy(
      y_atom = ppc_data$y_atom,
      arm_full = ppc_data$arm_full,
      rho_draw = rho_draw
    )
    atom_rep_discrepancy[s] <- bayes_ppc_atom_discrepancy(
      y_atom = ppc_data$yrep_atom[s, ],
      arm_full = ppc_data$arm_full,
      rho_draw = rho_draw
    )

    continuous_obs_discrepancy[s] <- bayes_ppc_continuous_discrepancy(
      y_cont = ppc_data$y_cont,
      arm_obs = ppc_data$arm_obs,
      weights_draw = parameters$weights[draw, , ],
      means_draw = if(parameters$support %in% c("real_line", "positive_real")) parameters$means[draw, , ] else NULL,
      continuous_support = parameters$support,
      sds_draw = if(identical(parameters$support, "real_line")) parameters$sds[draw, , ] else NULL,
      shapes_draw = if(identical(parameters$support, "positive_real")) parameters$shapes[draw, , ] else NULL,
      m_comp_draw = if(parameters$support %in% c("bounded_continuous", "bounded_score")) parameters$m_comp[draw, , ] else NULL,
      phi_comp_draw = if(parameters$support %in% c("bounded_continuous", "bounded_score")) parameters$phi_comp[draw, , ] else NULL,
      eta_draw = if(identical(parameters$support, "bounded_score")) bayes_ppc_eta_draw(parameters$eta, draw) else NULL,
      score_min = parameters$score_min,
      score_max = parameters$score_max,
      score_values = parameters$score_values,
      eta_group_by_arm = parameters$eta_group_by_arm,
      bin_lower = parameters$bin_lower,
      bin_upper = parameters$bin_upper,
      bin_valid = parameters$bin_valid
    )
    continuous_rep_discrepancy[s] <- bayes_ppc_continuous_discrepancy(
      y_cont = ppc_data$yrep_cont[s, ],
      arm_obs = ppc_data$arm_obs,
      weights_draw = parameters$weights[draw, , ],
      means_draw = if(parameters$support %in% c("real_line", "positive_real")) parameters$means[draw, , ] else NULL,
      continuous_support = parameters$support,
      sds_draw = if(identical(parameters$support, "real_line")) parameters$sds[draw, , ] else NULL,
      shapes_draw = if(identical(parameters$support, "positive_real")) parameters$shapes[draw, , ] else NULL,
      m_comp_draw = if(parameters$support %in% c("bounded_continuous", "bounded_score")) parameters$m_comp[draw, , ] else NULL,
      phi_comp_draw = if(parameters$support %in% c("bounded_continuous", "bounded_score")) parameters$phi_comp[draw, , ] else NULL,
      eta_draw = if(identical(parameters$support, "bounded_score")) bayes_ppc_eta_draw(parameters$eta, draw) else NULL,
      score_min = parameters$score_min,
      score_max = parameters$score_max,
      score_values = parameters$score_values,
      eta_group_by_arm = parameters$eta_group_by_arm,
      bin_lower = parameters$bin_lower,
      bin_upper = parameters$bin_upper,
      bin_valid = parameters$bin_valid
    )
  }

  list(
    table = bayes_ppc_table(
      atom_obs_discrepancy = atom_obs_discrepancy,
      atom_rep_discrepancy = atom_rep_discrepancy,
      continuous_obs_discrepancy = continuous_obs_discrepancy,
      continuous_rep_discrepancy = continuous_rep_discrepancy,
      continuous_support = parameters$support
    ),
    settings = list(
      ndraws = as.integer(n_rep),
      seed = settings$seed
    ),
    details = list(
      atom_observed_discrepancy = atom_obs_discrepancy,
      atom_replicated_discrepancy = atom_rep_discrepancy,
      continuous_observed_discrepancy = continuous_obs_discrepancy,
      continuous_replicated_discrepancy = continuous_rep_discrepancy
    )
  )
}

bayes_ppc_score_plot <- function(ppc_data, score_values, conf.level, p_value) {
  arm_labels <- levels(ppc_data$group_cont)
  rows <- vector("list", length(arm_labels))

  for(r in seq_along(arm_labels)) {
    arm_index <- which(ppc_data$group_cont == arm_labels[[r]])
    observed <- ppc_data$y_cont[arm_index]
    replicated <- ppc_data$yrep_cont[, arm_index, drop = FALSE]

    replicated_mass <- t(apply(
      replicated,
      1L,
      function(draw_values) {
        vapply(score_values, function(value) mean(draw_values == value), numeric(1))
      }
    ))
    intervals <- t(vapply(
      seq_along(score_values),
      function(index) bayes_equal_tail_interval(replicated_mass[, index], conf.level),
      numeric(2)
    ))

    rows[[r]] <- data.frame(
      arm_label = factor(arm_labels[[r]], levels = arm_labels),
      x = score_values,
      observed = vapply(score_values, function(value) mean(observed == value), numeric(1)),
      replicated_mean = colMeans(replicated_mass),
      conf.low = intervals[, 1],
      conf.high = intervals[, 2],
      row.names = NULL
    )
  }

  plot_data <- do.call(rbind, rows)

  ggplot2::ggplot(plot_data, ggplot2::aes(x = x)) +
    ggplot2::geom_col(
      ggplot2::aes(y = observed),
      fill = "grey82",
      width = min(diff(unique(score_values)), na.rm = TRUE) * 0.8
    ) +
    ggplot2::geom_linerange(
      ggplot2::aes(ymin = conf.low, ymax = conf.high),
      linewidth = 0.45,
      color = "#1f78b4"
    ) +
    ggplot2::geom_point(
      ggplot2::aes(y = replicated_mean),
      size = 1.5,
      color = "#1f78b4"
    ) +
    ggplot2::facet_wrap(ggplot2::vars(arm_label), scales = "fixed") +
    ggplot2::labs(
      title = "Posterior predictive check for the continuous part",
      subtitle = sprintf("Posterior predictive p = %.3f", p_value),
      x = "Reported score",
      y = "Probability"
    ) +
    ggplot2::theme_minimal()
}

#' Posterior predictive p-values for a Bayesian truncated-comparison fit
#'
#' Computes discrepancy-based posterior predictive p-values for the atom and
#' continuous parts of a successful Bayesian two-part fit. The atom p-value is
#' based on the maximum armwise absolute deviation between the observed atom
#' rate and the fitted atom probability. The continuous p-value is based on the
#' maximum armwise Kolmogorov-Smirnov discrepancy between the empirical CDF of
#' the observed non-atom outcomes and the fitted continuous mixture CDF. These
#' are posterior predictive model-checking summaries, not frequentist
#' hypothesis-test p-values.
#'
#' @param object A successful `"trunc_comp_bayes_fit"` object returned by
#'   [trunc_comp_bayes()].
#' @param ndraws Optional number of posterior draws to use when simulating
#'   predictive replications. If `NULL`, uses the cached default summary based
#'   on `min(200, n_posterior_draws)` draws.
#' @param seed Optional non-negative integer used to make the predictive
#'   p-values reproducible. If `NULL`, defaults to the Stan seed stored in the
#'   fit when available.
#' @return A data frame with one row for the atom PPC and one row for the
#'   continuous PPC. The returned columns are `p_value`, `statistic`, `scale`,
#'   `ndraws`, `mean_observed_discrepancy`, and
#'   `mean_replicated_discrepancy`.
#' @examples
#' \dontrun{
#' data("trunc_comp_example", package = "TruncComp2")
#' fit <- trunc_comp_bayes(
#'   Y ~ R,
#'   atom = 0,
#'   data = trunc_comp_example,
#'   chains = 4,
#'   iter_warmup = 500,
#'   iter_sampling = 1000,
#'   refresh = 0
#' )
#'
#' posterior_predictive_pvalues(fit)
#' }
#' @export
posterior_predictive_pvalues <- function(object, ndraws = NULL, seed = NULL) {
  bayes_ppc_summary(object = object, ndraws = ndraws, seed = seed)$table
}

#' Posterior predictive checks for a Bayesian truncated-comparison fit
#'
#' Uses [bayesplot](https://mc-stan.org/bayesplot/) to visualize posterior
#' predictive checks for the experimental Bayesian two-part model. The atom PPC
#' checks the Bernoulli atom model across all subjects, while the continuous PPC
#' checks the conditional non-atom mixture model among the observed outcomes.
#' For `continuous_support = "real_line"` this uses Gaussian predictive draws,
#' while `continuous_support = "positive_real"` uses Gamma predictive draws on
#' the original positive outcome scale and displays the continuous PPC as a
#' density overlay on `log(Y)` to avoid boundary artifacts at zero. Bounded
#' continuous fits use Beta-mixture predictive draws on the bounded outcome
#' scale. Bounded-score fits use discrete reported-score replications and a
#' score-mass PPC rather than a smooth density overlay.
#'
#' @param object A successful `"trunc_comp_bayes_fit"` object returned by
#'   [trunc_comp_bayes()].
#' @param type Which posterior predictive checks to return: both plots, only the
#'   atom PPC, or only the continuous PPC.
#' @param ndraws Number of posterior draws to use when simulating predictive
#'   replications. If `ndraws` exceeds the available number of posterior draws,
#'   it is clipped to the available size.
#' @param seed Optional non-negative integer used to make the predictive checks
#'   reproducible.
#' @return If `type = "both"`, a named list with `atom` and `continuous`
#'   `ggplot2` objects. Otherwise returns the requested `ggplot2` object. Each
#'   plot subtitle reports the corresponding posterior predictive p-value.
#' @examples
#' \dontrun{
#' data("trunc_comp_example", package = "TruncComp2")
#' fit <- trunc_comp_bayes(
#'   Y ~ R,
#'   atom = 0,
#'   data = trunc_comp_example,
#'   chains = 4,
#'   iter_warmup = 500,
#'   iter_sampling = 1000,
#'   refresh = 0
#' )
#'
#' ppc <- posterior_predictive_check(fit, seed = 1)
#' ppc$atom
#' ppc$continuous
#' }
#' @export
posterior_predictive_check <- function(object,
                                       type = c("both", "atom", "continuous"),
                                       ndraws = 50,
                                       seed = NULL) {
  type <- match.arg(type)
  bayes_validate_ppc_object(object)
  resolved_settings <- bayes_ppc_resolve_settings(object, ndraws = ndraws, seed = seed)
  parameters <- bayes_ppc_component_parameters(object)
  ppc_data <- bayes_ppc_data(
    object = object,
    ndraws = resolved_settings$ndraws,
    seed = resolved_settings$seed,
    parameters = parameters
  )
  ppc_summary <- bayes_ppc_summary(
    object = object,
    ndraws = resolved_settings$ndraws,
    seed = resolved_settings$seed,
    ppc_data = ppc_data,
    parameters = parameters
  )
  ppc_table <- ppc_summary$table
  continuous_support <- bayes_fit_continuous_support(object)
  continuous_plot_inputs <- bayes_ppc_continuous_plot_inputs(
    ppc_data,
    continuous_support = continuous_support
  )

  atom_plot <- bayesplot::ppc_bars_grouped(
    y = ppc_data$y_atom,
    yrep = ppc_data$yrep_atom,
    group = ppc_data$group_atom,
    prob = object$conf.level,
    freq = FALSE
  ) +
    ggplot2::labs(
      title = "Posterior predictive check for the atom model",
      subtitle = sprintf(
        "Posterior predictive p = %.3f",
        ppc_table["atom", "p_value"]
      ),
      x = paste0("Indicator for Y = atom (", object$atom, ")"),
      y = "Proportion"
    )

  continuous_plot <- if(identical(continuous_support, "bounded_score")) {
    bayes_ppc_score_plot(
      ppc_data = ppc_data,
      score_values = parameters$score_values,
      conf.level = object$conf.level,
      p_value = ppc_table["continuous", "p_value"]
    )
  } else {
    out <- bayesplot::ppc_dens_overlay_grouped(
      y = continuous_plot_inputs$y,
      yrep = continuous_plot_inputs$yrep,
      group = ppc_data$group_cont
    ) +
      ggplot2::labs(
        title = "Posterior predictive check for the continuous part",
        subtitle = sprintf(
          "Posterior predictive p = %.3f",
          ppc_table["continuous", "p_value"]
        ),
        x = continuous_plot_inputs$x_label,
        y = "Density"
      )

    if(identical(continuous_support, "bounded_continuous")) {
      out <- out + ggplot2::coord_cartesian(xlim = c(parameters$score_min, parameters$score_max))
    }

    out
  }

  if(identical(type, "atom")) {
    return(atom_plot)
  }

  if(identical(type, "continuous")) {
    return(continuous_plot)
  }

  list(
    atom = atom_plot,
    continuous = continuous_plot
  )
}

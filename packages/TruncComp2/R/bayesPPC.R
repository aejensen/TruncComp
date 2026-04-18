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

bayes_ppc_continuous_matrix <- function(weights, means, arm_obs, draw_indices,
                                        continuous_support = c("real_line", "positive_real"),
                                        sds = NULL, shapes = NULL) {
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
      } else {
        component_shapes <- shapes[draw, arm, component_index]
        component_means <- means[draw, arm, component_index]
        yrep[s, arm_index] <- stats::rgamma(
          length(arm_index),
          shape = component_shapes,
          rate = component_shapes / component_means
        )
      }
    }
  }

  yrep
}

bayes_ppc_data <- function(object, ndraws = 50L, seed = NULL) {
  if(!inherits(object, "trunc_comp_bayes_fit")) {
    stop("object must be a trunc_comp_bayes_fit returned by trunc_comp_bayes().")
  }

  if(!isTRUE(object$success)) {
    stop("Estimation failed. Cannot compute posterior predictive checks.")
  }

  ndraws <- validate_bayes_positive_integer(ndraws, "ndraws", min_value = 1L)
  seed <- validate_bayes_seed(seed)

  data <- object$data
  parameters <- bayes_ppc_component_parameters(object)
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
      shapes = parameters$shapes
    )

    list(
      draw_indices = draw_indices,
      y_atom = y_atom,
      yrep_atom = yrep_atom,
      group_atom = factor(arm_labels[arm_full], levels = arm_labels),
      y_cont = y_cont,
      yrep_cont = yrep_cont,
      group_cont = factor(arm_labels[arm_obs], levels = arm_labels)
    )
  })
}

bayes_ppc_continuous_plot_inputs <- function(ppc_data,
                                             continuous_support = c("real_line", "positive_real")) {
  continuous_support <- bayes_continuous_support(continuous_support)

  if(identical(continuous_support, "positive_real")) {
    return(list(
      y = log(pmax(ppc_data$y_cont, .Machine$double.xmin)),
      yrep = log(pmax(ppc_data$yrep_cont, .Machine$double.xmin)),
      x_label = "log(Outcome)"
    ))
  }

  list(
    y = ppc_data$y_cont,
    yrep = ppc_data$yrep_cont,
    x_label = "Outcome"
  )
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
#' density overlay on `log(Y)` to avoid boundary artifacts at zero.
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
#'   `ggplot2` objects. Otherwise returns the requested `ggplot2` object.
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
  ppc_data <- bayes_ppc_data(object = object, ndraws = ndraws, seed = seed)
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
      x = paste0("Indicator for Y = atom (", object$atom, ")"),
      y = "Proportion"
    )

  continuous_plot <- bayesplot::ppc_dens_overlay_grouped(
    y = continuous_plot_inputs$y,
    yrep = continuous_plot_inputs$yrep,
    group = ppc_data$group_cont
  ) +
    ggplot2::labs(
      title = "Posterior predictive check for the continuous part",
      x = continuous_plot_inputs$x_label,
      y = "Density"
    )

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

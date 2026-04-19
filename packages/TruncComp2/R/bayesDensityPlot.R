utils::globalVariables(
  c(
    "arm_label", "x", "density_mean", "conf.low", "conf.high",
    "atom", "x_label", "y_label", "label", "y_start", "y_end", "hjust"
  )
)

bayes_fit_continuous_support <- function(object) {
  bayes_continuous_support(object$settings$continuous_support)
}

bayes_positive_grid_lower <- function(object, x_limits) {
  observed <- object$data$Y[object$data$A == 1]
  min_positive <- min(observed[observed > 0], na.rm = TRUE)

  if(!is.finite(min_positive) || min_positive <= 0) {
    min_positive <- 1e-4
  }

  lower <- if(x_limits[1] > 0) x_limits[1] else min_positive * 1e-3
  max(lower, .Machine$double.xmin^0.25)
}

bayes_density_grid <- function(object,
                               x = NULL,
                               n = 200L,
                               continuous_support = bayes_fit_continuous_support(object),
                               x_limits = NULL) {
  if(!is.null(x)) {
    if(!(is.numeric(x) && length(x) >= 2L && all(is.finite(x)))) {
      stop("x must be NULL or a finite numeric grid with length at least 2.")
    }

    return(as.numeric(x))
  }

  n <- validate_bayes_positive_integer(n, "n", min_value = 2L)
  observed <- object$data$Y[object$data$A == 1]

  if(!(length(observed) > 0L && all(is.finite(observed)))) {
    stop("Observed non-atom outcomes are required to build the default grid.")
  }

  continuous_support <- bayes_continuous_support(continuous_support)

  if(identical(continuous_support, "bounded_continuous")) {
    score_min <- object$settings$score_min
    score_max <- object$settings$score_max
    if(!(length(score_min) == 1L && length(score_max) == 1L &&
         is.finite(score_min) && is.finite(score_max) && score_min < score_max)) {
      stop("Bounded Bayesian fit settings are missing score_min and score_max.")
    }
    eps <- max(.Machine$double.eps^0.25, (score_max - score_min) * 1e-6)
    return(seq(score_min + eps, score_max - eps, length.out = n))
  }

  if(identical(continuous_support, "bounded_score")) {
    score_values <- object$settings$score_values
    if(!(is.numeric(score_values) && length(score_values) >= 1L && all(is.finite(score_values)))) {
      stop("Bounded-score Bayesian fit settings are missing score_values.")
    }
    return(as.numeric(score_values))
  }

  if(identical(continuous_support, "positive_real") && !is.null(x_limits)) {
    lower <- bayes_positive_grid_lower(object, x_limits)
    upper <- x_limits[2]

    if(!(is.finite(upper) && upper > lower)) {
      upper <- max(observed)
    }

    return(seq(lower, upper, length.out = n))
  }

  bounds <- range(observed)
  span <- diff(bounds)
  scale_hint <- object$settings$y_scale
  if(!(length(scale_hint) == 1L && is.finite(scale_hint) && scale_hint > 0)) {
    scale_hint <- 1
  }

  pad <- if(is.finite(span) && span > 0) {
    0.08 * span
  } else {
    0.5 * scale_hint
  }

  seq(bounds[1] - pad, bounds[2] + pad, length.out = n)
}

bayes_density_default_limits_real_line <- function(object, means, sds, weights) {
  observed <- object$data$Y[object$data$A == 1]
  observed_bounds <- range(c(observed, object$atom))
  scale_hint <- object$settings$y_scale
  if(!(length(scale_hint) == 1L && is.finite(scale_hint) && scale_hint > 0)) {
    scale_hint <- 1
  }

  component_threshold <- max(0.01, 0.1 / dim(weights)[3])
  lower_draw <- numeric(dim(weights)[1] * dim(weights)[2])
  upper_draw <- numeric(dim(weights)[1] * dim(weights)[2])
  index <- 1L

  for(draw in seq_len(dim(weights)[1])) {
    for(arm in seq_len(dim(weights)[2])) {
      keep <- weights[draw, arm, ] >= component_threshold
      if(!any(keep)) {
        keep[which.max(weights[draw, arm, ])] <- TRUE
      }

      lower_draw[index] <- min(means[draw, arm, keep] - 4 * sds[draw, arm, keep])
      upper_draw[index] <- max(means[draw, arm, keep] + 4 * sds[draw, arm, keep])
      index <- index + 1L
    }
  }

  lower <- min(
    observed_bounds[1],
    stats::quantile(lower_draw, probs = 0.01, names = FALSE, na.rm = TRUE)
  )
  upper <- max(
    observed_bounds[2],
    stats::quantile(upper_draw, probs = 0.99, names = FALSE, na.rm = TRUE)
  )

  span <- upper - lower
  pad <- if(is.finite(span) && span > 0) {
    0.04 * span
  } else {
    0.5 * scale_hint
  }

  c(lower - pad, upper + pad)
}

bayes_density_default_limits_positive_real <- function(object, means, shapes, weights) {
  observed <- object$data$Y[object$data$A == 1]
  observed_bounds <- range(c(observed, object$atom))
  scale_hint <- object$settings$y_scale
  if(!(length(scale_hint) == 1L && is.finite(scale_hint) && scale_hint > 0)) {
    scale_hint <- 1
  }

  component_threshold <- max(0.01, 0.1 / dim(weights)[3])
  lower_draw <- numeric(dim(weights)[1] * dim(weights)[2])
  upper_draw <- numeric(dim(weights)[1] * dim(weights)[2])
  index <- 1L

  for(draw in seq_len(dim(weights)[1])) {
    for(arm in seq_len(dim(weights)[2])) {
      keep <- weights[draw, arm, ] >= component_threshold
      if(!any(keep)) {
        keep[which.max(weights[draw, arm, ])] <- TRUE
      }

      lower_draw[index] <- min(stats::qgamma(
        0.001,
        shape = shapes[draw, arm, keep],
        rate = shapes[draw, arm, keep] / means[draw, arm, keep]
      ))
      upper_draw[index] <- max(stats::qgamma(
        0.999,
        shape = shapes[draw, arm, keep],
        rate = shapes[draw, arm, keep] / means[draw, arm, keep]
      ))
      index <- index + 1L
    }
  }

  lower <- min(
    observed[observed > 0],
    stats::quantile(lower_draw, probs = 0.01, names = FALSE, na.rm = TRUE)
  )
  upper <- max(
    observed_bounds[2],
    stats::quantile(upper_draw, probs = 0.99, names = FALSE, na.rm = TRUE)
  )

  span <- upper - lower
  pad <- if(is.finite(span) && span > 0) {
    0.04 * span
  } else {
    0.5 * scale_hint
  }

  c(min(object$atom, max(0, lower - pad)), upper + pad)
}

bayes_density_default_limits_bounded <- function(object) {
  score_min <- object$settings$score_min
  score_max <- object$settings$score_max

  if(!(length(score_min) == 1L &&
       length(score_max) == 1L &&
       is.finite(score_min) &&
       is.finite(score_max) &&
       score_min < score_max)) {
    stop("Bounded Bayesian fit settings are missing score_min and score_max.")
  }

  c(score_min, score_max)
}

bayes_density_component_draws <- function(object) {
  if(!inherits(object, "trunc_comp_bayes_fit")) {
    stop("object must be a trunc_comp_bayes_fit returned by trunc_comp_bayes().")
  }

  if(!isTRUE(object$success)) {
    stop("Estimation failed. Cannot plot posterior densities.")
  }

  if(is.null(object$fit)) {
    stop("Raw Stan mixture parameters are not available.")
  }

  continuous_support <- bayes_fit_continuous_support(object)
  pars <- switch(
    continuous_support,
    real_line = c("w", "mu_comp", "sigma_comp", "rho"),
    positive_real = c("w", "mean_comp", "shape_comp", "rho"),
    bounded_continuous = c("w", "m_comp", "phi_comp", "rho"),
    bounded_score = c("w", "m_comp", "phi_comp", "eta", "rho")
  )

  extracted <- tryCatch(
    rstan::extract(
      object$fit,
      pars = pars,
      permuted = TRUE,
      inc_warmup = FALSE
    ),
    error = identity
  )

  if(inherits(extracted, "error")) {
    stop("Raw Stan mixture parameters are not available.")
  }

  required <- pars
  missing <- required[!vapply(required, function(name) !is.null(extracted[[name]]), logical(1))]
  if(length(missing) > 0L) {
    stop(
      "Raw Stan mixture parameters are missing required variables: ",
      paste(missing, collapse = ", "),
      "."
    )
  }

  center <- object$settings$y_center
  scale <- object$settings$y_scale
  if(!(length(scale) == 1L && is.finite(scale) && scale > 0)) {
    stop("Bayesian fit settings are missing the outcome scale needed for density plotting.")
  }

  if(identical(continuous_support, "real_line")) {
    if(!(length(center) == 1L && is.finite(center))) {
      stop("Bayesian fit settings are missing the outcome center needed for density plotting.")
    }

    return(list(
      support = continuous_support,
      weights = extracted$w,
      rho = extracted$rho,
      observed_prob = 1 - extracted$rho,
      means = center + scale * extracted$mu_comp,
      sds = scale * extracted$sigma_comp
    ))
  }

  if(identical(continuous_support, "bounded_continuous")) {
    return(list(
      support = continuous_support,
      weights = extracted$w,
      rho = extracted$rho,
      observed_prob = 1 - extracted$rho,
      m_comp = extracted$m_comp,
      phi_comp = extracted$phi_comp,
      score_min = object$settings$score_min,
      score_max = object$settings$score_max,
      score_range = object$settings$score_range
    ))
  }

  if(identical(continuous_support, "bounded_score")) {
    return(list(
      support = continuous_support,
      weights = extracted$w,
      rho = extracted$rho,
      observed_prob = 1 - extracted$rho,
      m_comp = extracted$m_comp,
      phi_comp = extracted$phi_comp,
      eta = extracted$eta,
      score_min = object$settings$score_min,
      score_max = object$settings$score_max,
      score_range = object$settings$score_range,
      score_values = object$settings$score_values,
      heaping_grids = object$settings$heaping_grids,
      heaping = object$settings$heaping,
      eta_groups = object$settings$eta_groups,
      eta_group_by_arm = object$settings$eta_group_by_arm,
      bin_lower = object$settings$bin_lower,
      bin_upper = object$settings$bin_upper,
      bin_valid = object$settings$bin_valid
    ))
  }

  list(
    support = continuous_support,
    weights = extracted$w,
    rho = extracted$rho,
    observed_prob = 1 - extracted$rho,
    means = scale * extracted$mean_comp,
    shapes = extracted$shape_comp
  )
}

bayes_positive_kernel_density <- function(x, mean, shape) {
  x_eval <- pmax(x, .Machine$double.xmin)
  out <- stats::dgamma(x_eval, shape = shape, rate = shape / mean)
  out[x < 0] <- 0
  out
}

bayes_bounded_kernel_density <- function(x, m, phi, score_min, score_max) {
  score_range <- score_max - score_min
  x_unit <- (x - score_min) / score_range
  in_bounds <- x_unit >= 0 & x_unit <= 1
  x_eval <- pmin(pmax(x_unit, .Machine$double.eps^0.25), 1 - .Machine$double.eps^0.25)
  out <- stats::dbeta(x_eval, shape1 = m * phi, shape2 = (1 - m) * phi) / score_range
  out[!in_bounds] <- 0
  out
}

bayes_observed_density <- function(x, weights, means = NULL, pi,
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
        return(stats::dnorm(x, mean = means[h], sd = sds[h]))
      }

      if(identical(continuous_support, "bounded_continuous")) {
        return(bayes_bounded_kernel_density(
          x = x,
          m = m_comp[h],
          phi = phi_comp[h],
          score_min = score_min,
          score_max = score_max
        ))
      }

      bayes_positive_kernel_density(x, mean = means[h], shape = shapes[h])
    },
    numeric(length(x))
  )

  as.numeric(pi * (kernels %*% weights))
}

bayes_score_density_plot_data <- function(object, extracted, conf.level) {
  score_values <- extracted$score_values
  n_draws <- dim(extracted$weights)[1]
  arm_labels <- c("Control", "Treatment")
  density_rows <- vector("list", 2)
  max_upper <- numeric(2)

  for(r in 1:2) {
    mass_matrix <- t(vapply(
      seq_len(n_draws),
      function(draw) {
        eta_group <- extracted$eta_group_by_arm[[r]]
        pmf <- bayes_score_pmf(
          weights = extracted$weights[draw, r, ],
          m_comp = extracted$m_comp[draw, r, ],
          phi_comp = extracted$phi_comp[draw, r, ],
          eta = extracted$eta[draw, eta_group, ],
          bin_lower = extracted$bin_lower,
          bin_upper = extracted$bin_upper,
          bin_valid = extracted$bin_valid
        )
        extracted$observed_prob[draw, r] * pmf
      },
      numeric(length(score_values))
    ))

    intervals <- t(vapply(
      seq_along(score_values),
      function(index) bayes_equal_tail_interval(mass_matrix[, index], conf.level),
      numeric(2)
    ))

    density_rows[[r]] <- data.frame(
      arm = r - 1L,
      arm_label = factor(arm_labels[r], levels = arm_labels),
      x = score_values,
      density_mean = colMeans(mass_matrix),
      conf.low = intervals[, 1],
      conf.high = intervals[, 2],
      row.names = NULL
    )

    max_upper[r] <- max(intervals[, 2], na.rm = TRUE)
  }

  density_data <- do.call(rbind, density_rows)
  x_limits <- range(c(score_values, object$atom))
  x_span <- diff(x_limits)
  if(!(is.finite(x_span) && x_span > 0)) {
    x_span <- 1
  }
  x_limits <- x_limits + c(-0.03, 0.03) * x_span

  atom_mass <- colMeans(extracted$rho)
  spike_height <- pmax(atom_mass, 1.05 * max_upper, 0.05)
  panel_upper <- spike_height + 0.14 * pmax(spike_height, 0.05)
  midpoint <- mean(x_limits)
  label_offset <- 0.035 * diff(x_limits)
  place_right <- object$atom <= midpoint
  atom_data <- data.frame(
    arm = 0:1,
    arm_label = factor(arm_labels, levels = arm_labels),
    atom = rep(object$atom, 2),
    atom_mass = atom_mass,
    y_start = rep(0, 2),
    y_end = spike_height,
    x_label = ifelse(place_right, object$atom + label_offset, object$atom - label_offset),
    y_label = pmin(spike_height + 0.06 * panel_upper, 0.96 * panel_upper),
    hjust = ifelse(place_right, 0, 1),
    label = sprintf("Atom mass = %.2f", atom_mass),
    row.names = NULL
  )

  list(
    plot_type = "mass",
    density_data = density_data,
    atom_data = atom_data,
    x = score_values,
    x_limits = x_limits,
    y_limits = c(0, max(panel_upper, na.rm = TRUE)),
    conf.level = conf.level
  )
}

bayes_density_plot_data <- function(object,
                                    x = NULL,
                                    n = 200L,
                                    conf.level = object$conf.level) {
  if(!inherits(object, "trunc_comp_bayes_fit")) {
    stop("object must be a trunc_comp_bayes_fit returned by trunc_comp_bayes().")
  }

  if(!isTRUE(object$success)) {
    stop("Estimation failed. Cannot plot posterior densities.")
  }

  conf.level <- validateConfidenceLevel(conf.level)
  extracted <- bayes_density_component_draws(object)
  continuous_support <- extracted$support

  if(identical(continuous_support, "bounded_score")) {
    if(!is.null(x)) {
      stop("x is not used for bounded-score posterior mass plots.")
    }
    validate_bayes_positive_integer(n, "n", min_value = 2L)
    return(bayes_score_density_plot_data(
      object = object,
      extracted = extracted,
      conf.level = conf.level
    ))
  }

  weights <- extracted$weights
  means <- extracted$means
  observed_prob <- extracted$observed_prob
  n_grid <- if(is.null(x)) validate_bayes_positive_integer(n, "n", min_value = 2L) else NULL

  if(is.null(x)) {
    x_limits <- if(identical(continuous_support, "real_line")) {
      bayes_density_default_limits_real_line(
        object = object,
        means = means,
        sds = extracted$sds,
        weights = weights
      )
    } else if(identical(continuous_support, "positive_real")) {
      bayes_density_default_limits_positive_real(
        object = object,
        means = means,
        shapes = extracted$shapes,
        weights = weights
      )
    } else {
      bayes_density_default_limits_bounded(object)
    }

    grid <- bayes_density_grid(
      object = object,
      x = NULL,
      n = n_grid,
      continuous_support = continuous_support,
      x_limits = x_limits
    )
  } else {
    grid <- bayes_density_grid(
      object,
      x = x,
      n = n,
      continuous_support = continuous_support
    )
    x_limits <- range(c(grid, object$atom))
  }

  density_grid <- grid
  if(identical(continuous_support, "positive_real")) {
    density_grid <- grid[grid > 0]

    if(length(density_grid) < 2L) {
      stop(
        "For positive-support Bayesian density plots, x must contain at least two positive grid points."
      )
    }
  }
  if(identical(continuous_support, "bounded_continuous")) {
    density_grid <- grid[
      grid >= extracted$score_min &
        grid <= extracted$score_max
    ]

    if(length(density_grid) < 2L) {
      stop(
        "For bounded-continuous Bayesian density plots, x must contain at least two grid points inside [score_min, score_max]."
      )
    }
    x_limits <- bayes_density_default_limits_bounded(object)
  }

  n_draws <- dim(weights)[1]
  arm_labels <- c("Control", "Treatment")
  density_rows <- vector("list", 2)
  max_upper <- numeric(2)

  for(r in 1:2) {
    density_matrix <- t(vapply(
      seq_len(n_draws),
      function(draw) {
        bayes_observed_density(
          x = density_grid,
          weights = weights[draw, r, ],
          means = if(identical(continuous_support, "bounded_continuous")) NULL else means[draw, r, ],
          pi = observed_prob[draw, r],
          continuous_support = continuous_support,
          sds = if(identical(continuous_support, "real_line")) extracted$sds[draw, r, ] else NULL,
          shapes = if(identical(continuous_support, "positive_real")) extracted$shapes[draw, r, ] else NULL,
          m_comp = if(identical(continuous_support, "bounded_continuous")) extracted$m_comp[draw, r, ] else NULL,
          phi_comp = if(identical(continuous_support, "bounded_continuous")) extracted$phi_comp[draw, r, ] else NULL,
          score_min = if(identical(continuous_support, "bounded_continuous")) extracted$score_min else NULL,
          score_max = if(identical(continuous_support, "bounded_continuous")) extracted$score_max else NULL
        )
      },
      numeric(length(density_grid))
    ))

    intervals <- t(vapply(
      seq_along(density_grid),
      function(index) bayes_equal_tail_interval(density_matrix[, index], conf.level),
      numeric(2)
    ))

    density_rows[[r]] <- data.frame(
      arm = r - 1L,
      arm_label = factor(arm_labels[r], levels = arm_labels),
      x = density_grid,
      density_mean = colMeans(density_matrix),
      conf.low = intervals[, 1],
      conf.high = intervals[, 2],
      row.names = NULL
    )

    max_upper[r] <- max(intervals[, 2], na.rm = TRUE)
  }

  density_data <- do.call(rbind, density_rows)
  x_span <- diff(x_limits)
  if(!(is.finite(x_span) && x_span > 0)) {
    x_span <- 1
  }

  atom_mass <- colMeans(extracted$rho)
  spike_height <- pmax(atom_mass, 1.05 * max_upper, 0.05)
  panel_upper <- spike_height + 0.14 * pmax(spike_height, 0.05)
  midpoint <- mean(x_limits)
  label_offset <- 0.035 * x_span
  place_right <- object$atom <= midpoint
  atom_data <- data.frame(
    arm = 0:1,
    arm_label = factor(arm_labels, levels = arm_labels),
    atom = rep(object$atom, 2),
    atom_mass = atom_mass,
    y_start = rep(0, 2),
    y_end = spike_height,
    x_label = ifelse(place_right, object$atom + label_offset, object$atom - label_offset),
    y_label = pmin(spike_height + 0.06 * panel_upper, 0.96 * panel_upper),
    hjust = ifelse(place_right, 0, 1),
    label = sprintf("Atom mass = %.2f", atom_mass),
    row.names = NULL
  )

  list(
    plot_type = "density",
    density_data = density_data,
    atom_data = atom_data,
    x = grid,
    x_limits = x_limits,
    y_limits = c(0, max(panel_upper, na.rm = TRUE)),
    conf.level = conf.level
  )
}

#' Plot posterior outcome densities or masses from a Bayesian fit
#'
#' Visualizes the arm-specific posterior mean continuous densities implied by a
#' fitted `"trunc_comp_bayes_fit"` object, or posterior mean reported-score
#' masses for `continuous_support = "bounded_score"`. The plotted survivor
#' distribution is weighted by the posterior probability of being observed away
#' from the atom. The atom itself is shown separately as a solid vertical spike
#' with an arrowhead, annotated with the posterior mean atom probability in each
#' treatment arm. Fits created with `continuous_support = "real_line"` use
#' Gaussian kernels, `continuous_support = "positive_real"` uses Gamma kernels,
#' `continuous_support = "bounded_continuous"` uses Beta kernels with the
#' outcome-scale Jacobian, and `continuous_support = "bounded_score"` uses a
#' discrete posterior predictive mass plot.
#'
#' @param object A successful `"trunc_comp_bayes_fit"` object returned by
#'   [trunc_comp_bayes()].
#' @param x Optional numeric grid on the outcome scale where the posterior
#'   densities should be evaluated. When omitted, a default grid is constructed
#'   from the pooled non-atom outcomes and widened using the fitted posterior
#'   mixture support so the density tails are visible.
#' @param n Number of grid points used when `x` is omitted. Must be an integer
#'   greater than or equal to `2`.
#' @param conf.level Pointwise credible level for the density ribbon.
#' @return A `ggplot2` plot showing the arm-specific posterior mean densities
#'   or score masses, pointwise credible intervals, and a solid atom spike with
#'   its posterior mean atom mass.
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
#' posterior_density_plot(fit)
#' }
#' @export
posterior_density_plot <- function(object,
                                   x = NULL,
                                   n = 200L,
                                   conf.level = object$conf.level) {
  plot_data <- bayes_density_plot_data(
    object = object,
    x = x,
    n = n,
    conf.level = conf.level
  )

  arm_colors <- c(
    Control = "#1b9e77",
    Treatment = "#d95f02"
  )

  if(identical(plot_data$plot_type, "mass")) {
    return(ggplot2::ggplot() +
      ggplot2::geom_linerange(
        data = plot_data$density_data,
        ggplot2::aes(x = x, ymin = conf.low, ymax = conf.high, color = arm_label),
        linewidth = 0.45,
        show.legend = FALSE
      ) +
      ggplot2::geom_point(
        data = plot_data$density_data,
        ggplot2::aes(x = x, y = density_mean, color = arm_label),
        size = 1.5,
        show.legend = FALSE
      ) +
      ggplot2::geom_segment(
        data = plot_data$atom_data,
        ggplot2::aes(x = atom, xend = atom, y = y_start, yend = y_end),
        inherit.aes = FALSE,
        linewidth = 0.7,
        lineend = "round",
        arrow = grid::arrow(
          type = "closed",
          angle = 20,
          length = grid::unit(0.11, "inches")
        ),
        color = "grey30"
      ) +
      ggplot2::geom_text(
        data = plot_data$atom_data,
        ggplot2::aes(x = x_label, y = y_label, label = label, hjust = hjust),
        inherit.aes = FALSE,
        vjust = 0,
        size = 3.3
      ) +
      ggplot2::facet_wrap(ggplot2::vars(arm_label), scales = "fixed") +
      ggplot2::labs(
        x = "Reported score",
        y = "Posterior predictive probability"
      ) +
      ggplot2::scale_color_manual(values = arm_colors) +
      ggplot2::coord_cartesian(xlim = plot_data$x_limits, ylim = plot_data$y_limits) +
      ggplot2::theme_minimal())
  }

  ggplot2::ggplot() +
    ggplot2::geom_ribbon(
      data = plot_data$density_data,
      ggplot2::aes(x = x, y = density_mean, ymin = conf.low, ymax = conf.high),
      alpha = 0.25,
      linewidth = 0,
      fill = "darkgray",
      show.legend = FALSE
    ) +
    ggplot2::geom_line(
      data = plot_data$density_data,
      ggplot2::aes(x = x, y = density_mean, color = arm_label),
      linewidth = 0.9,
      show.legend = FALSE
    ) +
    ggplot2::geom_segment(
      data = plot_data$atom_data,
      ggplot2::aes(x = atom, xend = atom, y = y_start, yend = y_end),
      inherit.aes = FALSE,
      linewidth = 0.7,
      lineend = "round",
      arrow = grid::arrow(
        type = "closed",
        angle = 20,
        length = grid::unit(0.11, "inches")
      ),
      color = "grey30"
    ) +
    ggplot2::geom_text(
      data = plot_data$atom_data,
      ggplot2::aes(x = x_label, y = y_label, label = label, hjust = hjust),
      inherit.aes = FALSE,
      vjust = 0,
      size = 3.3
    ) +
    ggplot2::facet_wrap(ggplot2::vars(arm_label), scales = "fixed") +
    ggplot2::labs(
      x = "Outcome",
      y = "Posterior predictive density"
    ) +
    ggplot2::scale_color_manual(values = arm_colors) +
    ggplot2::coord_cartesian(xlim = plot_data$x_limits, ylim = plot_data$y_limits) +
    ggplot2::theme_minimal()
}

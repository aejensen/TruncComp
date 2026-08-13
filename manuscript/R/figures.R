.surface_component <- function(surface, name) {
  if (!is.null(surface[[name]])) {
    return(surface[[name]])
  }
  stop(sprintf("Joint likelihood surface is missing component `%s`.", name), call. = FALSE)
}

.placeholder_plot <- function(message) {
  graphics::plot.new()
  graphics::text(0.5, 0.5, message, cex = 1)
}

.figure_title_size <- 12
.base_figure_title_cex <- 1

build_application_figure <- function(application_results, path) {
  d <- application_results$data
  metadata <- application_results$metadata
  s_ci <- application_results$surface
  m_splrt <- application_results$model_splrt
  group_labels <- metadata$group_labels
  mu_axis <- .surface_component(s_ci, "mu_delta")
  log_or_axis <- .surface_component(s_ci, "log_or_delta")
  surface <- .surface_component(s_ci, "surface")
  non_atom <- d[d$A == 1L, , drop = FALSE]

  death_alive <- rbind(
    "Death by 6 months" = tapply(d$A == 0L, factor(d$R, levels = 0:1), sum),
    "Alive with observed EQ-VAS" = tapply(d$A == 1L, factor(d$R, levels = 0:1), sum)
  )
  colnames(death_alive) <- group_labels
  endpoint_y_top <- ceiling(max(colSums(death_alive, na.rm = TRUE)) / 200) * 200
  endpoint_y_ticks <- seq(0, endpoint_y_top, by = 200)

  save_pdf_plot(path, width = 8.4, height = 3.2, expr = {
    op <- graphics::par(
      mfrow = c(1, 3),
      mgp = c(2.2, 0.8, 0),
      mar = c(3.4, 3.5, 2.2, 0.8),
      font.main = 2,
      cex.main = .base_figure_title_cex
    )
    on.exit(graphics::par(op), add = TRUE)

    graphics::barplot(
      death_alive,
      beside = FALSE,
      col = c("gray35", "gray78"),
      border = "white",
      ylim = c(0, endpoint_y_top),
      axes = FALSE,
      ylab = "Participants",
      main = "Endpoint status counts",
      cex.names = 0.9,
      legend.text = rownames(death_alive),
      args.legend = list(x = "topright", bty = "n", cex = 0.66)
    )
    graphics::axis(2, at = endpoint_y_ticks, cex.axis = 1)

    hist0 <- graphics::hist(non_atom$Y[non_atom$R == 0L], breaks = seq(0, 100, by = 10), plot = FALSE)
    hist1 <- graphics::hist(non_atom$Y[non_atom$R == 1L], breaks = seq(0, 100, by = 10), plot = FALSE)
    y_max <- max(hist0$density, hist1$density)
    graphics::plot(
      hist0,
      freq = FALSE,
      col = grDevices::rgb(0.35, 0.35, 0.35, 0.35),
      border = "white",
      xlim = c(0, 100),
      ylim = c(0, y_max * 1.15),
      main = "EQ-VAS among alive participants",
      xlab = "6-month EQ-VAS",
      ylab = "Density"
    )
    graphics::plot(
      hist1,
      freq = FALSE,
      col = grDevices::rgb(0.05, 0.37, 0.55, 0.35),
      border = "white",
      add = TRUE
    )
    graphics::lines(stats::density(non_atom$Y[non_atom$R == 0L], from = 0, to = 100), col = "gray20", lwd = 1.6)
    graphics::lines(stats::density(non_atom$Y[non_atom$R == 1L], from = 0, to = 100), col = "#0b5d7a", lwd = 1.6)
    graphics::legend(
      "topleft",
      legend = group_labels,
      col = c("gray20", "#0b5d7a"),
      lwd = 1.6,
      bty = "n",
      cex = 0.8
    )

    graphics::image(
      mu_axis,
      log_or_axis,
      surface,
      col = rev(grDevices::hcl.colors(256, "YlOrRd", rev = FALSE)),
      useRaster = FALSE,
      xlab = expression(mu[delta]),
      ylab = expression(log(alpha[delta])),
      main = "Semiparametric LR region"
    )
    graphics::points(m_splrt$mu_delta, log(m_splrt$alpha_delta), pch = 19, cex = 1.1)
    graphics::abline(v = 0, h = 0, lty = 2, col = "gray30")
    for (level in c(0.5, 0.8, 0.95, 0.99)) {
      graphics::contour(
        mu_axis,
        log_or_axis,
        surface,
        add = TRUE,
        levels = stats::qchisq(level, 2),
        labels = level,
        drawlabels = TRUE,
        lwd = 0.8
      )
    }
  })
}

build_application_posterior_figure <- function(application_results, path) {
  bayes <- application_results$bayes

  save_pdf_plot(path, width = 7.2, height = 2.5, expr = {
    op <- graphics::par(
      mfrow = c(1, 3),
      mgp = c(2.0, 0.8, 0),
      mar = c(4.8, 3.1, 2.1, 0.7),
      font.main = 2,
      cex.main = .base_figure_title_cex
    )
    on.exit(graphics::par(op), add = TRUE)

    if (!isTRUE(bayes$success)) {
      .placeholder_plot(paste("Bayesian fit unavailable:", bayes$error))
      return(invisible(NULL))
    }

    draws <- bayes$fit$draws
    panels <- list(
      list(x = draws$mu_delta, ref = 0, main = expression(mu[delta]), xlab = "Mean EQ-VAS difference\namong alive with observed score", xlab_line = 2.8),
      list(x = log(draws$alpha_delta), ref = 0, main = expression(log(alpha[delta])), xlab = "Log odds ratio for being alive\nwith observed EQ-VAS", xlab_line = 2.8),
      list(x = draws$delta, ref = 0, main = "Coded mean contrast", xlab = "Atom coded as -1", xlab_line = 2.5)
    )

    for (panel in panels) {
      den <- stats::density(panel$x)
      graphics::plot(
          den,
          main = panel$main,
          xlab = "",
          ylab = "Posterior density",
          col = "black",
          lwd = 1.6,
          bty = "l"
        )
        graphics::abline(v = panel$ref, lty = 2, col = "black")
        graphics::title(xlab = panel$xlab, line = panel$xlab_line)
      }
  })
}

build_application_ppc_figure <- function(application_results, path) {
  bayes <- application_results$bayes
  metadata <- application_results$metadata

  if (!isTRUE(bayes$success) || !requireNamespace("gridExtra", quietly = TRUE)) {
    save_pdf_plot(path, width = 7, height = 3, expr = {
      .placeholder_plot(if (isTRUE(bayes$success)) "gridExtra is unavailable." else paste("PPC unavailable:", bayes$error))
    })
    return(invisible(path))
  }

  ppc_theme <- ggplot2::theme_classic() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = .figure_title_size),
      plot.subtitle = ggplot2::element_text(hjust = 0.5),
      strip.background = ggplot2::element_blank(),
      legend.key = ggplot2::element_blank()
    )
  old_bayesplot_theme <- bayesplot::bayesplot_theme_set(ppc_theme)
  on.exit(bayesplot::bayesplot_theme_set(old_bayesplot_theme), add = TRUE)
  group_labels <- c(Control = metadata$group_labels[[1]], Treatment = metadata$group_labels[[2]])

  plots <- tryCatch(
    TruncComp::posterior_predictive_check(bayes$fit, type = "both", ndraws = 200, seed = bayes$settings$seed + 1L),
    error = function(e) e
  )

  if (inherits(plots, "error")) {
    save_pdf_plot(path, width = 7, height = 3, expr = {
      .placeholder_plot(paste("PPC plot failed:", conditionMessage(plots)))
    })
    return(invisible(path))
  }
  atom_data <- plots$atom$data
  atom_data$group <- factor(group_labels[as.character(atom_data$group)], levels = unname(group_labels))
  atom_plot <- ggplot2::ggplot(atom_data, ggplot2::aes(x = factor(x))) +
    ggplot2::geom_col(
      ggplot2::aes(y = y_obs),
      fill = "grey82",
      width = 0.75
    ) +
    ggplot2::geom_pointrange(
      ggplot2::aes(y = m, ymin = l, ymax = h),
      color = "black",
      linewidth = 0.45,
      size = 0.45
    ) +
    ggplot2::facet_wrap(ggplot2::vars(group), scales = "fixed") +
    ggplot2::labs(
      title = "Posterior predictive check for death",
      subtitle = sub("Posterior predictive p =", "Posterior predictive check, p =", plots$atom$labels$subtitle, fixed = TRUE),
      x = "Death indicator",
      y = plots$atom$labels$y
    ) +
    ppc_theme

  continuous_data <- plots$continuous$data
  continuous_data$arm_label <- factor(group_labels[as.character(continuous_data$arm_label)], levels = unname(group_labels))
  continuous_plot <- ggplot2::ggplot(continuous_data, ggplot2::aes(x = x)) +
    ggplot2::geom_col(
      ggplot2::aes(y = observed),
      fill = "grey82",
      width = 0.8
    ) +
    ggplot2::geom_linerange(
      ggplot2::aes(ymin = conf.low, ymax = conf.high),
      color = "black",
      linewidth = 0.35
    ) +
    ggplot2::geom_point(
      ggplot2::aes(y = replicated_mean),
      color = "black",
      size = 1.1
    ) +
    ggplot2::facet_wrap(ggplot2::vars(arm_label), scales = "fixed") +
    ggplot2::labs(
      title = "Posterior predictive check for non-atom scores",
      subtitle = sub("Posterior predictive p =", "Posterior predictive check, p =", plots$continuous$labels$subtitle, fixed = TRUE),
      x = plots$continuous$labels$x,
      y = plots$continuous$labels$y
    ) +
    ppc_theme

  grDevices::pdf(path, width = 9.2, height = 3.6, useDingbats = FALSE)
  on.exit(grDevices::dev.off(), add = TRUE)
  gridExtra::grid.arrange(atom_plot, continuous_plot, ncol = 2)
  invisible(path)
}

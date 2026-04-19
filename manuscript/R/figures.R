.surface_component <- function(surface, name, legacy_name = NULL) {
  if (!is.null(surface[[name]])) {
    return(surface[[name]])
  }
  if (!is.null(legacy_name) && !is.null(surface[[legacy_name]])) {
    return(surface[[legacy_name]])
  }
  stop(sprintf("Joint likelihood surface is missing component `%s`.", name), call. = FALSE)
}

.placeholder_plot <- function(message) {
  graphics::plot.new()
  graphics::text(0.5, 0.5, message, cex = 1)
}

build_application_figure <- function(application_results, path) {
  d <- application_results$data
  metadata <- application_results$metadata
  s_ci <- application_results$surface
  m_splrt <- application_results$model_splrt
  group_labels <- metadata$group_labels
  mu_axis <- .surface_component(s_ci, "mu_delta", "muDelta")
  log_or_axis <- .surface_component(s_ci, "log_or_delta", "logORdelta")
  surface <- .surface_component(s_ci, "surface")
  non_atom <- d[d$A == 1L, , drop = FALSE]

  death_alive <- rbind(
    "Death" = tapply(d$A == 0L, factor(d$R, levels = 0:1), sum),
    "Alive with EQ-VAS" = tapply(d$A == 1L, factor(d$R, levels = 0:1), sum)
  )
  colnames(death_alive) <- group_labels

  save_pdf_plot(path, width = 8.4, height = 3.2, expr = {
    op <- graphics::par(mfrow = c(1, 3), mgp = c(2.2, 0.8, 0), mar = c(3.4, 3.5, 2.2, 0.8))
    on.exit(graphics::par(op), add = TRUE)

    graphics::barplot(
      death_alive,
      beside = FALSE,
      col = c("gray35", "gray78"),
      border = "white",
      ylab = "Participants",
      main = "Endpoint components",
      cex.axis = 1,
      cex.names = 0.9,
      legend.text = rownames(death_alive),
      args.legend = list(x = "topright", bty = "n", cex = 0.75)
    )

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
      main = "Survivor EQ-VAS",
      xlab = "EQ-VAS at 6 months",
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
      main = "SPLRT region"
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
    op <- graphics::par(mfrow = c(1, 3), mgp = c(2.0, 0.8, 0), mar = c(3.2, 3.1, 2.1, 0.7))
    on.exit(graphics::par(op), add = TRUE)

    if (!isTRUE(bayes$success)) {
      .placeholder_plot(paste("Bayesian fit unavailable:", bayes$error))
      return(invisible(NULL))
    }

    draws <- bayes$fit$draws
    panels <- list(
      list(x = draws$mu_delta, ref = 0, main = expression(mu[delta]), xlab = "Mean EQ-VAS difference"),
      list(x = log(draws$alpha_delta), ref = 0, main = expression(log(alpha[delta])), xlab = "Log odds ratio alive"),
      list(x = draws$delta, ref = 0, main = expression(Delta), xlab = "Combined mean contrast")
    )

    for (panel in panels) {
      den <- stats::density(panel$x)
      graphics::plot(
        den,
        main = panel$main,
        xlab = panel$xlab,
        ylab = "Posterior density",
        col = "#0b5d7a",
        lwd = 1.6
      )
      graphics::polygon(den, col = grDevices::rgb(0.05, 0.37, 0.55, 0.20), border = NA)
      graphics::lines(den, col = "#0b5d7a", lwd = 1.6)
      graphics::abline(v = panel$ref, lty = 2, col = "gray25")
    }
  })
}

build_application_ppc_figure <- function(application_results, path) {
  bayes <- application_results$bayes

  if (!isTRUE(bayes$success) || !requireNamespace("gridExtra", quietly = TRUE)) {
    save_pdf_plot(path, width = 7, height = 3, expr = {
      .placeholder_plot(if (isTRUE(bayes$success)) "gridExtra is unavailable." else paste("PPC unavailable:", bayes$error))
    })
    return(invisible(path))
  }

  plots <- tryCatch(
    TruncComp2::posterior_predictive_check(bayes$fit, type = "both", ndraws = 50, seed = bayes$settings$seed + 2L),
    error = function(e) e
  )

  if (inherits(plots, "error")) {
    save_pdf_plot(path, width = 7, height = 3, expr = {
      .placeholder_plot(paste("PPC plot failed:", conditionMessage(plots)))
    })
    return(invisible(path))
  }

  grDevices::pdf(path, width = 8.2, height = 3.6, useDingbats = FALSE)
  on.exit(grDevices::dev.off(), add = TRUE)
  gridExtra::grid.arrange(plots$atom, plots$continuous, ncol = 2)
  invisible(path)
}

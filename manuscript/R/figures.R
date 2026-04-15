build_example_histogram <- function(example_results, path) {
  d <- example_results$data

  save_pdf_plot(path, width = 6, height = 3, expr = {
    op <- graphics::par(mfrow = c(1, 2), mgp = c(0.5, 0.5, 0), mar = c(3.5, 2, 1, 0.5))
    on.exit(graphics::par(op), add = TRUE)

    graphics::hist(
      subset(d, R == 0)$Y,
      xlab = "",
      xlim = c(0, 10),
      main = expression(Y ~ "|" ~ R == 0),
      ylim = c(0, 14),
      col = "gray70",
      border = "gray90",
      cex.axis = 0.7,
      cex.main = 0.8,
      cex.lab = 0.8,
      ylab = ""
    )
    graphics::hist(
      subset(d, R == 1)$Y,
      xlab = "",
      xlim = c(0, 10),
      main = expression(Y ~ "|" ~ R == 1),
      ylim = c(0, 14),
      col = "gray70",
      border = "gray90",
      cex.axis = 0.7,
      cex.main = 0.8,
      cex.lab = 0.8,
      ylab = ""
    )
  })
}

build_example_surface_figure <- function(example_results, path) {
  save_pdf_plot(path, width = 5, height = 4.5, expr = {
    op <- graphics::par(mgp = c(1.6, 0.6, 0), cex.axis = 0.8, cex.lab = 0.8, bty = "n")
    on.exit(graphics::par(op), add = TRUE)

    suppressMessages(stats::confint(
      example_results$model,
      type = "simultaneous",
      plot = TRUE,
      offset = 1.4,
      resolution = 50
    ))
    graphics::abline(v = 0, lty = 2)
    graphics::abline(h = 0, lty = 2)
  })
}

build_power_curves_figure <- function(simulation_results, path) {
  power1 <- set_method_colnames(simulation_results$power1)
  power2 <- set_method_colnames(simulation_results$power2)
  power3 <- set_method_colnames(simulation_results$power3)
  power4 <- set_method_colnames(simulation_results$power4)
  n_seq <- simulation_results$nSeq
  cols <- c("firebrick1", "forestgreen", "cornflowerblue", "black")

  save_pdf_plot(path, width = 7, height = 4.6, expr = {
    op <- graphics::par(mfrow = c(2, 2), bty = "n", mgp = c(1.6, 0.6, 0), mar = c(3, 3, 2, 0))
    on.exit(graphics::par(op), add = TRUE)

    graphics::plot(
      n_seq, power1[, "Wilcoxon"], ylim = c(0, 1), type = "l", col = cols[4],
      xlab = "# observations pr group", ylab = "Power", lwd = 2,
      cex.axis = 0.85, cex.lab = 1.1
    )
    graphics::title("Scenario 1")
    graphics::lines(n_seq, power1[, "T-test"], col = cols[3], lwd = 2)
    graphics::lines(n_seq, power1[, "LRT"], col = cols[1], lwd = 2)
    graphics::lines(n_seq, power1[, "SPLRT"], col = cols[2], lwd = 2)
    graphics::legend(
      "topleft",
      c("Parametric LRT", "Semi-parametric LRT", "t-test", "Wilcoxon"),
      lwd = 2,
      lty = 1,
      col = cols,
      bty = "n",
      cex = 0.75
    )

    graphics::plot(
      n_seq, power2[, "Wilcoxon"], ylim = c(0, 1), type = "l", col = cols[4],
      xlab = "#observations pr group", ylab = "Power", lwd = 2,
      cex.axis = 0.85, cex.lab = 1.1
    )
    graphics::title("Scenario 2")
    graphics::lines(n_seq, power2[, "T-test"], col = cols[3], lwd = 2)
    graphics::lines(n_seq, power2[, "LRT"], col = cols[1], lwd = 2)
    graphics::lines(n_seq, power2[, "SPLRT"], col = cols[2], lwd = 2)

    graphics::plot(
      n_seq, power3[, "Wilcoxon"], ylim = c(0, 1), type = "l", col = cols[4],
      xlab = "#observations pr group", ylab = "Power", lwd = 2,
      cex.axis = 0.85, cex.lab = 1.1
    )
    graphics::title("Scenario 3")
    graphics::lines(n_seq, power3[, "T-test"], col = cols[3], lwd = 2)
    graphics::lines(n_seq, power3[, "LRT"], col = cols[1], lwd = 2)
    graphics::lines(n_seq, power3[, "SPLRT"], col = cols[2], lwd = 2)

    graphics::plot(
      n_seq, power4[, "Wilcoxon"], ylim = c(0, 1), type = "l", col = cols[4],
      xlab = "#observations pr group", ylab = "Power", lwd = 2,
      cex.axis = 0.85, cex.lab = 1.1
    )
    graphics::title("Scenario 4")
    graphics::lines(n_seq, power4[, "T-test"], col = cols[3], lwd = 2)
    graphics::lines(n_seq, power4[, "LRT"], col = cols[1], lwd = 2)
    graphics::lines(n_seq, power4[, "SPLRT"], col = cols[2], lwd = 2)
  })
}

build_application_figure <- function(application_results, path) {
  d <- application_results$data
  metadata <- application_results$metadata
  s_ci <- application_results$surface
  m_splrt <- application_results$model_splrt
  display_y <- if (is.null(metadata$histogram_cap)) d$Y else pmin(d$Y, metadata$histogram_cap)
  display_breaks <- metadata$histogram_breaks
  count_max <- max(
    graphics::hist(display_y[d$R == 0], breaks = display_breaks, plot = FALSE)$counts,
    graphics::hist(display_y[d$R == 1], breaks = display_breaks, plot = FALSE)$counts
  )

  save_pdf_plot(path, width = 8, height = 3, expr = {
    op <- graphics::par(mfrow = c(1, 3), mgp = c(2.2, 1, 0), mar = c(3.2, 3.5, 2, 0))
    on.exit(graphics::par(op), add = TRUE)

    graphics::hist(
      display_y[d$R == 0],
      breaks = display_breaks,
      xlim = metadata$histogram_xlim,
      main = metadata$group_labels[[1]],
      ylim = c(0, count_max),
      col = "gray70",
      border = "gray90",
      cex.axis = 1.2,
      cex.main = 1.5,
      cex.lab = 1.5,
      ylab = "",
      xlab = metadata$display_x_label
    )
    graphics::hist(
      display_y[d$R == 1],
      breaks = display_breaks,
      xlim = metadata$histogram_xlim,
      main = metadata$group_labels[[2]],
      ylim = c(0, count_max),
      col = "gray70",
      border = "gray90",
      cex.axis = 1.2,
      cex.main = 1.5,
      cex.lab = 1.5,
      ylab = "",
      xlab = metadata$display_x_label
    )
    graphics::image(
      s_ci$muDelta,
      s_ci$logORdelta,
      s_ci$surface,
      col = rev(fields::tim.colors(256)),
      useRaster = FALSE,
      xlab = "Difference in means",
      ylab = "log OR",
      cex.axis = 1.2,
      cex.main = 1.5,
      cex.lab = 1.5
    )
    graphics::points(m_splrt$muDelta, log(m_splrt$alphaDelta), pch = 19, cex = 1.5)
    graphics::contour(
      s_ci$muDelta, s_ci$logORdelta, s_ci$surface,
      add = TRUE, levels = stats::qchisq(0.99, 2), lwd = 1, labels = 0.99, lty = 1
    )
    graphics::contour(
      s_ci$muDelta, s_ci$logORdelta, s_ci$surface,
      add = TRUE, levels = stats::qchisq(0.95, 2), lwd = 1, labels = 0.95, lty = 1
    )
    graphics::contour(
      s_ci$muDelta, s_ci$logORdelta, s_ci$surface,
      add = TRUE, levels = stats::qchisq(0.8, 2), lwd = 1, labels = 0.8, lty = 1
    )
    graphics::contour(
      s_ci$muDelta, s_ci$logORdelta, s_ci$surface,
      add = TRUE, levels = stats::qchisq(0.5, 2), lwd = 1, labels = 0.5, lty = 1
    )
    graphics::title("Confidence regions", cex.main = 1.5)
    graphics::abline(v = 0, lty = 2)
    graphics::abline(h = 0, lty = 2)
  })
}

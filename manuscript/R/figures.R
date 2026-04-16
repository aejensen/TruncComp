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

.simulation_method_palette <- function() {
  c(
    "Parametric LRT" = "#C0392B",
    "Semi-parametric LRT" = "#1F7A8C",
    "T-test" = "#3B6EA5",
    "Wilcoxon" = "#111111"
  )
}

.simulation_figure_theme <- function() {
  ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.title = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold"),
      plot.title.position = "plot"
    )
}

.simulation_main_power_data <- function(simulation_results) {
  strongest_h <- max(simulation_results$config$effect_levels)
  simulation_results$method_metrics[
    simulation_results$method_metrics$scenario_id %in% simulation_results$config$main_text_scenarios &
      simulation_results$method_metrics$h == strongest_h,
    ,
    drop = FALSE
  ]
}

build_power_curves_figure <- function(simulation_results, path) {
  plot_data <- .simulation_main_power_data(simulation_results)
  plot_data$short_label <- factor(
    plot_data$short_label,
    levels = unique(plot_data$short_label[order(plot_data$scenario_order)])
  )

  plot_object <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = n, y = reject_rate, color = method_label)
  ) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 1.2) +
    ggplot2::facet_wrap(~ short_label, ncol = 2) +
    ggplot2::scale_x_continuous(breaks = simulation_results$config$n_seq) +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::scale_color_manual(values = .simulation_method_palette()) +
    ggplot2::labs(
      x = "Observations per group",
      y = "Estimated power",
      title = "Power by sample size at the strongest non-null effect level"
    ) +
    .simulation_figure_theme()

  ggplot2::ggsave(path, plot = plot_object, width = 8.5, height = 5.3, device = grDevices::pdf)
}

build_power_effect_figure <- function(simulation_results, path) {
  plot_data <- simulation_results$method_metrics[
    simulation_results$method_metrics$scenario_id %in% simulation_results$config$main_text_scenarios &
      simulation_results$method_metrics$n == simulation_results$config$power_effect_n,
    ,
    drop = FALSE
  ]
  plot_data$short_label <- factor(
    plot_data$short_label,
    levels = unique(plot_data$short_label[order(plot_data$scenario_order)])
  )

  plot_object <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = h, y = reject_rate, color = method_label)
  ) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 1.2) +
    ggplot2::facet_wrap(~ short_label, ncol = 2) +
    ggplot2::scale_x_continuous(breaks = simulation_results$config$effect_levels) +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::scale_color_manual(values = .simulation_method_palette()) +
    ggplot2::labs(
      x = "Effect level h",
      y = "Estimated power",
      title = sprintf("Power by effect level at n = %d per group", simulation_results$config$power_effect_n)
    ) +
    .simulation_figure_theme()

  ggplot2::ggsave(path, plot = plot_object, width = 8.5, height = 5.3, device = grDevices::pdf)
}

build_supplementary_power_figure <- function(simulation_results, path) {
  strongest_h <- max(simulation_results$config$effect_levels)
  plot_data <- simulation_results$method_metrics[
    simulation_results$method_metrics$scenario_id %in% simulation_results$config$supplementary_scenarios &
      simulation_results$method_metrics$h == strongest_h,
    ,
    drop = FALSE
  ]
  plot_data$short_label <- factor(
    plot_data$short_label,
    levels = unique(plot_data$short_label[order(plot_data$scenario_order)])
  )

  plot_object <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = n, y = reject_rate, color = method_label)
  ) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 1.2) +
    ggplot2::facet_wrap(~ short_label, ncol = 2) +
    ggplot2::scale_x_continuous(breaks = simulation_results$config$n_seq) +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::scale_color_manual(values = .simulation_method_palette()) +
    ggplot2::labs(
      x = "Observations per group",
      y = "Estimated power",
      title = "Supplementary power curves for the reference and concordant-effect scenarios"
    ) +
    .simulation_figure_theme()

  ggplot2::ggsave(path, plot = plot_object, width = 8.5, height = 4.4, device = grDevices::pdf)
}

build_type1_curves_figure <- function(simulation_results, path) {
  plot_data <- simulation_results$method_metrics[simulation_results$method_metrics$is_null, , drop = FALSE]
  plot_data$short_label <- factor(
    plot_data$short_label,
    levels = unique(plot_data$short_label[order(plot_data$scenario_order)])
  )

  plot_object <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = n, y = reject_rate, color = method_label)
  ) +
    ggplot2::geom_hline(yintercept = 0.05, linetype = "dashed", color = "gray35", linewidth = 0.4) +
    ggplot2::geom_hline(yintercept = c(0.04, 0.06), linetype = "dotted", color = "gray55", linewidth = 0.35) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 1.1) +
    ggplot2::facet_wrap(~ short_label, ncol = 2) +
    ggplot2::scale_x_continuous(breaks = simulation_results$config$n_seq) +
    ggplot2::scale_y_continuous(limits = c(0, 0.12)) +
    ggplot2::scale_color_manual(values = .simulation_method_palette()) +
    ggplot2::labs(
      x = "Observations per group",
      y = "Empirical size at h = 0",
      title = "Scenario-matched null calibration across all six designs"
    ) +
    .simulation_figure_theme()

  ggplot2::ggsave(path, plot = plot_object, width = 8.5, height = 7.2, device = grDevices::pdf)
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

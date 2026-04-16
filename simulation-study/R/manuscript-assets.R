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

simulation_study_build_power_curves_figure <- function(simulation_results, path) {
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

simulation_study_build_power_effect_figure <- function(simulation_results, path) {
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

simulation_study_build_supplementary_power_figure <- function(simulation_results, path) {
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

simulation_study_build_type1_curves_figure <- function(simulation_results, path) {
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

simulation_study_write_kable_table <- function(table_object, path) {
  writeLines(as.character(table_object), con = path, useBytes = TRUE)
  invisible(path)
}

.simulation_power_table_data <- function(simulation_results, scenario_id) {
  power_wide <- simulation_study_metric_wide(simulation_results, "reject_rate")
  fail_wide <- simulation_study_metric_wide(simulation_results, "fail_rate")

  power_wide <- power_wide[power_wide$scenario_id == scenario_id, , drop = FALSE]
  fail_wide <- fail_wide[fail_wide$scenario_id == scenario_id, c("cell_id", simulation_study_methods()), drop = FALSE]
  merged <- merge(power_wide, fail_wide, by = "cell_id", suffixes = c("_power", "_fail"), sort = FALSE)
  merged$max_fail <- apply(
    merged[, paste0(simulation_study_methods(), "_fail"), drop = FALSE],
    1,
    function(x) max(x, na.rm = TRUE)
  )

  out <- merged[, c(
    "n", "h", "effect_text", "expected_combined_delta",
    paste0(simulation_study_methods(), "_power"), "max_fail"
  )]
  names(out) <- c(
    "n", "h", "Effect", "Expected $\\Delta$",
    unname(simulation_study_method_labels()),
    "Max fail"
  )
  out <- out[order(out$h, out$n), , drop = FALSE]
  out[, unname(simulation_study_method_labels())] <- out[, unname(simulation_study_method_labels())] * 100
  out$`Max fail` <- out$`Max fail` * 100
  out
}

.simulation_type1_table_data <- function(simulation_results) {
  power_wide <- simulation_study_metric_wide(simulation_results, "reject_rate")
  fail_wide <- simulation_study_metric_wide(simulation_results, "fail_rate")

  power_wide <- power_wide[power_wide$is_null, , drop = FALSE]
  fail_wide <- fail_wide[fail_wide$is_null, c("cell_id", simulation_study_methods()), drop = FALSE]
  merged <- merge(power_wide, fail_wide, by = "cell_id", suffixes = c("_power", "_fail"), sort = FALSE)
  merged$max_fail <- apply(
    merged[, paste0(simulation_study_methods(), "_fail"), drop = FALSE],
    1,
    function(x) max(x, na.rm = TRUE)
  )

  out <- merged[, c(
    "scenario_id", "short_label", "n",
    paste0(simulation_study_methods(), "_power"), "max_fail"
  )]
  names(out) <- c(
    "Scenario", "Label", "n",
    unname(simulation_study_method_labels()),
    "Max fail"
  )
  out <- out[order(out$Scenario, out$n), , drop = FALSE]
  out[, unname(simulation_study_method_labels())] <- out[, unname(simulation_study_method_labels())] * 100
  out$`Max fail` <- out$`Max fail` * 100
  out
}

simulation_study_build_scenario_table <- function(simulation_results, output_dir) {
  if (!requireNamespace("kableExtra", quietly = TRUE)) {
    stop("The kableExtra package is required to build manuscript tables.", call. = FALSE)
  }

  scenario_table <- simulation_results$scenarios[, c("scenario_id", "title", "atom_latex", "survivor_latex")]
  names(scenario_table) <- c("Scenario", "Configuration", "$P(A = 1 \\mid R)$", "$Y \\mid A = 1, R$")

  table_object <- kableExtra::kbl(
    scenario_table,
    format = "latex",
    booktabs = TRUE,
    escape = FALSE,
    linesep = rep("", nrow(scenario_table)),
    caption = "Simulation scenarios for the TruncComp2 power study"
  )
  table_object <- kableExtra::kable_styling(
    table_object,
    latex_options = c("hold_position", "scale_down")
  )
  simulation_study_write_kable_table(table_object, file.path(output_dir, "simulation-scenarios.tex"))
}

simulation_study_build_supplementary_tables <- function(simulation_results, output_dir) {
  if (!requireNamespace("kableExtra", quietly = TRUE)) {
    stop("The kableExtra package is required to build manuscript tables.", call. = FALSE)
  }

  scenario_rows <- split(simulation_results$scenarios, simulation_results$scenarios$scenario_id)
  for (scenario_id in names(scenario_rows)) {
    scenario_info <- scenario_rows[[scenario_id]][1, , drop = FALSE]
    power_table <- .simulation_power_table_data(simulation_results, scenario_id = scenario_id)

    table_object <- kableExtra::kbl(
      power_table,
      format = "latex",
      digits = c(0, 0, 0, 3, 2, 2, 2, 2, 2),
      booktabs = TRUE,
      escape = FALSE,
      longtable = TRUE,
      linesep = rep("", nrow(power_table)),
      caption = sprintf("%s: estimated power by sample size and effect level", scenario_info$title)
    )
    table_object <- kableExtra::kable_styling(
      table_object,
      latex_options = "repeat_header"
    )

    simulation_study_write_kable_table(
      table_object,
      file.path(output_dir, sprintf("supplementary-power-%s.tex", tolower(scenario_id)))
    )
  }

  type1_table <- .simulation_type1_table_data(simulation_results)
  type1_table$Label <- NULL

  type1_object <- kableExtra::kbl(
    type1_table,
    format = "latex",
    digits = c(0, 0, 2, 2, 2, 2, 2),
    booktabs = TRUE,
    longtable = TRUE,
    linesep = rep("", nrow(type1_table)),
    caption = "Scenario-matched null calibration at effect level $h = 0$"
  )
  type1_object <- kableExtra::kable_styling(
    type1_object,
    latex_options = "repeat_header"
  )
  simulation_study_write_kable_table(type1_object, file.path(output_dir, "supplementary-type1.tex"))

  invisible(output_dir)
}

simulation_study_build_manuscript_assets <- function(repo_root, simulation_results) {
  manuscript_build_dir <- ensure_dir(file.path(repo_root, "manuscript", "build"))
  figures_dir <- ensure_dir(file.path(manuscript_build_dir, "figures"))
  tables_dir <- ensure_dir(file.path(manuscript_build_dir, "tables"))

  simulation_study_build_power_curves_figure(
    simulation_results,
    file.path(figures_dir, "power-curves.pdf")
  )
  simulation_study_build_power_effect_figure(
    simulation_results,
    file.path(figures_dir, "power-effects.pdf")
  )
  simulation_study_build_supplementary_power_figure(
    simulation_results,
    file.path(figures_dir, "supplementary-power-curves.pdf")
  )
  simulation_study_build_type1_curves_figure(
    simulation_results,
    file.path(figures_dir, "type1-curves.pdf")
  )

  simulation_study_build_scenario_table(simulation_results, tables_dir)
  simulation_study_build_supplementary_tables(simulation_results, tables_dir)

  invisible(list(figures_dir = figures_dir, tables_dir = tables_dir))
}

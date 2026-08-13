.simulation_method_palette <- function() {
  c(
    "Parametric LRT" = "#D55E00",
    "Semiparametric LRT" = "#0072B2",
    "Welch t-test" = "#009E73",
    "Wilcoxon" = "#CC79A7"
  )
}

.simulation_figure_method_label <- function(x) {
  x <- as.character(x)
  x[x == paste0("T", "-", "test")] <- "Welch t-test"
  x[x == paste0("Semi", "-", "parametric LRT")] <- "Semiparametric LRT"
  factor(x, levels = names(.simulation_method_palette()))
}

.simulation_figure_theme <- function() {
  ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.title = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 12),
      plot.title.position = "plot"
    )
}

.simulation_table_method_labels <- function() {
  c("Wilcoxon (\\%)", "Welch $t$-test (\\%)", "LRT (\\%)", "SPLRT (\\%)")
}

.simulation_current_short_labels <- function(scenario_ids) {
  scenarios <- simulation_study_scenarios()
  vapply(
    as.character(scenario_ids),
    function(scenario_id) scenarios[[scenario_id]]$short_label,
    character(1)
  )
}

.simulation_power_data_all_scenarios <- function(simulation_results) {
  strongest_h <- max(simulation_results$config$effect_levels)
  simulation_results$method_metrics[
    simulation_results$method_metrics$h == strongest_h,
    ,
    drop = FALSE
  ]
}

simulation_study_build_power_curves_figure <- function(simulation_results, path) {
  plot_data <- .simulation_power_data_all_scenarios(simulation_results)
  plot_data <- plot_data[is.finite(plot_data$reject_rate), , drop = FALSE]
  plot_data$method_label <- .simulation_figure_method_label(plot_data$method_label)
  plot_data$short_label <- .simulation_current_short_labels(plot_data$scenario_id)
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
    ggplot2::facet_wrap(~ short_label, ncol = 3) +
    ggplot2::scale_x_continuous(breaks = simulation_results$config$n_seq) +
    ggplot2::scale_y_continuous() +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::scale_color_manual(values = .simulation_method_palette()) +
    ggplot2::labs(
      x = "Observations per group",
      y = "Estimated rejection rate",
      title = "Rejection rates by sample size at effect level h = 3 across all six scenarios"
    ) +
    .simulation_figure_theme()

  ggplot2::ggsave(path, plot = plot_object, width = 10, height = 6.8, device = grDevices::pdf)
}

simulation_study_build_power_effect_figure <- function(simulation_results, path) {
  plot_data <- simulation_results$method_metrics[
    simulation_results$method_metrics$n == simulation_results$config$power_effect_n,
    ,
    drop = FALSE
  ]
  plot_data <- plot_data[is.finite(plot_data$reject_rate), , drop = FALSE]
  plot_data$method_label <- .simulation_figure_method_label(plot_data$method_label)
  plot_data$short_label <- .simulation_current_short_labels(plot_data$scenario_id)
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
    ggplot2::facet_wrap(~ short_label, ncol = 3) +
    ggplot2::scale_x_continuous(breaks = simulation_results$config$effect_levels) +
    ggplot2::scale_y_continuous() +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::scale_color_manual(values = .simulation_method_palette()) +
    ggplot2::labs(
      x = "Effect level h",
      y = "Estimated rejection rate",
      title = sprintf(
        "Rejection rates by effect level at n = %d per group across all six scenarios",
        simulation_results$config$power_effect_n
      )
    ) +
    .simulation_figure_theme()

  ggplot2::ggsave(path, plot = plot_object, width = 10, height = 6.8, device = grDevices::pdf)
}

simulation_study_build_type1_curves_figure <- function(simulation_results, path) {
  plot_data <- simulation_results$method_metrics[simulation_results$method_metrics$is_null, , drop = FALSE]
  plot_data <- plot_data[is.finite(plot_data$reject_rate), , drop = FALSE]
  plot_data$method_label <- .simulation_figure_method_label(plot_data$method_label)
  plot_data$short_label <- .simulation_current_short_labels(plot_data$scenario_id)
  plot_data$short_label <- factor(
    plot_data$short_label,
    levels = unique(plot_data$short_label[order(plot_data$scenario_order)])
  )

  plot_object <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = n, y = reject_rate, color = method_label)
  ) +
    ggplot2::geom_hline(yintercept = 0.05, linetype = "dashed", color = "gray35", linewidth = 0.4) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 1.1) +
    ggplot2::facet_wrap(~ short_label, ncol = 2, scales = "free_y") +
    ggplot2::scale_x_continuous(breaks = simulation_results$config$n_seq) +
    ggplot2::scale_y_continuous() +
    ggplot2::scale_color_manual(values = .simulation_method_palette()) +
    ggplot2::labs(
      x = "Observations per group",
      y = "Null rejection rate at h = 0",
      title = "Null rejection rates at h = 0 across all six scenarios"
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

  power_wide <- power_wide[power_wide$scenario_id == scenario_id, , drop = FALSE]
  scenario <- simulation_study_scenarios()[[scenario_id]]
  power_wide$effect_text <- vapply(
    power_wide$h,
    function(h) scenario$effect_parameters(as.integer(h))$effect_text,
    character(1)
  )

  out <- power_wide[, c(
    "n", "h", "effect_text", "expected_combined_delta",
    simulation_study_methods()
  )]
  out <- out[order(out$h, out$n), , drop = FALSE]
  names(out) <- c(
    "$n_{\\mathrm{arm}}$", "$h$", "Effect", "Expected $\\Delta^{(0)}$",
    .simulation_table_method_labels()
  )
  out$Effect <- gsub("\\\\", intToUtf8(92), out$Effect, fixed = TRUE)
  out$Effect <- sprintf("\\makecell[l]{%s}", out$Effect)
  out[, .simulation_table_method_labels()] <- out[, .simulation_table_method_labels()] * 100
  row.names(out) <- NULL
  out
}

.simulation_type1_table_data <- function(simulation_results) {
  power_wide <- simulation_study_metric_wide(simulation_results, "reject_rate")

  power_wide <- power_wide[power_wide$is_null, , drop = FALSE]
  power_wide$short_label <- .simulation_current_short_labels(power_wide$scenario_id)

  out <- power_wide[, c(
    "scenario_id", "short_label", "n",
    simulation_study_methods()
  )]
  out <- out[order(out$scenario_id, out$n), , drop = FALSE]
  names(out) <- c(
    "Scenario", "Label", "$n_{\\mathrm{arm}}$",
    .simulation_table_method_labels()
  )
  out[, .simulation_table_method_labels()] <- out[, .simulation_table_method_labels()] * 100
  row.names(out) <- NULL
  out
}

simulation_study_build_scenario_table <- function(simulation_results, output_dir) {
  if (!requireNamespace("kableExtra", quietly = TRUE)) {
    stop("The kableExtra package is required to build manuscript tables.", call. = FALSE)
  }

  scenario_table <- simulation_results$scenarios[, c("scenario_id", "title")]
  names(scenario_table) <- c("Scenario", "Configuration")
  configuration_labels <- c(
    S1 = "Normal non-atom shift,\\\\ low atom probability",
    S2 = "Effect only on probability\\\\ outside atom",
    S3 = "Antagonistic\\\\ cancellation",
    S4 = "Concordant joint\\\\ improvement",
    S5 = "Gamma shape change\\\\ with mean shift",
    S6 = "Rare outlier\\\\ contamination"
  )
  scenario_table$Configuration <-
    sprintf("\\makecell[l]{%s}", configuration_labels[scenario_table$Scenario])

  nonatom_probability <- c(
    S1 = "$\\pi_{0,h} = \\pi_{1,h} = 0.85$",
    S2 = paste(
      "$\\pi_{0,h} = 0.45$",
      "$(\\pi_{1,0},\\ldots,\\pi_{1,3}) = (0.45, 0.50, 0.55, 0.60)$",
      sep = "\\\\ "
    ),
    S3 = paste(
      "$\\pi_{0,h} = 0.50$",
      "$(\\pi_{1,0},\\ldots,\\pi_{1,3}) = (0.50, 0.525, 0.55, 0.60)$",
      sep = "\\\\ "
    ),
    S4 = paste(
      "$\\pi_{0,h} = 0.50$",
      "$(\\pi_{1,0},\\ldots,\\pi_{1,3}) = (0.50, 0.55, 0.60, 0.65)$",
      sep = "\\\\ "
    ),
    S5 = "$\\pi_{0,h} = \\pi_{1,h} = 0.65$",
    S6 = "$\\pi_{0,h} = \\pi_{1,h} = 0.80$"
  )
  scenario_table[["\\makecell[l]{Probability outside atom\\\\ $\\pi_{r,h}=\\Prob(A=1\\mid R=r)$}"]] <-
    sprintf("\\makecell[l]{%s}", nonatom_probability[scenario_table$Scenario])

  non_atom_law <- c(
    S1 = paste(
      "$F_{0,h}^c = N(3, 1)$",
      "$F_{1,h}^c = N(3 + \\delta_h, 1)$",
      "$(\\delta_0,\\ldots,\\delta_3) = (0, 0.20, 0.35, 0.50)$",
      sep = "\\\\ "
    ),
    S2 = "$F_{0,h}^c = F_{1,h}^c = N(3, 1)$",
    S3 = paste(
      "$F_{0,h}^c = N(3, 1)$",
      "$F_{1,h}^c = N(1.5 / \\pi_{1,h}, 1)$",
      sep = "\\\\ "
    ),
    S4 = paste(
      "$F_{0,h}^c = N(3, 1)$",
      "$F_{1,h}^c = N(3 + 0.15h, 1)$",
      sep = "\\\\ "
    ),
    S5 = paste(
      "$F_{0,h}^c = \\Gamma(9, 3/9)$",
      "$F_{1,h}^c = \\Gamma(k_h, m_h / k_h)$",
      "$(k_0,\\ldots,k_3) = (9, 6, 4, 2.5)$",
      "$(m_0,\\ldots,m_3) = (3.0, 3.2, 3.4, 3.6)$",
      sep = "\\\\ "
    ),
    S6 = paste(
      "$F_{r,h}^c = 0.95 N(3 + r\\delta_h, 0.35^2)$",
      "$\\hspace{2.6em}{}+ 0.05 N(10 + r\\delta_h, 0.5^2)$",
      "$(\\delta_0,\\ldots,\\delta_3) = (0, 0.15, 0.25, 0.35)$",
      sep = "\\\\ "
    )
  )
  scenario_table[["\\makecell[l]{Non-atom law\\\\ $F_{r,h}^c$}"]] <-
    sprintf("\\makecell[l]{%s}", non_atom_law[scenario_table$Scenario])

  table_object <- kableExtra::kbl(
    scenario_table,
    format = "latex",
    booktabs = TRUE,
    escape = FALSE,
    linesep = rep("", nrow(scenario_table)),
    caption = "Simulation scenarios S1--S6 for effect levels $h=0,1,2,3$",
    label = "simulation"
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

  scenario_captions <- c(
    S1 = "Scenario S1 (Normal non-atom shift, low atom probability)",
    S2 = "Scenario S2 (effect only on probability outside atom)",
    S3 = "Scenario S3 (antagonistic cancellation)",
    S4 = "Scenario S4 (concordant joint improvement)",
    S5 = "Scenario S5 (Gamma shape change with non-atom mean shift)",
    S6 = "Scenario S6 (rare outlier contamination)"
  )

  scenario_rows <- split(simulation_results$scenarios, simulation_results$scenarios$scenario_id)
  for (scenario_id in names(scenario_rows)) {
    power_table <- .simulation_power_table_data(simulation_results, scenario_id = scenario_id)
    caption <- scenario_captions[[scenario_id]]
    if (is.null(caption)) {
      caption <- sprintf("Scenario %s", scenario_id)
    }

    table_object <- kableExtra::kbl(
      power_table,
      format = "latex",
      digits = c(0, 0, 0, 3, 2, 2, 2, 2),
      booktabs = TRUE,
      escape = FALSE,
      longtable = TRUE,
      linesep = rep("", nrow(power_table)),
      caption = sprintf("%s: estimated rejection percentages by sample size and effect level", caption)
    )
    table_object <- kableExtra::kable_styling(
      table_object,
      latex_options = "repeat_header",
      font_size = 7
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
    digits = c(0, 0, 2, 2, 2, 2),
    booktabs = TRUE,
    escape = FALSE,
    longtable = TRUE,
    linesep = rep("", nrow(type1_table)),
    caption = "Scenarios S1--S6: scenario-matched null calibration at effect level $h = 0$"
  )
  type1_object <- kableExtra::kable_styling(
    type1_object,
    latex_options = "repeat_header",
    font_size = 7
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
  simulation_study_build_type1_curves_figure(
    simulation_results,
    file.path(figures_dir, "type1-curves.pdf")
  )

  simulation_study_build_scenario_table(simulation_results, tables_dir)
  simulation_study_build_supplementary_tables(simulation_results, tables_dir)

  invisible(list(figures_dir = figures_dir, tables_dir = tables_dir))
}

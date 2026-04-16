write_kable_table <- function(table_object, path) {
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

build_simulation_scenario_table <- function(simulation_results, output_dir) {
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
  write_kable_table(table_object, file.path(output_dir, "simulation-scenarios.tex"))
}

build_supplementary_tables <- function(simulation_results, output_dir) {
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

    write_kable_table(
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
  write_kable_table(type1_object, file.path(output_dir, "supplementary-type1.tex"))

  invisible(output_dir)
}

write_kable_table <- function(table_object, path) {
  writeLines(as.character(table_object), con = path, useBytes = TRUE)
  invisible(path)
}

build_supplementary_tables <- function(simulation_results, output_dir) {
  if (!requireNamespace("kableExtra", quietly = TRUE)) {
    stop("The kableExtra package is required to build manuscript tables.", call. = FALSE)
  }

  power1 <- set_method_colnames(simulation_results$power1)
  power2 <- set_method_colnames(simulation_results$power2)
  power3 <- set_method_colnames(simulation_results$power3)
  power4 <- set_method_colnames(simulation_results$power4)
  null1 <- lapply(simulation_results$null1, set_method_colnames)
  null2 <- lapply(simulation_results$null2, set_method_colnames)

  tbl_power_1_2 <- kableExtra::kbl(
    cbind(power1, power2) * 100,
    format = "latex",
    digits = 2,
    booktabs = TRUE,
    linesep = rep("", nrow(power1)),
    caption = "Power simulation - scenarios 1 and 2"
  )
  tbl_power_1_2 <- kableExtra::kable_styling(tbl_power_1_2, latex_options = "hold_position")
  tbl_power_1_2 <- kableExtra::add_header_above(
    tbl_power_1_2,
    c(" " = 1, "Scenario 1" = 4, "Scenario 2" = 4),
    bold = TRUE
  )
  write_kable_table(tbl_power_1_2, file.path(output_dir, "supplementary-power-1-2.tex"))

  tbl_power_3_4 <- kableExtra::kbl(
    cbind(power3, power4) * 100,
    format = "latex",
    digits = 2,
    booktabs = TRUE,
    linesep = rep("", nrow(power3)),
    caption = "Power simulation - scenarios 3 and 4"
  )
  tbl_power_3_4 <- kableExtra::kable_styling(tbl_power_3_4, latex_options = "hold_position")
  tbl_power_3_4 <- kableExtra::add_header_above(
    tbl_power_3_4,
    c(" " = 1, "Scenario 3" = 4, "Scenario 4" = 4),
    bold = TRUE
  )
  write_kable_table(tbl_power_3_4, file.path(output_dir, "supplementary-power-3-4.tex"))

  tbl_type1_1 <- kableExtra::kbl(
    do.call("cbind", null1) * 100,
    format = "latex",
    digits = 2,
    booktabs = TRUE,
    linesep = rep("", nrow(null1[[1]])),
    caption = "Type I error simulation - scenario 1"
  )
  tbl_type1_1 <- kableExtra::kable_styling(
    tbl_type1_1,
    latex_options = c("hold_position", "scale_down")
  )
  tbl_type1_1 <- kableExtra::add_header_above(
    tbl_type1_1,
    c(" " = 1, "$\\pi = 0.2$" = 4, "$\\pi = 0.4$" = 4, "$\\pi = 0.6$" = 4, "$\\pi = 0.8$" = 4),
    bold = TRUE,
    escape = FALSE
  )
  write_kable_table(tbl_type1_1, file.path(output_dir, "supplementary-type1-1.tex"))

  tbl_type1_2 <- kableExtra::kbl(
    do.call("cbind", null2) * 100,
    format = "latex",
    digits = 2,
    booktabs = TRUE,
    linesep = rep("", nrow(null2[[1]])),
    caption = "Type I error simulation - scenario 2"
  )
  tbl_type1_2 <- kableExtra::kable_styling(
    tbl_type1_2,
    latex_options = c("hold_position", "scale_down")
  )
  tbl_type1_2 <- kableExtra::add_header_above(
    tbl_type1_2,
    c(" " = 1, "$\\pi = 0.2$" = 4, "$\\pi = 0.4$" = 4, "$\\pi = 0.6$" = 4, "$\\pi = 0.8$" = 4),
    bold = TRUE,
    escape = FALSE
  )
  write_kable_table(tbl_type1_2, file.path(output_dir, "supplementary-type1-2.tex"))

  invisible(output_dir)
}

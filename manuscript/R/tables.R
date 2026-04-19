write_kable_table <- function(table_object, path) {
  writeLines(as.character(table_object), con = path, useBytes = TRUE)
  invisible(path)
}

write_latex_table <- function(path, caption, label, header, rows, align = NULL, size = "\\small") {
  if (is.null(align)) {
    align <- paste0("l", paste(rep("c", length(header) - 1L), collapse = ""))
  }
  body <- apply(rows, 1, function(x) paste(x, collapse = " & "))
  lines <- c(
    "\\begin{table}[htbp]",
    "\\centering",
    sprintf("\\caption{%s}", caption),
    sprintf("\\label{%s}", label),
    size,
    sprintf("\\begin{tabular}{%s}", align),
    "\\toprule",
    paste(header, collapse = " & "),
    "\\\\",
    "\\midrule",
    paste0(body, "\\\\"),
    "\\bottomrule",
    "\\end{tabular}",
    "\\end{table}"
  )

  writeLines(lines, con = path, useBytes = TRUE)
  invisible(path)
}

write_application_descriptive_table <- function(application_results, path) {
  s <- application_results$summary_stats
  m <- application_results$metadata
  control <- s[s$R == 0L, , drop = FALSE]
  treatment <- s[s$R == 1L, , drop = FALSE]

  by_group <- function(var) {
    c(control[[var]], treatment[[var]])
  }
  n_pct <- function(n, p) sprintf("%s (%s)", format_count(n), format_percent(p))

  rows <- rbind(
    c("Randomized participants", format_count(m$randomized_group_n[[1]]), format_count(m$randomized_group_n[[2]])),
    c("Included in analysis", format_count(by_group("n")[[1]]), format_count(by_group("n")[[2]])),
    c(
      "Alive but missing 6-month EQ-VAS",
      format_count(m$excluded_alive_missing_euroqol6_by_group[[1]]),
      format_count(m$excluded_alive_missing_euroqol6_by_group[[2]])
    ),
    c(
      "Death atom, n (\\%)",
      n_pct(control$atom_n, control$atom_prop),
      n_pct(treatment$atom_n, treatment$atom_prop)
    ),
    c(
      "Alive with EQ-VAS, n (\\%)",
      n_pct(control$non_atom_n, control$non_atom_prop),
      n_pct(treatment$non_atom_n, treatment$non_atom_prop)
    ),
    c(
      "Survivor EQ-VAS, mean (SD)",
      sprintf("%s (%s)", format_number(control$survivor_mean, 1), format_number(control$survivor_sd, 1)),
      sprintf("%s (%s)", format_number(treatment$survivor_mean, 1), format_number(treatment$survivor_sd, 1))
    ),
    c(
      "Survivor EQ-VAS, median (IQR)",
      sprintf(
        "%s (%s, %s)",
        format_number(control$survivor_median, 0),
        format_number(control$survivor_q1, 0),
        format_number(control$survivor_q3, 0)
      ),
      sprintf(
        "%s (%s, %s)",
        format_number(treatment$survivor_median, 0),
        format_number(treatment$survivor_q1, 0),
        format_number(treatment$survivor_q3, 0)
      )
    ),
    c(
      "Survivor EQ-VAS, range",
      sprintf("%s to %s", format_number(control$survivor_min, 0), format_number(control$survivor_max, 0)),
      sprintf("%s to %s", format_number(treatment$survivor_min, 0), format_number(treatment$survivor_max, 0))
    ),
    c(
      "Combined endpoint mean",
      format_number(control$combined_mean, 1),
      format_number(treatment$combined_mean, 1)
    )
  )

  write_latex_table(
    path = path,
    caption = "IST-3 analysis sample and endpoint summaries.",
    label = "tab:ist3-descriptives",
    header = c("Summary", m$group_labels[[1]], m$group_labels[[2]]),
    rows = rows,
    align = "lcc"
  )
}

write_application_frequentist_table <- function(application_results, path) {
  tests <- application_results$standard_tests
  lrt <- application_results$model_lrt
  splrt <- application_results$model_splrt

  std_row <- function(analysis, method, estimate = "") {
    row <- tests[tests$analysis == analysis & tests$method == method, , drop = FALSE]
    c(
      analysis,
      method,
      estimate,
      ifelse(is.na(row$statistic), "", format_number(row$statistic, 2)),
      format_p_value(row$p_value)
    )
  }
  model_row <- function(label, fit) {
    c(
      label,
      "Joint likelihood-ratio test",
      sprintf(
        "$\\mu_\\delta = %s$, $\\alpha_\\delta = %s$, $\\Delta = %s$",
        format_number(fit$mu_delta, 2),
        format_number(fit$alpha_delta, 3),
        format_number(fit$delta, 2)
      ),
      format_number(fit$statistic, 2),
      format_p_value(fit$p.value)
    )
  }

  delta_hat <- application_results$contrasts$estimate[application_results$contrasts$contrast == "Combined mean contrast Delta"]
  mu_hat <- application_results$contrasts$estimate[application_results$contrasts$contrast == "Survivor EQ-VAS mean difference"]

  rows <- rbind(
    std_row("Combined endpoint", "Welch t-test", sprintf("$\\Delta = %s$", format_number(delta_hat, 2))),
    std_row("Combined endpoint", "Wilcoxon rank-sum", "Death tied as worst value"),
    std_row("Survivors only", "Welch t-test", sprintf("$\\mu_\\delta = %s$", format_number(mu_hat, 2))),
    std_row("Survivors only", "Wilcoxon rank-sum", "Conditional on survival"),
    std_row("Death atom", "Fisher exact test", "Death proportions"),
    model_row("Parametric model", lrt),
    model_row("Semiparametric model", splrt)
  )

  write_latex_table(
    path = path,
    caption = "Standard and likelihood-ratio analyses of the IST-3 endpoint.",
    label = "tab:ist3-frequentist",
    header = c("Analysis", "Method", "Estimate/target", "Statistic", "$p$-value"),
    rows = rows,
    align = "lllcc",
    size = "\\scriptsize"
  )
}

write_application_delta_ci_table <- function(application_results, path) {
  intervals <- application_results$delta_intervals
  method_label <- c(welch = "Welch", profile = "Profile likelihood", projected = "Projected joint region")
  rows <- t(apply(intervals, 1, function(row) {
    estimated <- identical(row[["status"]], "estimated")
    c(
      method_label[[row[["method"]]]],
      if (estimated) format_number(as.numeric(row[["lower"]]), 2) else "--",
      if (estimated) format_number(as.numeric(row[["upper"]]), 2) else "--",
      if (estimated) "Estimated" else row[["status"]]
    )
  }))

  write_latex_table(
    path = path,
    caption = "Confidence intervals for the combined mean contrast $\\Delta$ from the semiparametric IST-3 fit.",
    label = "tab:ist3-delta-ci",
    header = c("Interval construction", "Lower", "Upper", "Status"),
    rows = rows,
    align = "lccc"
  )
}

write_application_bayes_summary_table <- function(application_results, path) {
  bayes <- application_results$bayes
  if (!isTRUE(bayes$success)) {
    rows <- rbind(c("Bayesian model", "--", "--", latex_escape(paste("Unavailable:", bayes$error))))
    return(write_latex_table(
      path = path,
      caption = "Bayesian IST-3 posterior summaries.",
      label = "tab:ist3-bayes",
      header = c("Quantity", "Posterior median", "95\\% CrI", "Probability"),
      rows = rows,
      align = "llll",
      size = "\\small"
    ))
  }

  arm <- bayes$arm_table
  contrast <- bayes$summary_table
  probs <- bayes$probabilities
  diagnostics <- bayes$diagnostics
  interval <- function(tab, row) {
    sprintf("%s to %s", format_number(tab[row, "conf.low"], 2), format_number(tab[row, "conf.high"], 2))
  }
  diagnostic_text <- sprintf(
    "Divergences %s; max $\\hat R$ %s; min bulk ESS %s",
    format_count(diagnostics$divergences),
    format_number(diagnostics$max_rhat, 3),
    format_number(diagnostics$min_bulk_ess, 0)
  )
  truncation_text <- if (!is.null(diagnostics$truncation)) {
    sprintf(
      "Tail q95 %s/%s; passed = %s",
      format_number(diagnostics$truncation$q95[["q95_tail_0"]], 3),
      format_number(diagnostics$truncation$q95[["q95_tail_1"]], 3),
      ifelse(isTRUE(diagnostics$truncation_ok), "yes", "no")
    )
  } else {
    "Not assessed"
  }

  rows <- rbind(
    c("Death probability, control", format_number(arm["rho_0", "estimate"], 3), interval(arm, "rho_0"), ""),
    c("Death probability, rt-PA", format_number(arm["rho_1", "estimate"], 3), interval(arm, "rho_1"), ""),
    c("Alive/EQ-VAS probability, control", format_number(arm["pi_0", "estimate"], 3), interval(arm, "pi_0"), ""),
    c("Alive/EQ-VAS probability, rt-PA", format_number(arm["pi_1", "estimate"], 3), interval(arm, "pi_1"), ""),
    c("Survivor mean EQ-VAS, control", format_number(arm["mu_0_c", "estimate"], 2), interval(arm, "mu_0_c"), ""),
    c("Survivor mean EQ-VAS, rt-PA", format_number(arm["mu_1_c", "estimate"], 2), interval(arm, "mu_1_c"), ""),
    c("Death risk difference", format_number(contrast["delta_atom", "estimate"], 3), interval(contrast, "delta_atom"), sprintf("$\\Pr(<0) = %s$", format_probability(probs[["death_risk_reduction"]]))),
    c("$\\mu_\\delta$", format_number(contrast["mu_delta", "estimate"], 2), interval(contrast, "mu_delta"), sprintf("$\\Pr(>0) = %s$", format_probability(probs[["mu_delta_gt_0"]]))),
    c("$\\alpha_\\delta$", format_number(contrast["alpha_delta", "estimate"], 3), interval(contrast, "alpha_delta"), sprintf("$\\Pr(>1) = %s$", format_probability(probs[["alpha_delta_gt_1"]]))),
    c("$\\Delta$", format_number(contrast["delta", "estimate"], 2), interval(contrast, "delta"), sprintf("$\\Pr(>0) = %s$", format_probability(probs[["delta_gt_0"]]))),
    c("Sampler diagnostics", "--", "--", diagnostic_text),
    c("Mixture truncation check", "--", "--", truncation_text)
  )

  write_latex_table(
    path = path,
    caption = "Bayesian IST-3 posterior summaries from the two-part mixture model.",
    label = "tab:ist3-bayes",
    header = c("Quantity", "Posterior median", "95\\% CrI", "Probability"),
    rows = rows,
    align = "llll",
    size = "\\scriptsize"
  )
}

write_application_ppc_table <- function(application_results, path) {
  bayes <- application_results$bayes
  if (!isTRUE(bayes$success) || is.null(bayes$ppc_table)) {
    rows <- rbind(c(
      "Posterior predictive check",
      "--",
      "--",
      latex_escape(ifelse(is.null(bayes$ppc_error), "Unavailable", bayes$ppc_error))
    ))
    return(write_latex_table(
      path = path,
      caption = "Posterior predictive checks for the IST-3 Bayesian model.",
      label = "tab:ist3-ppc",
      header = c("Component", "Statistic", "Posterior predictive $p$", "Notes"),
      rows = rows,
      align = "llll"
    ))
  }

  ppc <- bayes$ppc_table
  component_label <- c(atom = "Death atom", continuous = "Survivor EQ-VAS")
  rows <- t(vapply(rownames(ppc), function(name) {
    c(
      component_label[[name]],
      latex_escape(ppc[name, "statistic"]),
      format_p_value(ppc[name, "p_value"]),
      sprintf("%s draws; %s scale", format_count(ppc[name, "ndraws"]), latex_escape(ppc[name, "scale"]))
    )
  }, character(4)))

  write_latex_table(
    path = path,
    caption = "Posterior predictive checks for the IST-3 Bayesian model.",
    label = "tab:ist3-ppc",
    header = c("Component", "Statistic", "Posterior predictive $p$", "Notes"),
    rows = rows,
    align = "p{0.16\\textwidth}p{0.31\\textwidth}p{0.16\\textwidth}p{0.18\\textwidth}",
    size = "\\scriptsize"
  )
}

write_application_tables <- function(application_results, tables_dir) {
  write_application_descriptive_table(application_results, file.path(tables_dir, "application-descriptives.tex"))
  write_application_frequentist_table(application_results, file.path(tables_dir, "application-frequentist.tex"))
  write_application_delta_ci_table(application_results, file.path(tables_dir, "application-delta-ci.tex"))
  write_application_bayes_summary_table(application_results, file.path(tables_dir, "application-bayes-summary.tex"))
  write_application_ppc_table(application_results, file.path(tables_dir, "application-ppc.tex"))
  invisible(TRUE)
}

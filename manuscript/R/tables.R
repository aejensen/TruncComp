write_kable_table <- function(table_object, path) {
  writeLines(as.character(table_object), con = path, useBytes = TRUE)
  invisible(path)
}

write_latex_table <- function(path, caption, label, header, rows, align = NULL, size = "\\small", placement = "htbp") {
  if (is.null(align)) {
    align <- paste0("l", paste(rep("c", length(header) - 1L), collapse = ""))
  }
  body <- apply(rows, 1, function(x) paste(x, collapse = " & "))
  lines <- c(
    sprintf("\\begin{table}[%s]", placement),
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
    c(
      "Alive but missing 6-month EQ-VAS",
      format_count(m$excluded_alive_missing_euroqol6_by_group[[1]]),
      format_count(m$excluded_alive_missing_euroqol6_by_group[[2]])
    ),
    c("Included in analysis", format_count(by_group("n")[[1]]), format_count(by_group("n")[[2]])),
    c(
      "Death atom, n (\\%)",
      n_pct(control$atom_n, control$atom_prop),
      n_pct(treatment$atom_n, treatment$atom_prop)
    ),
    c(
      "Alive with observed EQ-VAS, n (\\%)",
      n_pct(control$non_atom_n, control$non_atom_prop),
      n_pct(treatment$non_atom_n, treatment$non_atom_prop)
    ),
    c(
      "EQ-VAS, alive with observed score, mean (SD)",
      sprintf("%s (%s)", format_number(control$survivor_mean, 1), format_number(control$survivor_sd, 1)),
      sprintf("%s (%s)", format_number(treatment$survivor_mean, 1), format_number(treatment$survivor_sd, 1))
    ),
    c(
      "EQ-VAS, alive with observed score, median (IQR)",
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
      "EQ-VAS, alive with observed score, range",
      sprintf("%s to %s", format_number(control$survivor_min, 0), format_number(control$survivor_max, 0)),
      sprintf("%s to %s", format_number(treatment$survivor_min, 0), format_number(treatment$survivor_max, 0))
    ),
    c(
      "Coded endpoint mean",
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
    display_analysis <- if (identical(analysis, "Alive with EQ-VAS")) {
      "Alive with observed EQ-VAS"
    } else {
      analysis
    }
    display_method <- switch(
      method,
      "Welch t-test" = "Welch $t$-test",
      "Fisher exact test" = "Fisher's exact test",
      method
    )
    c(
      display_analysis,
      display_method,
      estimate,
      ifelse(is.na(row$statistic), "", format_number(row$statistic, 2)),
      format_p_value(row$p_value)
    )
  }
  interval_text <- function(interval, digits) {
    if (length(interval) == 2L && all(is.finite(interval))) {
      return(sprintf(
        "%s to %s",
        format_number(interval[[1]], digits),
        format_number(interval[[2]], digits)
      ))
    }
    "not estimated"
  }
  model_row <- function(label, fit) {
    c(
      label,
      "Joint likelihood ratio test",
      sprintf(
        "\\makecell[l]{$\\widehat{\\mu}_\\delta = %s$ (95\\%% CI %s)\\\\$\\widehat{\\alpha}_\\delta = %s$ (95\\%% CI %s)}",
        format_number(fit$mu_delta, 2),
        interval_text(fit$mu_delta_ci, 2),
        format_number(fit$alpha_delta, 3),
        interval_text(fit$alpha_delta_ci, 3)
      ),
      format_number(fit$statistic, 2),
      format_p_value(fit$p.value)
    )
  }

  delta_hat <- application_results$contrasts$estimate[application_results$contrasts$contrast == "Coded mean contrast Delta"]
  mu_hat <- application_results$contrasts$estimate[application_results$contrasts$contrast == "Alive-with-EQ-VAS mean difference"]

  rows <- rbind(
    std_row("Coded endpoint", "Welch t-test", sprintf("$\\widehat{\\Delta}^{(-1)} = %s$", format_number(delta_hat, 2))),
    std_row("Coded endpoint", "Wilcoxon rank-sum", "Death tied as worst value"),
    std_row("Alive with EQ-VAS", "Welch t-test", sprintf("$\\widehat{\\mu}_\\delta = %s$", format_number(mu_hat, 2))),
    std_row("Alive with EQ-VAS", "Wilcoxon rank-sum", "Conditional on being alive with observed EQ-VAS"),
    std_row("Death atom", "Fisher exact test", "Death proportions"),
    model_row("Parametric model", lrt),
    model_row("Semiparametric model", splrt)
  )

  write_latex_table(
    path = path,
    caption = "Standard and likelihood ratio analyses of the IST-3 endpoint.",
    label = "tab:ist3-frequentist",
    header = c("Analysis", "Method", "Estimate/target and 95\\% CI", "Statistic", "$p$-value"),
    rows = rows,
    align = "p{0.15\\textwidth}p{0.18\\textwidth}p{0.38\\textwidth}cc",
    size = "\\scriptsize",
    placement = "H"
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
      header = c("Quantity", "Posterior median", "95\\% CrI", "Posterior probability"),
      rows = rows,
      align = "llll",
      size = "\\small"
    ))
  }

  arm <- bayes$arm_table
  contrast <- bayes$summary_table
  probs <- bayes$probabilities
  interval <- function(tab, row) {
    sprintf("%s to %s", format_number(tab[row, "conf.low"], 2), format_number(tab[row, "conf.high"], 2))
  }

  rows <- rbind(
    c("Death probability, placebo/control", format_number(arm["rho_0", "estimate"], 3), interval(arm, "rho_0"), ""),
    c("Death probability, rt-PA", format_number(arm["rho_1", "estimate"], 3), interval(arm, "rho_1"), ""),
    c("Alive with observed EQ-VAS probability, placebo/control", format_number(arm["pi_0", "estimate"], 3), interval(arm, "pi_0"), ""),
    c("Alive with observed EQ-VAS probability, rt-PA", format_number(arm["pi_1", "estimate"], 3), interval(arm, "pi_1"), ""),
    c("Mean EQ-VAS, placebo/control", format_number(arm["mu_0_c", "estimate"], 2), interval(arm, "mu_0_c"), ""),
    c("Mean EQ-VAS, rt-PA", format_number(arm["mu_1_c", "estimate"], 2), interval(arm, "mu_1_c"), ""),
    c("$\\Delta_{\\mathrm{atom}}$", format_number(contrast["delta_atom", "estimate"], 3), interval(contrast, "delta_atom"), sprintf("$\\Prob(\\Delta_{\\mathrm{atom}}<0\\mid\\mathcal{D}_n) = %s$", format_probability(probs[["death_risk_reduction"]]))),
    c("$\\mu_\\delta$", format_number(contrast["mu_delta", "estimate"], 2), interval(contrast, "mu_delta"), sprintf("$\\Prob(\\mu_\\delta>0\\mid\\mathcal{D}_n) = %s$", format_probability(probs[["mu_delta_gt_0"]]))),
    c("$\\alpha_\\delta$", format_number(contrast["alpha_delta", "estimate"], 3), interval(contrast, "alpha_delta"), sprintf("$\\Prob(\\alpha_\\delta>1\\mid\\mathcal{D}_n) = %s$", format_probability(probs[["alpha_delta_gt_1"]]))),
    c("$\\Delta^{(-1)}$", format_number(contrast["delta", "estimate"], 2), interval(contrast, "delta"), sprintf("$\\Prob(\\Delta^{(-1)}>0\\mid\\mathcal{D}_n) = %s$", format_probability(probs[["delta_gt_0"]])))
  )

  write_latex_table(
    path = path,
    caption = "Bayesian IST-3 posterior summaries from the bounded reported score logit-normal two-part mixture model.",
    label = "tab:ist3-bayes",
    header = c("Quantity", "Posterior median", "95\\% CrI", "Posterior probability"),
    rows = rows,
    align = "p{0.34\\textwidth}p{0.14\\textwidth}p{0.17\\textwidth}p{0.17\\textwidth}",
    size = "\\scriptsize"
  )
}

write_application_tables <- function(application_results, tables_dir) {
  write_application_descriptive_table(application_results, file.path(tables_dir, "application-descriptives.tex"))
  write_application_frequentist_table(application_results, file.path(tables_dir, "application-frequentist.tex"))
  write_application_bayes_summary_table(application_results, file.path(tables_dir, "application-bayes-summary.tex"))
  invisible(TRUE)
}

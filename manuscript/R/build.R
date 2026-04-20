simulation_study_manuscript_results_path <- function(repo_root) {
  file.path(repo_root, "manuscript", "simulation-study-results.rds")
}

load_manuscript_simulation_results <- function(repo_root) {
  results_path <- simulation_study_manuscript_results_path(repo_root)
  if (!file.exists(results_path)) {
    stop(
      paste(
        "Missing manuscript/simulation-study-results.rds.",
        paste(
          "Run `Rscript simulation-study/scripts/collect-simulation-study-results.R`",
          "or `Rscript simulation-study/scripts/run-simulation-study.R` on the machine",
          "with the full raw simulation outputs to publish the manuscript-ready results file."
        )
      ),
      call. = FALSE
    )
  }

  simulation_results <- readRDS(results_path)
  if (!isTRUE(simulation_results$complete)) {
    stop(
      paste(
        "The published manuscript simulation results are incomplete.",
        "Regenerate the full study before building manuscript simulation figures and tables."
      ),
      call. = FALSE
    )
  }

  simulation_results
}

build_simulation_manuscript_assets <- function(repo_root) {
  study_dir <- file.path(repo_root, "simulation-study")

  source(file.path(study_dir, "R", "simulation-study.R"), local = globalenv())
  source(file.path(study_dir, "R", "manuscript-assets.R"), local = globalenv())

  simulation_results <- load_manuscript_simulation_results(repo_root)
  simulation_study_build_manuscript_assets(repo_root, simulation_results)

  invisible(simulation_results)
}

prune_manuscript_build_outputs <- function(output_dir) {
  figures_dir <- file.path(output_dir, "figures")
  tables_dir <- file.path(output_dir, "tables")

  expected_figure_files <- c(
    "application.pdf",
    "application-posterior.pdf",
    "application-ppc.pdf",
    "appendix-liver.pdf",
    "appendix-liver-posterior.pdf",
    "appendix-liver-ppc.pdf",
    "appendix-licorice.pdf",
    "appendix-licorice-posterior.pdf",
    "appendix-licorice-ppc.pdf",
    "power-curves.pdf",
    "power-effects.pdf",
    "type1-curves.pdf"
  )
  expected_table_files <- c(
    "application-descriptives.tex",
    "application-frequentist.tex",
    "application-delta-ci.tex",
    "application-bayes-summary.tex",
    "application-ppc.tex",
    "appendix-liver-descriptives.tex",
    "appendix-liver-frequentist.tex",
    "appendix-liver-delta-ci.tex",
    "appendix-liver-bayes-summary.tex",
    "appendix-liver-ppc.tex",
    "appendix-licorice-descriptives.tex",
    "appendix-licorice-frequentist.tex",
    "appendix-licorice-delta-ci.tex",
    "appendix-licorice-bayes-summary.tex",
    "appendix-licorice-ppc.tex",
    "simulation-scenarios.tex",
    sprintf("supplementary-power-s%d.tex", 1:6),
    "supplementary-type1.tex"
  )

  unlink(
    file.path(figures_dir, setdiff(list.files(figures_dir, full.names = FALSE), expected_figure_files)),
    force = TRUE
  )
  unlink(
    file.path(tables_dir, setdiff(list.files(tables_dir, full.names = FALSE), expected_table_files)),
    force = TRUE
  )
}

.application_test_p <- function(application_results, analysis, method) {
  tests <- application_results$standard_tests
  row <- tests[tests$analysis == analysis & tests$method == method, , drop = FALSE]
  if (!nrow(row)) {
    return(NA_real_)
  }
  row$p_value[[1]]
}

.format_macro_p_value <- function(x, digits = 3) {
  stopifnot(length(x) == 1, is.numeric(x), is.finite(x))
  if (x < 10^-digits) {
    return(sprintf("less than %s", formatC(10^-digits, format = "f", digits = digits)))
  }
  formatC(round(x, digits), format = "f", digits = digits)
}

.application_summary_row <- function(application_results, r) {
  application_results$summary_stats[application_results$summary_stats$R == r, , drop = FALSE]
}

.application_value_macros <- function(application_results) {
  metadata <- application_results$metadata
  control <- .application_summary_row(application_results, 0L)
  treatment <- .application_summary_row(application_results, 1L)
  lrt <- application_results$model_lrt
  splrt <- application_results$model_splrt
  bayes <- application_results$bayes

  bayes_probability <- function(name) {
    if (isTRUE(bayes$success)) {
      format_probability(bayes$probabilities[[name]])
    } else {
      "not estimated"
    }
  }

  c(
    list(
      ApplicationDatasetName = metadata$dataset_name,
      ApplicationDatasetDescription = metadata$dataset_description,
      ApplicationOutcomeLabel = metadata$outcome_label,
      ApplicationOutcomeShort = metadata$outcome_short,
      ApplicationAtomMeaning = metadata$atom_meaning,
      ApplicationObservedMeaning = metadata$observed_event_label,
      ApplicationGroupControl = metadata$group_labels[[1]],
      ApplicationGroupTreatment = metadata$group_labels[[2]],
      ApplicationN = format_count(nrow(application_results$data)),
      ApplicationRandomizedN = format_count(metadata$randomized_n),
      ApplicationControlN = format_count(control$n),
      ApplicationTreatmentN = format_count(treatment$n),
      ApplicationControlDeaths = format_count(control$atom_n),
      ApplicationTreatmentDeaths = format_count(treatment$atom_n),
      ApplicationControlDeathPct = format_percent(control$atom_prop),
      ApplicationTreatmentDeathPct = format_percent(treatment$atom_prop),
      ApplicationControlNonAtom = format_count(control$non_atom_n),
      ApplicationTreatmentNonAtom = format_count(treatment$non_atom_n),
      ApplicationControlMean = format_number(control$survivor_mean, 1),
      ApplicationTreatmentMean = format_number(treatment$survivor_mean, 1),
      ApplicationMuDelta = format_number(splrt$mu_delta, 2),
      ApplicationAlphaDelta = format_number(splrt$alpha_delta, 3),
      ApplicationDelta = format_number(splrt$delta, 2),
      ApplicationTTestP = format_p_value(.application_test_p(application_results, "Combined endpoint", "Welch t-test")),
      ApplicationWilcoxP = format_p_value(.application_test_p(application_results, "Combined endpoint", "Wilcoxon rank-sum")),
      ApplicationSurvivorTTestP = format_p_value(.application_test_p(application_results, "Alive with EQ-VAS", "Welch t-test")),
      ApplicationSurvivorWilcoxP = format_p_value(.application_test_p(application_results, "Alive with EQ-VAS", "Wilcoxon rank-sum")),
      ApplicationLRTP = format_p_value(lrt$p.value),
      ApplicationSPLRTP = format_p_value(splrt$p.value),
      ApplicationLRTStatistic = format_number(lrt$statistic, 2),
      ApplicationSPLRTStatistic = format_number(splrt$statistic, 2),
      ApplicationDeltaBenefitProbability = bayes_probability("delta_gt_0"),
      ApplicationMuBenefitProbability = bayes_probability("mu_delta_gt_0"),
      ApplicationAlphaBenefitProbability = bayes_probability("alpha_delta_gt_1"),
      ApplicationDeathBenefitProbability = bayes_probability("death_risk_reduction"),
      ApplicationBayesCache = "manuscript/application-data/ist3-bayes-cache.rds",
      ApplicationISTThreeDoi = metadata$doi
    )
  )
}

.liver_appendix_test_p <- function(liver_results, analysis, method) {
  tests <- liver_results$standard_tests
  row <- tests[tests$analysis == analysis & tests$method == method, , drop = FALSE]
  if (!nrow(row)) {
    return(NA_real_)
  }
  row$p_value[[1]]
}

.liver_appendix_summary_row <- function(liver_results, r) {
  liver_results$summary_stats[liver_results$summary_stats$R == r, , drop = FALSE]
}

.liver_appendix_value_macros <- function(liver_results) {
  metadata <- liver_results$metadata
  control <- .liver_appendix_summary_row(liver_results, 0L)
  treatment <- .liver_appendix_summary_row(liver_results, 1L)
  lrt <- liver_results$model_lrt
  splrt <- liver_results$model_splrt
  bayes <- liver_results$bayes

  bayes_probability <- function(name) {
    if (isTRUE(bayes$success)) {
      format_probability(bayes$probabilities[[name]])
    } else {
      "not estimated"
    }
  }

  c(
    list(
      LiverAppendixDatasetName = metadata$dataset_name,
      LiverAppendixDatasetDescription = metadata$dataset_description,
      LiverAppendixOutcomeLabel = metadata$outcome_label,
      LiverAppendixOutcomeShort = metadata$outcome_short,
      LiverAppendixAtomMeaning = metadata$atom_meaning,
      LiverAppendixObservedMeaning = metadata$observed_event_label,
      LiverAppendixGroupControl = metadata$group_labels[[1]],
      LiverAppendixGroupTreatment = metadata$group_labels[[2]],
      LiverAppendixPackageVersion = metadata$package_version,
      LiverAppendixN = format_count(nrow(liver_results$data)),
      LiverAppendixRandomizedN = format_count(metadata$randomized_n),
      LiverAppendixControlN = format_count(control$n),
      LiverAppendixTreatmentN = format_count(treatment$n),
      LiverAppendixControlDeaths = format_count(control$atom_n),
      LiverAppendixTreatmentDeaths = format_count(treatment$atom_n),
      LiverAppendixControlDeathPct = format_percent(control$atom_prop),
      LiverAppendixTreatmentDeathPct = format_percent(treatment$atom_prop),
      LiverAppendixControlNonAtom = format_count(control$non_atom_n),
      LiverAppendixTreatmentNonAtom = format_count(treatment$non_atom_n),
      LiverAppendixControlMean = format_number(control$survivor_mean, 1),
      LiverAppendixTreatmentMean = format_number(treatment$survivor_mean, 1),
      LiverAppendixControlCombinedMean = format_number(control$combined_mean, 1),
      LiverAppendixTreatmentCombinedMean = format_number(treatment$combined_mean, 1),
      LiverAppendixMuDelta = format_number(splrt$mu_delta, 2),
      LiverAppendixAlphaDelta = format_number(splrt$alpha_delta, 3),
      LiverAppendixDelta = format_number(splrt$delta, 2),
      LiverAppendixTTestP = format_p_value(.liver_appendix_test_p(liver_results, "Combined endpoint", "Welch t-test")),
      LiverAppendixWilcoxP = format_p_value(.liver_appendix_test_p(liver_results, "Combined endpoint", "Wilcoxon rank-sum")),
      LiverAppendixNonAtomTTestP = format_p_value(.liver_appendix_test_p(liver_results, "Non-atom component", "Welch t-test")),
      LiverAppendixNonAtomWilcoxP = format_p_value(.liver_appendix_test_p(liver_results, "Non-atom component", "Wilcoxon rank-sum")),
      LiverAppendixLRTP = format_p_value(lrt$p.value),
      LiverAppendixSPLRTP = format_p_value(splrt$p.value),
      LiverAppendixLRTStatistic = format_number(lrt$statistic, 2),
      LiverAppendixSPLRTStatistic = format_number(splrt$statistic, 2),
      LiverAppendixDeltaBenefitProbability = bayes_probability("delta_gt_0"),
      LiverAppendixMuBenefitProbability = bayes_probability("mu_delta_gt_0"),
      LiverAppendixAlphaBenefitProbability = bayes_probability("alpha_delta_gt_1"),
      LiverAppendixDeathBenefitProbability = bayes_probability("death_risk_reduction"),
      LiverAppendixExcludedCensored = format_count(metadata$excluded_censored_before_landmark),
      LiverAppendixExcludedControl = format_count(metadata$excluded_censored_before_landmark_by_group[[1]]),
      LiverAppendixExcludedTreatment = format_count(metadata$excluded_censored_before_landmark_by_group[[2]]),
      LiverAppendixBayesCache = "manuscript/application-data/liver-appendix-bayes-cache.rds"
    )
  )
}

.licorice_appendix_test_p <- function(licorice_results, analysis, method) {
  tests <- licorice_results$standard_tests
  row <- tests[tests$analysis == analysis & tests$method == method, , drop = FALSE]
  if (!nrow(row)) {
    return(NA_real_)
  }
  row$p_value[[1]]
}

.licorice_appendix_summary_row <- function(licorice_results, r) {
  licorice_results$summary_stats[licorice_results$summary_stats$R == r, , drop = FALSE]
}

.licorice_appendix_value_macros <- function(licorice_results) {
  metadata <- licorice_results$metadata
  control <- .licorice_appendix_summary_row(licorice_results, 0L)
  treatment <- .licorice_appendix_summary_row(licorice_results, 1L)
  lrt <- licorice_results$model_lrt
  splrt <- licorice_results$model_splrt
  bayes <- licorice_results$bayes

  bayes_probability <- function(name) {
    if (isTRUE(bayes$success)) {
      format_probability(bayes$probabilities[[name]])
    } else {
      "not estimated"
    }
  }

  c(
    list(
      LicoriceAppendixDatasetName = metadata$dataset_name,
      LicoriceAppendixDatasetDescription = metadata$dataset_description,
      LicoriceAppendixOutcomeLabel = metadata$outcome_label,
      LicoriceAppendixOutcomeShort = metadata$outcome_short,
      LicoriceAppendixAtomMeaning = metadata$atom_meaning,
      LicoriceAppendixObservedMeaning = metadata$observed_event_label,
      LicoriceAppendixGroupControl = metadata$group_labels[[1]],
      LicoriceAppendixGroupTreatment = metadata$group_labels[[2]],
      LicoriceAppendixPackageVersion = metadata$package_version,
      LicoriceAppendixPackageLicense = metadata$package_license,
      LicoriceAppendixN = format_count(nrow(licorice_results$data)),
      LicoriceAppendixRandomizedN = format_count(metadata$randomized_n),
      LicoriceAppendixControlN = format_count(control$n),
      LicoriceAppendixTreatmentN = format_count(treatment$n),
      LicoriceAppendixControlNoPain = format_count(control$atom_n),
      LicoriceAppendixTreatmentNoPain = format_count(treatment$atom_n),
      LicoriceAppendixControlNoPainPct = format_percent(control$atom_prop),
      LicoriceAppendixTreatmentNoPainPct = format_percent(treatment$atom_prop),
      LicoriceAppendixControlPain = format_count(control$non_atom_n),
      LicoriceAppendixTreatmentPain = format_count(treatment$non_atom_n),
      LicoriceAppendixControlPainPct = format_percent(control$non_atom_prop),
      LicoriceAppendixTreatmentPainPct = format_percent(treatment$non_atom_prop),
      LicoriceAppendixControlMean = format_number(control$survivor_mean, 2),
      LicoriceAppendixTreatmentMean = format_number(treatment$survivor_mean, 2),
      LicoriceAppendixControlCombinedMean = format_number(control$combined_mean, 2),
      LicoriceAppendixTreatmentCombinedMean = format_number(treatment$combined_mean, 2),
      LicoriceAppendixMuDelta = format_number(splrt$mu_delta, 2),
      LicoriceAppendixAlphaDelta = format_number(splrt$alpha_delta, 3),
      LicoriceAppendixDelta = format_number(splrt$delta, 2),
      LicoriceAppendixTTestP = .format_macro_p_value(.licorice_appendix_test_p(licorice_results, "Combined endpoint", "Welch t-test")),
      LicoriceAppendixWilcoxP = .format_macro_p_value(.licorice_appendix_test_p(licorice_results, "Combined endpoint", "Wilcoxon rank-sum")),
      LicoriceAppendixNonAtomTTestP = .format_macro_p_value(.licorice_appendix_test_p(licorice_results, "Positive pain component", "Welch t-test")),
      LicoriceAppendixNonAtomWilcoxP = .format_macro_p_value(.licorice_appendix_test_p(licorice_results, "Positive pain component", "Wilcoxon rank-sum")),
      LicoriceAppendixAtomFisherP = .format_macro_p_value(.licorice_appendix_test_p(licorice_results, "No-pain atom", "Fisher exact test")),
      LicoriceAppendixLRTP = .format_macro_p_value(lrt$p.value),
      LicoriceAppendixSPLRTP = .format_macro_p_value(splrt$p.value),
      LicoriceAppendixLRTStatistic = format_number(lrt$statistic, 2),
      LicoriceAppendixSPLRTStatistic = format_number(splrt$statistic, 2),
      LicoriceAppendixDeltaBenefitProbability = bayes_probability("delta_lt_0"),
      LicoriceAppendixMuBenefitProbability = bayes_probability("mu_delta_lt_0"),
      LicoriceAppendixAlphaBenefitProbability = bayes_probability("alpha_delta_lt_1"),
      LicoriceAppendixNoPainBenefitProbability = bayes_probability("no_pain_probability_increase"),
      LicoriceAppendixAnyPainReductionProbability = bayes_probability("any_pain_probability_reduction"),
      LicoriceAppendixExcludedMissing = format_count(metadata$excluded_missing_outcome),
      LicoriceAppendixExcludedControl = format_count(metadata$excluded_missing_outcome_by_group[[1]]),
      LicoriceAppendixExcludedTreatment = format_count(metadata$excluded_missing_outcome_by_group[[2]]),
      LicoriceAppendixBayesCache = "manuscript/application-data/licorice-appendix-bayes-cache.rds"
    )
  )
}

build_manuscript_assets <- function(output_dir, repo_root) {
  manuscript_dir <- file.path(repo_root, "manuscript")
  output_dir <- ensure_dir(output_dir)
  figures_dir <- ensure_dir(file.path(output_dir, "figures"))
  tables_dir <- ensure_dir(file.path(output_dir, "tables"))

  load_local_trunccomp2(repo_root)
  simulation_results <- build_simulation_manuscript_assets(repo_root)
  application_results <- compute_application_results(manuscript_dir)
  liver_results <- compute_liver_appendix_results(manuscript_dir)
  licorice_results <- compute_licorice_appendix_results(manuscript_dir)
  prune_manuscript_build_outputs(output_dir)

  build_application_figure(application_results, file.path(figures_dir, "application.pdf"))
  build_application_posterior_figure(application_results, file.path(figures_dir, "application-posterior.pdf"))
  build_application_ppc_figure(application_results, file.path(figures_dir, "application-ppc.pdf"))
  write_application_tables(application_results, tables_dir)
  build_liver_appendix_figure(liver_results, file.path(figures_dir, "appendix-liver.pdf"))
  build_liver_appendix_posterior_figure(liver_results, file.path(figures_dir, "appendix-liver-posterior.pdf"))
  build_liver_appendix_ppc_figure(liver_results, file.path(figures_dir, "appendix-liver-ppc.pdf"))
  write_liver_appendix_tables(liver_results, tables_dir)
  build_licorice_appendix_figure(licorice_results, file.path(figures_dir, "appendix-licorice.pdf"))
  build_licorice_appendix_posterior_figure(licorice_results, file.path(figures_dir, "appendix-licorice-posterior.pdf"))
  build_licorice_appendix_ppc_figure(licorice_results, file.path(figures_dir, "appendix-licorice-ppc.pdf"))
  write_licorice_appendix_tables(licorice_results, tables_dir)

  write_generated_values_tex(
    file.path(output_dir, "generated-values.tex"),
    c(
      .application_value_macros(application_results),
      .liver_appendix_value_macros(liver_results),
      .licorice_appendix_value_macros(licorice_results)
    )
  )
  unlink(file.path(manuscript_dir, "Rplots.pdf"), force = TRUE)

  invisible(list(
    simulation_results = simulation_results,
    application_results = application_results,
    liver_results = liver_results,
    licorice_results = licorice_results
  ))
}

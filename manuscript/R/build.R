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
      ApplicationSurvivorTTestP = format_p_value(.application_test_p(application_results, "Survivors only", "Welch t-test")),
      ApplicationSurvivorWilcoxP = format_p_value(.application_test_p(application_results, "Survivors only", "Wilcoxon rank-sum")),
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

build_manuscript_assets <- function(output_dir, repo_root) {
  manuscript_dir <- file.path(repo_root, "manuscript")
  output_dir <- ensure_dir(output_dir)
  figures_dir <- ensure_dir(file.path(output_dir, "figures"))
  tables_dir <- ensure_dir(file.path(output_dir, "tables"))

  load_local_trunccomp2(repo_root)
  simulation_results <- build_simulation_manuscript_assets(repo_root)
  application_results <- compute_application_results(manuscript_dir)
  prune_manuscript_build_outputs(output_dir)

  build_application_figure(application_results, file.path(figures_dir, "application.pdf"))
  build_application_posterior_figure(application_results, file.path(figures_dir, "application-posterior.pdf"))
  build_application_ppc_figure(application_results, file.path(figures_dir, "application-ppc.pdf"))
  write_application_tables(application_results, tables_dir)

  write_generated_values_tex(
    file.path(output_dir, "generated-values.tex"),
    .application_value_macros(application_results)
  )
  unlink(file.path(manuscript_dir, "Rplots.pdf"), force = TRUE)

  invisible(list(
    simulation_results = simulation_results,
    application_results = application_results
  ))
}

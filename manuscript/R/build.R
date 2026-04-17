assert_simulation_assets_exist <- function(output_dir) {
  required_paths <- c(
    file.path(output_dir, "figures", "power-curves.pdf"),
    file.path(output_dir, "figures", "power-effects.pdf"),
    file.path(output_dir, "figures", "supplementary-power-curves.pdf"),
    file.path(output_dir, "figures", "type1-curves.pdf"),
    file.path(output_dir, "tables", "simulation-scenarios.tex"),
    file.path(output_dir, "tables", "supplementary-power-s1.tex"),
    file.path(output_dir, "tables", "supplementary-power-s2.tex"),
    file.path(output_dir, "tables", "supplementary-power-s3.tex"),
    file.path(output_dir, "tables", "supplementary-power-s4.tex"),
    file.path(output_dir, "tables", "supplementary-power-s5.tex"),
    file.path(output_dir, "tables", "supplementary-power-s6.tex"),
    file.path(output_dir, "tables", "supplementary-type1.tex")
  )

  missing <- required_paths[!file.exists(required_paths)]
  if (length(missing)) {
    stop(
      paste(
        "Missing prebuilt simulation-study manuscript assets under manuscript/build.",
        "Run `Rscript simulation-study/scripts/run-simulation-study.R` to generate them."
      ),
      call. = FALSE
    )
  }

  invisible(required_paths)
}

build_manuscript_assets <- function(output_dir, repo_root) {
  manuscript_dir <- file.path(repo_root, "manuscript")
  output_dir <- ensure_dir(output_dir)
  figures_dir <- ensure_dir(file.path(output_dir, "figures"))
  ensure_dir(file.path(output_dir, "tables"))

  load_local_trunccomp2(repo_root)

  example_results <- compute_example_results(repo_root)
  application_results <- compute_application_results(manuscript_dir)
  assert_simulation_assets_exist(output_dir)

  build_example_histogram(example_results, file.path(figures_dir, "example-histogram.pdf"))
  build_example_surface_figure(example_results, file.path(figures_dir, "example-simultaneous-confidence.pdf"))
  build_application_figure(application_results, file.path(figures_dir, "application.pdf"))

  values <- list(
    ExampleDelta = format_fixed(example_results$delta, digits = 3),
    ExampleTTestP = format_fixed(example_results$ttest_p, digits = 3),
    ExampleWilcoxP = format_fixed(example_results$wilcox_p, digits = 3),
    ApplicationDatasetName = application_results$metadata$dataset_name,
    ApplicationDatasetDescription = application_results$metadata$dataset_description,
    ApplicationOutcomeLabel = application_results$metadata$outcome_label,
    ApplicationOutcomeShort = application_results$metadata$outcome_short,
    ApplicationAtomMeaning = application_results$metadata$atom_meaning,
    ApplicationObservedMeaning = application_results$metadata$observed_event_label,
    ApplicationGroupControl = application_results$metadata$group_labels[[1]],
    ApplicationGroupTreatment = application_results$metadata$group_labels[[2]],
    ApplicationN = format_count(nrow(application_results$data)),
    ApplicationTTestP = format_fixed(application_results$ttest_p, digits = 3),
    ApplicationWilcoxP = format_fixed(application_results$wilcox_p, digits = 3),
    ApplicationLRTP = format_fixed(application_results$model_lrt$p, digits = 3),
    ApplicationSPLRTP = format_fixed(application_results$model_splrt$p, digits = 3),
    ApplicationMuDelta = format_number(application_results$model_splrt$muDelta, digits = 1),
    ApplicationAlphaDelta = format_fixed(unname(as.numeric(application_results$model_splrt$alphaDelta)), digits = 3),
    ApplicationDelta = format_number(application_results$model_splrt$Delta, digits = 1),
    ApplicationConclusionSummary = application_results$conclusion_summary,
    ApplicationComponentSummary = application_results$component_summary
  )
  write_generated_values_tex(file.path(output_dir, "generated-values.tex"), values)
  unlink(file.path(manuscript_dir, "Rplots.pdf"), force = TRUE)

  invisible(list(
    example_results = example_results,
    application_results = application_results
  ))
}

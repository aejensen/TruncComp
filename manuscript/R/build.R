build_manuscript_assets <- function(output_dir, repo_root) {
  manuscript_dir <- file.path(repo_root, "manuscript")
  output_dir <- ensure_dir(output_dir)
  figures_dir <- ensure_dir(file.path(output_dir, "figures"))
  tables_dir <- ensure_dir(file.path(output_dir, "tables"))

  load_local_trunccomp(repo_root)

  example_results <- compute_example_results(repo_root)
  application_results <- compute_application_results(manuscript_dir)
  simulation_results <- load_simulation_results(manuscript_dir)

  build_example_histogram(example_results, file.path(figures_dir, "example-histogram.pdf"))
  build_example_surface_figure(example_results, file.path(figures_dir, "example-simultaneous-confidence.pdf"))
  build_power_curves_figure(simulation_results, file.path(figures_dir, "power-curves.pdf"))
  build_application_figure(application_results, file.path(figures_dir, "application.pdf"))
  build_supplementary_tables(simulation_results, tables_dir)

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
    application_results = application_results,
    simulation_results = simulation_results
  ))
}

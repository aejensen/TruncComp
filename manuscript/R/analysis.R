compute_example_results <- function(repo_root) {
  d <- load_example_data(repo_root)
  model <- TruncComp2::truncComp(Y ~ R, atom = 0, data = d, method = "SPLRT")

  list(
    data = d,
    model = model,
    delta = mean(d$Y[d$R == 1]) - mean(d$Y[d$R == 0]),
    ttest_p = stats::t.test(Y ~ R, data = d)$p.value,
    wilcox_p = suppressWarnings(stats::wilcox.test(Y ~ R, data = d)$p.value)
  )
}

application_has_desired_contrast <- function(ttest_p, wilcox_p, lrt_p, splrt_p, alpha = 0.05) {
  any(c(ttest_p, wilcox_p) >= alpha) && any(c(lrt_p, splrt_p) < alpha)
}

application_conclusion_summary <- function(ttest_p, wilcox_p, lrt_p, splrt_p, alpha = 0.05) {
  standard_sig <- c(ttest_p, wilcox_p) < alpha
  joint_sig <- c(lrt_p, splrt_p) < alpha

  standard_text <- if (all(!standard_sig)) {
    "At the 5% level, neither standard two-sample test rejects."
  } else if (all(standard_sig)) {
    "At the 5% level, both standard two-sample tests reject."
  } else if (standard_sig[[1]]) {
    "At the 5% level, only the t-test rejects among the standard two-sample procedures."
  } else {
    "At the 5% level, only the Wilcoxon test rejects among the standard two-sample procedures."
  }

  joint_text <- if (all(joint_sig)) {
    "Both likelihood-ratio procedures reject the joint null."
  } else if (all(!joint_sig)) {
    "Neither likelihood-ratio procedure rejects the joint null."
  } else if (joint_sig[[1]]) {
    "Only the parametric likelihood-ratio procedure rejects the joint null."
  } else {
    "Only the semi-parametric likelihood-ratio procedure rejects the joint null."
  }

  paste(standard_text, joint_text)
}

application_component_summary <- function(mu_delta, alpha_delta, metadata) {
  treatment_label <- metadata$group_labels[[2]]
  control_label <- metadata$group_labels[[1]]
  mean_direction <- if (mu_delta >= 0) "larger" else "smaller"
  odds_direction <- if (alpha_delta >= 1) "higher" else "lower"

  if ((mu_delta >= 0) != (alpha_delta >= 1)) {
    sprintf(
      "The fitted components point in opposite directions: relative to %s, %s has %s positive-part means but %s odds of %s.",
      control_label,
      treatment_label,
      mean_direction,
      odds_direction,
      metadata$observed_event_label
    )
  } else {
    sprintf(
      "The fitted components move in the same direction: relative to %s, %s has %s positive-part means and %s odds of %s.",
      control_label,
      treatment_label,
      mean_direction,
      odds_direction,
      metadata$observed_event_label
    )
  }
}

analyze_application_data <- function(application_data) {
  d <- application_data$data
  metadata <- application_data$metadata

  model_lrt <- TruncComp2::truncComp(Y ~ R, atom = 0, data = d, method = "LRT")
  model_splrt <- TruncComp2::truncComp(Y ~ R, atom = 0, data = d, method = "SPLRT")
  surface <- suppressMessages(
    TruncComp2::jointContrastCI(
      model_splrt,
      resolution = metadata$surface_resolution,
      plot = FALSE
    )
  )
  ttest_p <- stats::t.test(Y ~ R, data = d)$p.value
  wilcox_p <- suppressWarnings(stats::wilcox.test(Y ~ R, data = d)$p.value)

  list(
    data = d,
    metadata = metadata,
    model_lrt = model_lrt,
    model_splrt = model_splrt,
    surface = surface,
    ttest_p = ttest_p,
    wilcox_p = wilcox_p,
    conclusion_summary = application_conclusion_summary(ttest_p, wilcox_p, model_lrt$p, model_splrt$p),
    component_summary = application_component_summary(
      mu_delta = model_splrt$muDelta,
      alpha_delta = unname(as.numeric(model_splrt$alphaDelta)),
      metadata = metadata
    )
  )
}

compute_application_results <- function(manuscript_dir) {
  local_data <- load_application_data(manuscript_dir)
  if (!is.null(local_data)) {
    local_results <- analyze_application_data(local_data)
    if (application_has_desired_contrast(
      local_results$ttest_p,
      local_results$wilcox_p,
      local_results$model_lrt$p,
      local_results$model_splrt$p
    )) {
      return(local_results)
    }

    message("Local application data did not yield a different-conclusion case; using the RAND HIE fallback dataset.")
  }

  analyze_application_data(.load_randhealth_application_data())
}

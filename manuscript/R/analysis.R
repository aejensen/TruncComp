compute_example_results <- function(repo_root) {
  d <- load_example_data(repo_root)
  model <- TruncComp::truncComp(Y ~ R, atom = 0, data = d, method = "SPLRT")

  list(
    data = d,
    model = model,
    delta = mean(d$Y[d$R == 1]) - mean(d$Y[d$R == 0]),
    ttest_p = stats::t.test(Y ~ R, data = d)$p.value,
    wilcox_p = suppressWarnings(stats::wilcox.test(Y ~ R, data = d)$p.value)
  )
}

compute_application_results <- function(manuscript_dir) {
  d <- load_application_data(manuscript_dir)
  model_lrt <- TruncComp::truncComp(dawols28 ~ allocation, atom = 0, data = d, method = "LRT")
  model_splrt <- TruncComp::truncComp(dawols28 ~ allocation, atom = 0, data = d, method = "SPLRT")
  surface <- suppressMessages(stats::confint(
    model_splrt,
    type = "simultaneous",
    offset = 2,
    resolution = 100,
    plot = FALSE
  ))

  list(
    data = d,
    model_lrt = model_lrt,
    model_splrt = model_splrt,
    surface = surface,
    wilcox_p = suppressWarnings(stats::wilcox.test(dawols28 ~ allocation, data = d)$p.value)
  )
}

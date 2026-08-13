test_that("internal Stan model overrides require an explicit model object", {
  d <- data.frame(
    Y = c(-1, 20, -1, 30),
    A = c(0L, 1L, 0L, 1L),
    R = c(0L, 0L, 1L, 1L)
  )
  arguments <- list(
    data = d,
    atom = -1,
    continuous_support = "bounded_score",
    bounded_kernel = "logit_normal",
    bounded_kernel_supplied = TRUE,
    score_min = 0,
    score_max = 100,
    score_step = 1,
    heaping_grids = c(1, 5, 10),
    heaping = "shared",
    support_supplied = bayes_support_supplied_defaults(
      score_min = TRUE,
      score_max = TRUE,
      score_step = TRUE,
      heaping_grids = TRUE,
      heaping = TRUE
    ),
    mixture_components = 2L,
    auto_select_mixture_components = FALSE,
    chains = 1L,
    iter_warmup = 1L,
    iter_sampling = 1L,
    prior = list(eta_prior = c(2, 2, 2))
  )

  expect_error(
    do.call(
      fit_trunc_comp_bayes,
      c(arguments, list(model_name_override = "application_model"))
    ),
    "requires an explicit model_object",
    fixed = TRUE
  )
  expect_error(
    do.call(
      fit_trunc_comp_bayes,
      c(
        arguments,
        list(model_object = structure(list(), class = "test_model"), model_name_override = "")
      )
    ),
    "single non-empty character string",
    fixed = TRUE
  )
})

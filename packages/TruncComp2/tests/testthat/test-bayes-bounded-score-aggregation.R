bounded_score_aggregation_fixture <- function(heaping = "arm_specific") {
  data <- data.frame(
    Y = c(-1, 0, 0, 2, 4, -1, 1, 1, 3, 4, 4),
    A = c(0L, 1L, 1L, 1L, 1L, 0L, 1L, 1L, 1L, 1L, 1L),
    R = c(0L, 0L, 0L, 0L, 0L, 1L, 1L, 1L, 1L, 1L, 1L)
  )
  support_options <- TruncComp2:::bayes_normalize_support_options(
    continuous_support = "bounded_score",
    score_min = 0,
    score_max = 4,
    score_step = 1,
    heaping_grids = c(1, 2),
    heaping = heaping
  )
  prior <- TruncComp2:::normalize_bayes_prior(
    NULL,
    continuous_support = "bounded_score"
  )
  observed <- data[data$A == 1L, , drop = FALSE]
  y_obs_index <- TruncComp2:::bayes_score_grid_index(
    observed$Y,
    score_min = support_options$score_min,
    score_max = support_options$score_max,
    score_step = support_options$score_step
  )

  list(
    data = data,
    support_options = support_options,
    arm_obs = as.integer(observed$R) + 1L,
    y_obs_index = y_obs_index,
    standata = TruncComp2:::build_bayes_standata(
      data = data,
      atom = -1,
      mixture_components = 2,
      prior = prior,
      continuous_support = "bounded_score",
      support_options = support_options
    )$stan_data
  )
}

test_that("bounded-score Stan data aggregates non-atom scores by arm and score", {
  fixture <- bounded_score_aggregation_fixture()
  standata <- fixture$standata

  expect_equal(standata$N_score_cells, 6L)
  expect_equal(standata$score_cell_arm, c(1L, 1L, 1L, 2L, 2L, 2L))
  expect_equal(standata$score_cell_index, c(1L, 3L, 5L, 2L, 4L, 5L))
  expect_equal(standata$score_cell_count, c(2L, 1L, 1L, 2L, 1L, 2L))
  expect_equal(sum(standata$score_cell_count), 9L)
  expect_false(any(standata$score_cell_arm == 1L & standata$score_cell_index == 2L))
  expect_false(any(standata$score_cell_arm == 2L & standata$score_cell_index == 1L))

  expect_equal(standata$n_arm, c(5L, 6L))
  expect_equal(standata$n_obs_arm, c(4L, 5L))
  expect_null(standata[["N_obs"]])
  expect_null(standata[["arm_obs"]])
  expect_null(standata[["y_obs_index"]])
  expect_null(standata[["A"]])
  expect_null(standata[["arm"]])
})

test_that("bounded-score aggregated likelihood matches observation-level arithmetic", {
  for(heaping in c("shared", "arm_specific")) {
    fixture <- bounded_score_aggregation_fixture(heaping = heaping)
    standata <- fixture$standata
    weights <- rbind(c(0.65, 0.35), c(0.40, 0.60))
    m_comp <- rbind(c(0.25, 0.72), c(0.38, 0.82))
    phi_comp <- rbind(c(7, 12), c(8, 10))
    eta <- if(standata$eta_groups == 1L) {
      matrix(c(0.65, 0.35), nrow = 1L)
    } else {
      rbind(c(0.65, 0.35), c(0.30, 0.70))
    }

    pmf_by_arm <- lapply(1:2, function(arm) {
      eta_group <- standata$eta_group_by_arm[[arm]]
      TruncComp2:::bayes_score_pmf(
        weights = weights[arm, ],
        m_comp = m_comp[arm, ],
        phi_comp = phi_comp[arm, ],
        eta = eta[eta_group, ],
        bin_lower = standata$bin_lower,
        bin_upper = standata$bin_upper,
        bin_valid = standata$bin_valid
      )
    })

    observation_lp <- sum(vapply(
      seq_along(fixture$y_obs_index),
      function(i) log(pmf_by_arm[[fixture$arm_obs[[i]]]][[fixture$y_obs_index[[i]]]]),
      numeric(1)
    ))
    cell_lp <- sum(vapply(
      seq_len(standata$N_score_cells),
      function(cell) {
        standata$score_cell_count[[cell]] *
          log(pmf_by_arm[[standata$score_cell_arm[[cell]]]][[standata$score_cell_index[[cell]]]])
      },
      numeric(1)
    ))

    expect_equal(cell_lp, observation_lp, tolerance = 1e-12)

    pi <- c(0.70, 0.45)
    arm_full <- as.integer(fixture$data$R) + 1L
    atom_observation_lp <- sum(log(ifelse(fixture$data$A == 1L, pi[arm_full], 1 - pi[arm_full])))
    atom_count_lp <- sum(
      standata$n_obs_arm * log(pi) +
        (standata$n_arm - standata$n_obs_arm) * log1p(-pi)
    )

    expect_equal(atom_count_lp, atom_observation_lp, tolerance = 1e-12)
  }
})

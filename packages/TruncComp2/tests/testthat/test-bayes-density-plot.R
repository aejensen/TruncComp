bayes_density_plot_fit <- bayes_formula_fit(seed = 505)
bayes_positive_density_plot_fit <- bayes_positive_formula_fit(seed = 1505)
bayes_bounded_continuous_density_plot_fit <- bayes_bounded_continuous_formula_fit(seed = 2505)
bayes_bounded_score_density_plot_fit <- bayes_bounded_score_formula_fit(seed = 2506)

manual_bayes_density_matrix <- function(fit, x, arm_index) {
  support <- fit$settings$continuous_support

  if(identical(support, "positive_real")) {
    extracted <- rstan::extract(
      fit$fit,
      pars = c("w", "mean_comp", "shape_comp", "rho"),
      permuted = TRUE,
      inc_warmup = FALSE
    )

    draws <- dim(extracted$w)[1]

    return(t(vapply(
      seq_len(draws),
      function(draw) {
        weights <- extracted$w[draw, arm_index, ]
        means <- fit$settings$y_scale * extracted$mean_comp[draw, arm_index, ]
        shapes <- extracted$shape_comp[draw, arm_index, ]
        pi <- 1 - extracted$rho[draw, arm_index]

        kernels <- vapply(
          seq_along(weights),
          function(h) TruncComp2:::bayes_positive_kernel_density(x, mean = means[h], shape = shapes[h]),
          numeric(length(x))
        )

        as.numeric(pi * (kernels %*% weights))
      },
      numeric(length(x))
    )))
  }

  extracted <- rstan::extract(
    fit$fit,
    pars = c("w", "mu_comp", "sigma_comp", "rho"),
    permuted = TRUE,
    inc_warmup = FALSE
  )

  center <- fit$settings$y_center
  scale <- fit$settings$y_scale
  draws <- dim(extracted$w)[1]

  t(vapply(
    seq_len(draws),
    function(draw) {
      weights <- extracted$w[draw, arm_index, ]
      means <- center + scale * extracted$mu_comp[draw, arm_index, ]
      sds <- scale * extracted$sigma_comp[draw, arm_index, ]
      pi <- 1 - extracted$rho[draw, arm_index]

      kernels <- vapply(
        seq_along(weights),
        function(h) stats::dnorm(x, mean = means[h], sd = sds[h]),
        numeric(length(x))
      )

      as.numeric(pi * (kernels %*% weights))
    },
    numeric(length(x))
  ))
}

test_that("posterior_density_plot returns a faceted ggplot with the expected layers", {
  plot <- posterior_density_plot(bayes_density_plot_fit)

  expect_s3_class(plot, "ggplot")
  geoms <- vapply(plot$layers, function(layer) class(layer$geom)[1], character(1))
  expect_true(all(c("GeomRibbon", "GeomLine", "GeomSegment", "GeomText") %in% geoms))

  built <- ggplot2::ggplot_build(plot)
  expect_equal(nrow(built$layout$layout), 2)
})

test_that("Bayesian density helper respects x, summarizes densities, and responds to conf.level", {
  x_grid <- seq(0.25, 2.25, length.out = 7)
  data_default <- TruncComp2:::bayes_density_plot_data(
    bayes_density_plot_fit,
    x = NULL,
    n = 120
  )
  data_80 <- TruncComp2:::bayes_density_plot_data(
    bayes_density_plot_fit,
    x = x_grid,
    conf.level = 0.80
  )
  data_95 <- TruncComp2:::bayes_density_plot_data(
    bayes_density_plot_fit,
    x = x_grid,
    conf.level = 0.95
  )

  density_80 <- data_80$density_data
  density_95 <- data_95$density_data

  expect_equal(split(density_80$x, density_80$arm_label), list(Control = x_grid, Treatment = x_grid))
  expect_equal(unname(as.integer(table(density_80$arm_label))), c(length(x_grid), length(x_grid)))

  widths_80 <- density_80$conf.high - density_80$conf.low
  widths_95 <- density_95$conf.high - density_95$conf.low
  expect_true(all(widths_95 >= widths_80 - 1e-10))

  manual_control <- manual_bayes_density_matrix(bayes_density_plot_fit, x_grid, arm_index = 1)
  manual_treatment <- manual_bayes_density_matrix(bayes_density_plot_fit, x_grid, arm_index = 2)

  expect_equal(
    density_80$density_mean[density_80$arm_label == "Control"],
    colMeans(manual_control),
    tolerance = 1e-10
  )
  expect_equal(
    density_80$density_mean[density_80$arm_label == "Treatment"],
    colMeans(manual_treatment),
    tolerance = 1e-10
  )

  expect_equal(data_80$y_limits[1], 0)
  expect_true(all(data_80$atom_data$y_end >= data_80$atom_data$atom_mass))
  expect_true(all(data_80$atom_data$y_end < data_80$y_limits[2]))
  expect_true(data_80$y_limits[2] > max(density_80$conf.high))
  expect_true(all(data_80$atom_data$x_label > data_80$atom_data$atom))
  expect_true(all(data_80$atom_data$x_label < data_80$x_limits[2]))
  expect_true(all(data_80$atom_data$y_label > data_80$atom_data$y_end))
  expect_true(all(data_80$atom_data$y_label < data_80$y_limits[2]))

  observed_range <- range(bayes_density_plot_fit$data$Y[bayes_density_plot_fit$data$A == 1])
  expect_lte(data_default$x_limits[1], observed_range[1])
  expect_gte(data_default$x_limits[2], observed_range[2])
  expect_true(diff(data_default$x_limits) > diff(observed_range))

  edge_density <- vapply(
    split(data_default$density_data, data_default$density_data$arm_label),
    function(df) {
      peak <- max(df$density_mean)
      max(df$density_mean[c(1, nrow(df))] / peak)
    },
    numeric(1)
  )
  expect_true(all(edge_density < 0.4))
})

test_that("posterior_density_plot validates failed fits and plotting inputs", {
  failed_fit <- bayes_density_plot_fit
  failed_fit$success <- FALSE
  expect_error(posterior_density_plot(failed_fit), "Estimation failed")

  missing_fit <- bayes_density_plot_fit
  missing_fit$fit <- NULL
  expect_error(posterior_density_plot(missing_fit), "Raw Stan mixture parameters")

  expect_error(posterior_density_plot(bayes_density_plot_fit, x = 1), "x must be NULL or a finite numeric grid")
  expect_error(posterior_density_plot(bayes_density_plot_fit, x = c(1, NA)), "x must be NULL or a finite numeric grid")
  expect_error(posterior_density_plot(bayes_density_plot_fit, n = 1), "n must be a single integer >= 2")
  expect_error(posterior_density_plot(bayes_density_plot_fit, conf.level = 1), "conf.level must be")
})

test_that("Bayesian density helper builds a usable default grid for degenerate observed ranges", {
  degenerate_fit <- bayes_density_plot_fit
  degenerate_fit$data$Y[degenerate_fit$data$A == 1] <- 1.5
  degenerate_fit$settings$y_scale <- 1

  plot_data <- TruncComp2:::bayes_density_plot_data(
    degenerate_fit,
    x = NULL,
    n = 5
  )

  expect_equal(length(plot_data$x), 5)
  expect_true(all(is.finite(plot_data$x)))
  expect_true(diff(range(plot_data$x)) > 0)
})

test_that("posterior_density_plot works for positive-support fits", {
  x_grid <- seq(0, 4, length.out = 9)
  plot <- posterior_density_plot(bayes_positive_density_plot_fit)
  default_data <- TruncComp2:::bayes_density_plot_data(
    bayes_positive_density_plot_fit,
    x = NULL,
    n = 80
  )
  plot_data <- TruncComp2:::bayes_density_plot_data(
    bayes_positive_density_plot_fit,
    x = x_grid,
    conf.level = 0.9
  )
  manual_x <- x_grid[x_grid > 0]
  manual_control <- manual_bayes_density_matrix(bayes_positive_density_plot_fit, manual_x, arm_index = 1)
  manual_treatment <- manual_bayes_density_matrix(bayes_positive_density_plot_fit, manual_x, arm_index = 2)

  expect_s3_class(plot, "ggplot")
  expect_true(all(is.finite(plot_data$density_data$density_mean)))
  expect_true(all(plot_data$density_data$x > 0))
  expect_equal(plot_data$x_limits[1], 0)
  expect_true(default_data$x_limits[2] >= max(default_data$x))
  expect_true(default_data$x_limits[1] <= bayes_positive_density_plot_fit$atom)
  expect_equal(
    plot_data$density_data$density_mean[plot_data$density_data$arm_label == "Control"],
    colMeans(manual_control),
    tolerance = 1e-10
  )
  expect_equal(
    plot_data$density_data$density_mean[plot_data$density_data$arm_label == "Treatment"],
    colMeans(manual_treatment),
    tolerance = 1e-10
  )
})

test_that("posterior_density_plot works for bounded-continuous fits", {
  plot <- posterior_density_plot(bayes_bounded_continuous_density_plot_fit)
  plot_data <- TruncComp2:::bayes_density_plot_data(
    bayes_bounded_continuous_density_plot_fit,
    x = NULL,
    n = 50
  )

  expect_s3_class(plot, "ggplot")
  expect_equal(plot_data$plot_type, "density")
  expect_true(all(plot_data$x > 0 & plot_data$x < 100))
  expect_equal(plot_data$x_limits, c(0, 100))
  expect_true(all(is.finite(plot_data$density_data$density_mean)))
})

test_that("posterior_density_plot works for bounded-score mass plots", {
  plot <- posterior_density_plot(bayes_bounded_score_density_plot_fit)
  plot_data <- TruncComp2:::bayes_density_plot_data(
    bayes_bounded_score_density_plot_fit,
    n = 50
  )

  expect_s3_class(plot, "ggplot")
  expect_equal(plot_data$plot_type, "mass")
  expect_equal(plot$labels$x, "Reported score")
  expect_equal(plot$labels$y, "Posterior predictive probability")
  expect_true(all(plot_data$density_data$x %in% bayes_bounded_score_density_plot_fit$settings$score_values))
  expect_true(all(plot_data$density_data$density_mean >= 0))
  expect_error(
    posterior_density_plot(bayes_bounded_score_density_plot_fit, x = c(0, 100)),
    "x is not used"
  )
})

test_that("Bayesian example vignette renders with posterior density plot", {
  skip_on_cran()
  skip_if_not_installed("rmarkdown")

  vignette_path <- testthat::test_path("..", "..", "vignettes", "bayesian-example-data.Rmd")
  output_file <- rmarkdown::render(
    vignette_path,
    output_dir = tempdir(),
    quiet = TRUE,
    envir = new.env(parent = globalenv())
  )

  expect_true(file.exists(output_file))
})

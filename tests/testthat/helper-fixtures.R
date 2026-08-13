fixture_path <- function(name) {
  candidates <- c()

  testthat_candidate <- tryCatch(
    testthat::test_path("fixtures", name),
    error = function(e) NULL
  )

  if(!is.null(testthat_candidate)) {
    candidates <- c(candidates, testthat_candidate)
  }

  candidates <- c(
    candidates,
    file.path("tests", "testthat", "fixtures", name),
    file.path("testthat", "fixtures", name),
    file.path("fixtures", name)
  )

  matches <- unique(candidates[file.exists(candidates)])

  if(length(matches) == 0L) {
    stop("Fixture not found: ", name, call. = FALSE)
  }

  matches[[1]]
}

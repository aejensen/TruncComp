.validateSimulationInputs <- function(n, f0, f1, pi0, pi1, atom = 0) {
  if(!(length(n) == 1 && is.numeric(n) && is.finite(n) && n > 0 && n == as.integer(n))) {
    stop("n must be a single positive integer.")
  }

  if(!is.function(f0) || !is.function(f1)) {
    stop("f0 and f1 must be functions.")
  }

  probabilities <- list(pi0 = pi0, pi1 = pi1)
  for(name in names(probabilities)) {
    value <- probabilities[[name]]
    if(!(length(value) == 1 && is.numeric(value) && is.finite(value) && value >= 0 && value <= 1)) {
      stop(name, " must be a single number between 0 and 1.")
    }
  }

  if(!(length(atom) == 1 && is.numeric(atom) && is.finite(atom))) {
    stop("atom must be a single finite numeric value.")
  }
}

.drawObservedOutcome <- function(generator, n, label) {
  values <- generator(n)

  if(is.numeric(values) && length(values) == 1 && is.finite(values)) {
    if(n == 1L) {
      return(as.numeric(values))
    }

    values <- c(as.numeric(values),
                vapply(seq_len(n - 1L), function(i) generator(1), numeric(1)))
  }

  if(!(is.numeric(values) && length(values) == n && all(is.finite(values)))) {
    stop(label, " must return either a numeric vector of length n or a single finite numeric value when called with 1.")
  }

  as.numeric(values)
}

.simulateTruncatedGroup <- function(n, r, generator, probability, label, atom = 0) {
  alive <- stats::rbinom(n, 1, probability)
  observed <- .drawObservedOutcome(generator, n, label)

  if(any(alive == 1 & observed == atom)) {
    stop(label, " must not return the atom value for non-atom outcomes.")
  }

  data.frame(R = rep.int(r, n),
             A = alive,
             Y = ifelse(alive == 1, observed, atom))
}

#' Simulate two-group truncated outcome data
#'
#' Simulates two treatment groups with a binary non-atom indicator and a
#' continuous outcome conditional on being outside the atom. The substantive
#' atom is represented by the supplied `atom` value in the returned `Y` column,
#' so generated non-atom outcomes must differ from `atom`.
#'
#' @param n Number of observations per group.
#' @param f0 Function that generates non-atom outcomes for group `R = 0`. It may
#'   return a finite numeric vector of length `n` when called with `n`, or a
#'   single finite numeric value when called repeatedly with `1`. Any non-atom
#'   values used in the returned dataset must differ from `atom`.
#' @param f1 Function that generates non-atom outcomes for group `R = 1`. It may
#'   return a finite numeric vector of length `n` when called with `n`, or a
#'   single finite numeric value when called repeatedly with `1`. Any non-atom
#'   values used in the returned dataset must differ from `atom`.
#' @param pi0 Probability of being outside the atom in group `R = 0`.
#' @param pi1 Probability of being outside the atom in group `R = 1`.
#' @param atom Numeric atom value inserted into `Y` when `A = 0`.
#' @return A data frame with columns `R` (binary treatment indicator), `A`
#'   (binary non-atom indicator), and `Y` (coded endpoint with atom observations
#'   represented by `atom`).
#' @examples
#' f0 <- function(n) stats::rnorm(n, 3, 1)
#' f1 <- function(n) stats::rnorm(n, 3.5, 1)
#' simulate_truncated_data(n = 25, f0 = f0, f1 = f1, pi0 = 0.6, pi1 = 0.5)
#' @export
simulate_truncated_data <- function(n, f0, f1, pi0, pi1, atom = 0) {
  .validateSimulationInputs(n, f0, f1, pi0, pi1, atom = atom)
  n <- as.integer(n)

  d0 <- .simulateTruncatedGroup(n, 0, f0, pi0, "f0", atom = atom)
  d1 <- .simulateTruncatedGroup(n, 1, f1, pi1, "f1", atom = atom)

  rbind(d0, d1)
}

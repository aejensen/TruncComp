.validateSimulationInputs <- function(n, f0, f1, pi0, pi1, atom = 0) {
  if(!(length(n) == 1 && is.numeric(n) && is.finite(n) && n > 0 && n == as.integer(n))) {
    stop("n must be a single positive integer.")
  }

  if(!is.function(f0) || !is.function(f1)) {
    stop("f0 and f1 must be functions.")
  }

  for(probability in list(pi0 = pi0, pi1 = pi1)) {
    value <- probability[[1]]
    name <- names(probability)
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
    values <- vapply(seq_len(n), function(i) generator(1), numeric(1))
  }

  if(!(is.numeric(values) && length(values) == n && all(is.finite(values)))) {
    stop(label, " must return either a numeric vector of length n or a single finite numeric value when called with 1.")
  }

  as.numeric(values)
}

.simulateTruncatedGroup <- function(n, r, generator, probability, label, atom = 0) {
  alive <- stats::rbinom(n, 1, probability)
  observed <- .drawObservedOutcome(generator, n, label)

  data.frame(R = rep.int(r, n),
             A = alive,
             Y = ifelse(alive == 1, observed, atom))
}

#' Simulate two-group truncated outcome data
#'
#' Simulates two treatment groups with a binary observation indicator and a
#' continuous outcome among observed subjects. Unobserved outcomes are encoded as
#' the supplied `atom` value in the returned `Y` column.
#'
#' @param n Number of observations per group.
#' @param f0 Function that generates observed outcomes for group `R = 0`. It may
#'   return a finite numeric vector of length `n` when called with `n`, or a
#'   single finite numeric value when called repeatedly with `1`.
#' @param f1 Function that generates observed outcomes for group `R = 1`. It may
#'   return a finite numeric vector of length `n` when called with `n`, or a
#'   single finite numeric value when called repeatedly with `1`.
#' @param pi0 Probability that the outcome is observed in group `R = 0`.
#' @param pi1 Probability that the outcome is observed in group `R = 1`.
#' @param atom Numeric atom value inserted into `Y` when `A = 0`.
#' @return A data frame with columns `R` (binary treatment indicator), `A`
#'   (binary observation indicator), and `Y` (combined outcome with unobserved
#'   values replaced by `atom`).
#' @examples
#' f0 <- function(n) stats::rnorm(n, 3, 1)
#' f1 <- function(n) stats::rnorm(n, 3.5, 1)
#' simulateTruncatedData(n = 25, f0 = f0, f1 = f1, pi0 = 0.6, pi1 = 0.5)
#' @export
simulateTruncatedData <- function(n, f0, f1, pi0, pi1, atom = 0) {
  .validateSimulationInputs(n, f0, f1, pi0, pi1, atom = atom)
  n <- as.integer(n)

  d0 <- .simulateTruncatedGroup(n, 0, f0, pi0, "f0", atom = atom)
  d1 <- .simulateTruncatedGroup(n, 1, f1, pi1, "f1", atom = atom)

  rbind(d0, d1)
}

simTruncData <- function(n, mu0, mu1, pi0, pi1, sigma = 1, dist = "norm", df=4,
                         atom = 0) {
  generator <- function(mu) {
    function(n) {
      if(dist == "norm") {
        stats::rnorm(n, mu, sigma)
      } else if(dist == "t-sq") {
        stats::rt(n, df = df) + mu
      } else {
        stop("Invalid distribution")
      }
    }
  }

  simulateTruncatedData(n = n,
                        f0 = generator(mu0),
                        f1 = generator(mu1),
                        pi0 = pi0,
                        pi1 = pi1,
                        atom = atom)
}

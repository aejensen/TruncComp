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

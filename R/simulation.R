simTruncData <- function(n, mu0, mu1, pi0, pi1, sigma = 1, dist = "norm", df=4) {
  #pi0 probability of observing the outcome for Z = 0
  #pi1 probability of observing the outcome for Z = 1
  #CV = coefficient of variation for gamma outcome

  #Number of random treatment allocations
  #nTreatment <- rbinom(1, n, 0.5)
  #nControl <- n - nTreatment
  nTreatment <- n
  nControl <- n

  #Control group
  d0 <- t(sapply(1:nControl, function(i) {
    alive <- stats::rbinom(1, 1, pi0)
    if(dist == "norm") {
      y <- stats::rnorm(1, mu0, sigma)
    } else if(dist == "t-sq") {
      y <- stats::rt(1, df = df) + mu0
    } else {
      stop("Invalid distribution")
    }
    c(alive, y * as.numeric(alive == 1))
  }))
  d0 <- as.data.frame(cbind(0, d0))

  #Intervention group
  d1 <- t(sapply(1:nTreatment, function(i) {
    alive <- stats::rbinom(1, 1, pi1)
    if(dist == "norm") {
      y <- stats::rnorm(1, mu1, sigma)
    } else if(dist == "t-sq") {
      y <- stats::rt(1, df = df) + mu1
    } else {
      stop("Invalid distribution")
    }
    c(alive, y * as.numeric(alive == 1))
  }))
  d1 <- as.data.frame(cbind(1, d1))

  d <- rbind(d0, d1)
  colnames(d) <- c("R", "A", "Y")

  d
}

isDataOkay <- function(d) {
  dAlive <- subset(d, d$A == 1)

  yAlive1 <- dAlive$Y[dAlive$R == 0]
  yAlive2 <- dAlive$Y[dAlive$R == 1]

  if(length(yAlive1) < 2 | length(yAlive2) < 2) {
    return(FALSE)
  }

  TRUE
}

isValid <- function(truncCompObj) {
  truncCompObj$success
}

returnErrorData <- function(error, method, conf.level) {
  out <- list(muDelta = NULL,
              muDeltaCI = NULL,
              alphaDelta = NULL,
              alphaDeltaCI = NULL,
              W = NULL,
              p = NULL,
              method = method,
              conf.level = conf.level,
              success = FALSE,
              error = error,
              init = NULL,
              data = NULL)
  class(out) <- append(class(out), "CompositeOutcomeEstimation")
  out
}
